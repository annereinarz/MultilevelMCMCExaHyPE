#include "MUQ/SamplingAlgorithms/SLMCMC.h"
#include "MUQ/SamplingAlgorithms/GreedyMLMCMC.h"
#include "MUQ/SamplingAlgorithms/MIMCMC.h"
#include "MUQ/SamplingAlgorithms/ParallelMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/ParallelFixedSamplesMIMCMC.h"

#include "MUQ/SamplingAlgorithms/DummyKernel.h"

#include <boost/property_tree/ptree.hpp>

#include <initandsoon.h>

#include "calculateLikelihood.hh"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

#include <fstream>
#include <stdio.h>
#include <ctime>

#include "problem.h"

// NOTE: We only need our own static load balancing in case 1 proc per node plus TBB won't be enough

class MyStaticLoadBalancer : public StaticLoadBalancer {
    public:
      void setup(std::shared_ptr<ParallelizableMIComponentFactory> componentFactory, uint availableRanks) override {
        ranks_remaining = availableRanks;
        spdlog::info("Balancing load across {} ranks", availableRanks);
        auto indices = MultiIndexFactory::CreateFullTensor(componentFactory->FinestIndex()->GetVector());
        models_remaining = indices->Size();
      }

      int numCollectors(std::shared_ptr<MultiIndex> modelIndex) override {
        ranks_remaining--;
        return 1;
      }
      WorkerAssignment numWorkers(std::shared_ptr<MultiIndex> modelIndex) override {
        WorkerAssignment assignment;
        assignment.numWorkersPerGroup = 4;
        assignment.numGroups = (ranks_remaining / models_remaining) / 4;

        spdlog::debug("Of {}, assigning {} to model {}", ranks_remaining, assignment.numGroups * assignment.numWorkersPerGroup, *modelIndex);

        assert (assignment.numGroups * assignment.numWorkersPerGroup > 0);

        models_remaining--;
        ranks_remaining -= assignment.numGroups * assignment.numWorkersPerGroup;

        assert (ranks_remaining >= 0);

        return assignment;
      }
    private:
      uint ranks_remaining;
      uint models_remaining;

};

int main(int argc, char** argv){
  saved_argc = argc;
  saved_argv = argv;

  int initThreadProvidedThreadLevelSupport;
  //bool result = MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &initThreadProvidedThreadLevelSupport );

  count = 0;

  std::time_t result = std::time(nullptr);
  std::string timestamp = std::asctime(std::localtime(&result));
  auto tracer = std::make_shared<OTF2Tracer>("trace",timestamp);

  auto comm = std::make_shared<parcer::Communicator>(MPI_COMM_WORLD);

{ // Inverse UQ
  auto componentFactory = std::make_shared<MyMIComponentFactory>(comm);

  pt::ptree pt;

  pt.put("NumSamples", 100); // number of samples for single level
  pt.put("NumInitialSamples", 3); //ignore// number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); //ignore// estimator variance to be achieved by greedy algorithm
  pt.put("verbosity", 1); // show some output
  pt.put("MCMC.BurnIn", 10);
  pt.put("NumSamples_0", 1e2);
  pt.put("NumSamples_1", 5e1);
  pt.put("MLMCMC.Scheduling", true);
  pt.put("MLMCMC.Subsampling", 10);


  StaticLoadBalancingMIMCMC parallelMIMCMC (pt, componentFactory, std::make_shared<MyStaticLoadBalancer>(), comm, tracer);
  if (comm->GetRank() == 0) {
    parallelMIMCMC.Run();
    Eigen::VectorXd meanQOI = parallelMIMCMC.MeanQOI();
    std::cout << "mean QOI: " << meanQOI.transpose() << std::endl;
  }
  parallelMIMCMC.WriteToFile("parallelMIMCMC.h5");
  parallelMIMCMC.Finalize();

  comm->Barrier();
  if (comm->GetRank() == 0)
    remove("FullParallelMLMCMC.hdf5");
  comm->Barrier();

  parallelMIMCMC.WriteToFile("FullParallelMLMCMC.hdf5");


}

  tracer->write();

  muq::finalize();
  MPI_Finalize();
  return 0;
}
