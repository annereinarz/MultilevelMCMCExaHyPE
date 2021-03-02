#include "MUQ/SamplingAlgorithms/SLMCMC.h"
#include "MUQ/SamplingAlgorithms/GreedyMLMCMC.h"
#include "MUQ/SamplingAlgorithms/MIMCMC.h"
#include "MUQ/SamplingAlgorithms/ParallelMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/ParallelFixedSamplesMIMCMC.h"

#include "MUQ/SamplingAlgorithms/DummyKernel.h"

#include <boost/property_tree/ptree.hpp>

#include <initialiseMUQ.h>

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
        assignment.numWorkersPerGroup = 1;
        assignment.numGroups = (ranks_remaining / models_remaining) / assignment.numWorkersPerGroup;

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
  spdlog::set_level(spdlog::level::debug);

  saved_argc = argc;
  saved_argv = argv;

  int initThreadProvidedThreadLevelSupport;
  //bool result = MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &initThreadProvidedThreadLevelSupport );
  muq::initParallelEnvironment(&argc,&argv);

  count = 0;

  auto comm = std::make_shared<parcer::Communicator>(MPI_COMM_WORLD);

  std::time_t result = std::time(nullptr);
  std::string timestamp = std::asctime(std::localtime(&result));
  comm->Bcast(timestamp, 0);
  auto tracer = std::make_shared<OTF2Tracer>("trace",timestamp);

{ // Inverse UQ
  auto componentFactory = std::make_shared<MyMIComponentFactory>(comm);

  pt::ptree pt;

  pt.put("verbosity", 1); // show some output
  pt.put("MCMC.BurnIn", 0);
  pt.put("NumSamples_0", 50);
  pt.put("NumSamples_1", 15);
  pt.put("MLMCMC.Scheduling",false);
  pt.put("MLMCMC.Subsampling", 0);


  StaticLoadBalancingMIMCMC parallelMIMCMC (pt, componentFactory, std::make_shared<MyStaticLoadBalancer>(), comm, tracer);
  if (comm->GetRank() == 0) {
    parallelMIMCMC.Run();
    Eigen::VectorXd meanQOI = parallelMIMCMC.MeanQOI();
    std::cout << "mean QOI: " << meanQOI.transpose() << std::endl;
  }

  if (comm->GetRank() == 0)
    remove("FullParallelMLMCMC.hdf5");

  parallelMIMCMC.WriteToFile("FullParallelMLMCMC.hdf5");
  parallelMIMCMC.Finalize();
}

  tracer->write();

  muq::finalize();
  MPI_Finalize();
  return 0;
}
