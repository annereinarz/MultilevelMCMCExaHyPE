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

#include "problem.h"

int main(int argc, char** argv){
  saved_argv = argv;
  saved_argc = argc;
  muq::initParallelEnvironment(&argc,&argv);
  count = 0;

  auto comm = std::make_shared<parcer::Communicator>(MPI_COMM_WORLD);

{ // Inverse UQ
  auto localFactory = std::make_shared<MyMIComponentFactory>(comm);

  pt::ptree pt;

  pt.put("NumSamples", 100); // number of samples for single level
  pt.put("NumInitialSamples", 3); //ignore// number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); //ignore// estimator variance to be achieved by greedy algorithm
  pt.put("verbosity", 1); // show some output
  pt.put("BurnIn", 10);
  pt.put("NumSamples_0", 1e2);
  pt.put("NumSamples_1", 5e1);

  localFactory->SetComm(comm);
  auto componentFactory = std::make_shared<ParallelMIComponentFactory>(comm, comm, localFactory);

  if (comm->GetRank() == 0) {

    {
      MIMCMC mimcmc (pt, componentFactory);
      mimcmc.Run();
      remove("mlmcmc.hdf5");
      mimcmc.WriteToFile("mlmcmc.hdf5");
      std::cout << "ML mean QOI: " << mimcmc.MeanQOI().transpose() << std::endl;
    }

  }
  componentFactory->finalize();

}

  muq::finalize();
  MPI_Finalize();
  return 0;
}
