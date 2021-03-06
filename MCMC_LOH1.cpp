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

#include "problem_loh1.h"

int main(int argc, char** argv){
  int initThreadProvidedThreadLevelSupport;
  //bool result = MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &initThreadProvidedThreadLevelSupport );

  muq::init(argc,argv);
  count = 0;

{ // Inverse UQ
  auto localFactory = std::make_shared<MyMIComponentFactory>();

  pt::ptree pt;

  pt.put("NumSamples", 10); // number of samples for single level
  pt.put("NumInitialSamples", 3); //ignore// number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); //ignore// estimator variance to be achieved by greedy algorithm
  pt.put("verbosity", 1); // show some output
  pt.put("BurnIn", 0);

  auto comm = std::make_shared<parcer::Communicator>();
  localFactory->SetComm(comm);
  auto componentFactory = std::make_shared<ParallelMIComponentFactory>(comm, comm, localFactory);

  if (comm->GetRank() == 0) {
  std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  {
    SLMCMC slmcmc (pt, componentFactory);
    std::shared_ptr<SampleCollection> samples = slmcmc.Run();

    std::cout << "SL mean Param: " << slmcmc.MeanParameter().transpose() << std::endl;
    std::cout << "SL mean QOI: " << slmcmc.MeanQOI().transpose() << std::endl;

    remove("slmcmc.hdf5");
    samples->WriteToFile("slmcmc.hdf5");
  }

  }
  componentFactory->finalize();

}

  muq::finalize();
  MPI_Finalize();
  return 0;
}
