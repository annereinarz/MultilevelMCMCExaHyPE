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

#include "problem.h"

int main(int argc, char** argv){
  int initThreadProvidedThreadLevelSupport;
  //bool result = MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &initThreadProvidedThreadLevelSupport );

  muq::init(argc,argv);
  count = 0;

{ // Inverse UQ
  auto localFactory = std::make_shared<MyMIComponentFactory>();

  pt::ptree pt;

  pt.put("NumSamples", 100); // number of samples for single level
  pt.put("NumInitialSamples", 3); //ignore// number of initial samples for greedy MLMCMC
  pt.put("GreedyTargetVariance", 0.05); //ignore// estimator variance to be achieved by greedy algorithm
  pt.put("verbosity", 1); // show some output
  pt.put("BurnIn", 10);
  pt.put("NumSamples_0", 1e2);
  pt.put("NumSamples_1", 5e1);

  /*std::cout << std::endl << "*************** greedy multillevel chain" << std::endl << std::endl;

  GreedyMLMCMC greedymlmcmc (pt, componentFactory);
  greedymlmcmc.Run();
  std::cout << "mean QOI: " << greedymlmcmc.MeanQOI().transpose() << std::endl;
  greedymlmcmc.Draw(false);*/

  auto comm = std::make_shared<parcer::Communicator>();
  localFactory->SetComm(comm);
  auto componentFactory = std::make_shared<ParallelMIComponentFactory>(comm, comm, localFactory);

  if (comm->GetRank() == 0) {

    {
      MIMCMC mimcmc (pt, componentFactory);
      mimcmc.Run();
      std::cout << "ML mean QOI: " << mimcmc.MeanQOI().transpose() << std::endl;
    }


  //Write mean to file
  /*std::ofstream file("Input/parameters.csv");
  file << (slmcmc.MeanQOI()).format(CSVFormat);
  file.close();*/
  }
  componentFactory->finalize();

  //plot mean
  //std::vector<double> param = {slmcmc.MeanQOI()(0), slmcmc.MeanQOI()(1)};
  //muq::run_exahype(param);
  //calculateLikelihood();

}

  muq::finalize();
  MPI_Finalize();
  return 0;
}
