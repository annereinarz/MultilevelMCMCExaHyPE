#include "MUQ/SamplingAlgorithms/SLMCMC.h"
#include "MUQ/SamplingAlgorithms/GreedyMLMCMC.h"
#include "MUQ/SamplingAlgorithms/MIMCMC.h"
#include "MUQ/SamplingAlgorithms/ParallelMIComponentFactory.h"
#include "MUQ/SamplingAlgorithms/ParallelFixedSamplesMIMCMC.h"

#include "MUQ/SamplingAlgorithms/DummyKernel.h"

#include <boost/property_tree/ptree.hpp>

#include "calculateLikelihood.hh"

namespace pt = boost::property_tree;
using namespace muq::Modeling;
using namespace muq::SamplingAlgorithms;
using namespace muq::Utilities;

#include <fstream>
#include <stdio.h>

#include "problem.h"

int main(int argc, char** argv){
  saved_argc = argc;
  saved_argv = argv;
  count = 0;

  MPI_Init(&argc, &argv);
  auto comm = std::make_shared<parcer::Communicator>(MPI_COMM_WORLD);

  pt::ptree pt;

  pt.put("NumSamples", 500); // number of samples for single level
  pt.put("verbosity", 1); // show some output
  pt.put("BurnIn", 100);

  auto localFactory = std::make_shared<MyMIComponentFactory>(pt, comm);

  /*auto comm = std::make_shared<parcer::Communicator>(MPI_COMM_WORLD);
  pt::ptree pt;
  auto localFactory = std::make_shared<MyMIComponentFactory>(pt, comm);*/
  localFactory->SetComm(comm);

  auto componentFactory = std::make_shared<ParallelMIComponentFactory>(comm, comm, localFactory);



  if (comm->GetRank() == 0) {
    Eigen::VectorXd testparam(2);
    testparam << 0.0,0.0;
    //testparam << 4.67169, 21.3935;
    for (int repeat = 0; repeat <= 1; repeat++) {
    auto index = localFactory->FinestIndex();
    while (true){
      std::cout << "Testing model " << *index << " for parameter " << std::endl << testparam << std::endl;
      std::cout << "============================================" << std::endl;
      auto problem = componentFactory->SamplingProblem(index);
      //double logdens = problem->LogDensity(std::make_shared<SamplingState>(localFactory->StartingPoint(index)));
      double logdens = problem->LogDensity(std::make_shared<SamplingState>(testparam));
      std::cout << "LogDensity: " << logdens << std::endl;
    
      std::cout << "Got QOI:" << std::endl;
      std::cout << problem->QOI()->state[0] << std::endl;
      if (index->GetValue(0) == 0)
        break;
      index->SetValue(0, index->GetValue(0) - 1);
    }
    }
  }
  componentFactory->finalize();
  MPI_Finalize();
  return 0;
}
