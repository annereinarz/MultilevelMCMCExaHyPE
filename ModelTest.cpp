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
  saved_argc = argc;
  saved_argv = argv;
  muq::initParallelEnvironment(&argc,&argv);
  count = 0;


  auto comm = std::make_shared<parcer::Communicator>(MPI_COMM_WORLD);
  auto localFactory = std::make_shared<MyMIComponentFactory>(comm);
  localFactory->SetComm(comm);

  auto index = localFactory->FinestIndex();

  do {
    std::cout << "Testing model " << *index << std::endl;
    std::cout << "============================================" << std::endl;
    auto problem = localFactory->SamplingProblem(index);
    double logdens = problem->LogDensity(std::make_shared<SamplingState>(localFactory->StartingPoint(index)));
    std::cout << "LogDensity: " << logdens << std::endl;
    problem->QOI();
    index->SetValue(0, index->GetValue(0) - 1);
  } while (index->GetValue(0) > 0);


  muq::finalize();
  MPI_Finalize();
  return 0;
}
