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

{ // Inverse UQ
  MPI_Init(&argc, &argv);
  auto comm = std::make_shared<parcer::Communicator>(MPI_COMM_WORLD);


  pt::ptree pt;

  pt.put("NumSamples", 100); // number of samples for single level
  pt.put("verbosity", 1); // show some output
  pt.put("BurnIn", 0);

  auto localFactory = std::make_shared<MyMIComponentFactory>(pt, comm);
  localFactory->SetComm(comm);
  auto componentFactory = std::make_shared<ParallelMIComponentFactory>(comm, comm, localFactory);

  if (comm->GetRank() == 0) {
  std::cout << std::endl << "*************** single chain reference" << std::endl << std::endl;

  {
    auto modelIndex = componentFactory->FinestIndex();
    modelIndex->SetValue(0,0);
    SLMCMC slmcmc (pt, componentFactory, modelIndex);
    std::shared_ptr<SampleCollection> samples = slmcmc.Run();

    std::cout << "SL mean Param: " << slmcmc.MeanParameter().transpose() << std::endl;
    std::cout << "SL mean QOI: " << slmcmc.MeanQOI().transpose() << std::endl;

    remove("slmcmc.hdf5");
    samples->WriteToFile("slmcmc.hdf5");
  }

  }
  componentFactory->finalize();

}
  MPI_Finalize();
  return 0;
}
