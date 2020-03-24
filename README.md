# MultilevelMCMCExaHyPE

This is the MUQ interface for ExaHyPE. To run a MCMC chain you need to install both ExaHyPE and MUQ first. 

In ExaHyPE it is important to use the reinarz/muq branch and to switch Peano to the muq branch as well.
You must then run the toolkit

  ../../../Toolkit/toolkit.sh ../SWE_MC_ADERDG.exahype2

build by calling

  make -j4

and then set the exports noted in export.sh to set options. Finally, call

  make link_muq

to build the muq integration.


Note: For MPI2 compatibility, also set 'export PROJECT_CFLAGS=-DMPI2' if needed!


In order to create the build directory type:

  cmake -DCMAKE_PREFIX_PATH=$HOME/[Path-to-MUQ]/build -DEXAHYPE_PATH=[Path-to-ExaHyPE]/ExaHyPE-Engine/ApplicationExamples/SWE/SWE_MC_ADERDG ..


The application can then be compiled using 

  make


Running your binary requires the .exahype2 file to be passed:

  mpirun -np 4 [Path-to-ExaHyPE]/ExaHyPE-Engine/ApplicationExamples/SWE/SWE_MC_ADC_ADERDG.exahype2
