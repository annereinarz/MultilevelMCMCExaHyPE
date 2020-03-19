# MultilevelMCMCExaHyPE

This is the MUQ interface for ExaHyPE. To run a MCMC chain you need to install both ExaHyPE and MUQ first. 
In ExaHyPE it is important to use the muq branch and to switch Peano to the muq branch as well.


In order to create the build directory type:

cmake -DCMAKE_PREFIX_PATH=$HOME/[Path-to-MUQ]/build -DEXAHYPE_PATH=[Path-to-ExaHyPE]e/ExaHyPE-Engine/ApplicationExamples/SWE/SWE_MC_ADERDG ..


The application can then be compiled using 

make
