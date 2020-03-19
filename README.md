# MultilevelMCMCExaHyPE

This is the MUQ interface for ExaHyPE. To run a MCMC chain you need to install both ExaHyPE and MUQ first. 
In ExaHyPE it is important to use the muq branch and to switch Peano to the muq branch as well.


In order to create the build directory type:

cmake -DCMAKE_PREFIX_PATH=$HOME/exahype/muq2/build -DEXAHYPE_PATH=/home/hd/hd_hd/hd_bl386/exahype/ExaHyPE-Engine/Applica
tionExamples/SWE/SWE_MC_ADERDG ..


The application can then be compiled using 

make
