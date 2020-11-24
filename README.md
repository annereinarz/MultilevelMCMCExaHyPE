# MultilevelMCMCExaHyPE

This is the MUQ interface for ExaHyPE. To run a MCMC chain you need to install both ExaHyPE and MUQ first.

In ExaHyPE it is important to use the reinarz/muq branch and to switch Peano to the muq branch as well.
Also, in muq, currently use the linus/mimcmc branch until that is merged in master.

You must then run the toolkit (here it is assumed that ExaHyPE-Engine lies one folder above MultilevelMCMCExaHyPE, if not supply the correct path)

  ../ExaHyPE-Engine/Toolkit/toolkit.sh SWE_MC_ADERDG.exahype2

Now set the exports noted in export.sh to set options. Finally, call

  make -j4 link_muq

to build the muq integration.


Note: For MPI2 compatibility, also set 'export PROJECT_CFLAGS=-DMPI2' if needed!


In order to create the build directory type:

  cmake -DCMAKE_PREFIX_PATH=$HOME/[Path-to-MUQ]/build -DEXAHYPE_PATH=[Path-to-ExaHyPE]/ExaHyPE-Engine/ApplicationExamples/SWE/SWE_MC_ADERDG ..


The application can then be compiled using

  make


Running your binary requires the .exahype2 file to be passed:

  mpirun -np 4 ./MultilevelExaHyPE [Path-to-ExaHyPE]/ExaHyPE-Engine/ApplicationExamples/SWE/SWE_MC_ADC_ADERDG.exahype2



# Notes for SuperMUC

##MUQ deps:
since gaining access to all external repositories is a lot of work: use the download_dependencies script to download them and push MUQ dependencies as archives onto supermuc

##Modules:
module load tbb
module load python/3.6_intel
module load cmake/3.10

##ExaHyPE:
git clone git@gitlab.lrz.de:exahype/ExaHyPE-Engine.git
git checkout reinarz/muq

getting submodules via SSH:

./updateSubmodules.sh -s
./updateSubmodules.sh

export COMPILER_CFLAGS=-DEXAHYPE_LATE_TAKEOVER
export EXAHYPE_LATE_TAKEOVER

In ApplicationExamples/SWE/SWE_MC_ADERDG_l.exahype2 set:
        "architecture": "skx",

../../Toolkit/toolkit.sh SWE_MC_ADERDG_l.exahype2
make -j40 link_muq

#MUQ:

git clone git@bitbucket.org:mituq/muq2.git

Adapt paths to your local archives!
cmake -DCMAKE_CXX_FLAGS="-std=c++0x" -DMUQ_USE_OPENMP=OFF -DMUQ_USE_MPI=ON -DCMAKE_CXX_COMPILER=mpiicpc -DMPI_CXX_COMPILER=mpiicpc -DCMAKE_INSTALL_PREFIX=$PWD/install -DPARCER_EXTERNAL_SOURCE=$HOME/MUQ/muq_dependencies/parcer.zip -DHDF5_EXTERNAL_SOURCE=$HOME/MUQ/muq_dependencies/CMake-hdf5-1.8.19.tar.gz -DEIGEN_EXTERNAL_SOURCE=$HOME/MUQ/muq_dependencies/eigen-3.3.7.tar.bz2 -DNLOPT_EXTERNAL_SOURCE=$HOME/MUQ/muq_dependencies/nlopt-2.4.2.tar.gz -DBOOST_EXTERNAL_SOURCE=$HOME/MUQ/muq_dependencies/boost_1_63_0.tar.gz -DNANOFLANN_EXTERNAL_SOURCE=$HOME/MUQ/muq_dependencies/nanoflann.zip -DSTANMATH_EXTERNAL_SOURCE=$HOME/MUQ/muq_dependencies/stanmath.zip ..

make -j40 install

#MUQ+ExaHyPE:

git clone git@github.com:annereinarz/MultilevelMCMCExaHyPE.git

cmake -DCMAKE_PREFIX_PATH=$HOME/MUQ/muq2/build -DEXAHYPE_PATH=$HOME/MUQ/ExaHyPE-Engine/ApplicationExamples/SWE/SWE_MC_ADERDG ..
