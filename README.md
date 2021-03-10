# MultilevelMCMCExaHyPE

This is the MUQ interface for ExaHyPE. To run a MCMC chain you need to install both ExaHyPE and MUQ first.

In ExaHyPE it is important to use the reinarz/muq branch and to switch Peano to the muq branch as well.

  export COMPILER_CFLAGS=-DEXAHYPE_LATE_TAKEOVER
  export EXAHYPE_LATE_TAKEOVER

You must then run the toolkit (here it is assumed that ExaHyPE-Engine lies one folder above MultilevelMCMCExaHyPE, if not supply the correct path)

  ../ExaHyPE-Engine/Toolkit/toolkit.sh SWE_asagi_limited.exahype2

Now move to the SWE_asagi_limited directory and call

  make -j98 link_muq

to build the muq integration.


Note: For MPI2 compatibility, also set 'export PROJECT_CFLAGS=-DMPI2' if needed!


Once ExaHyPE is compiled, move to the build directory. In order to compile the programs in this repo, call

  cmake -DCMAKE_PREFIX_PATH=$HOME/[Path-to-MUQ]/build ..


The applications can then be compiled using

  make


Running your binary requires the .exahype2 file to be passed:

  mpirun -np 4 ./MLMCMC ../SWE_asagi_limited.exahype2


# Notes for SuperMUC

##MUQ deps:
since gaining access to all external repositories is a lot of work: use the download_dependencies script to download them and push MUQ dependencies as archives onto supermuc

##Modules:
module load tbb/2019
module load python/3.6_intel
module load cmake/3.14.4

##ExaHyPE:
git clone git@gitlab.lrz.de:exahype/ExaHyPE-Engine.git
git checkout reinarz/muq

getting submodules via SSH:

./updateSubmodules.sh -s
./updateSubmodules.sh

export COMPILER_CFLAGS=-DEXAHYPE_LATE_TAKEOVER
export EXAHYPE_LATE_TAKEOVER

In ApplicationExamples/SWE/SWE_asagi_limited.exahype2 set:
        "architecture": "skx",

../../Toolkit/toolkit.sh SWE_asagi_limited.exahype2
make -j40 link_muq

#MUQ:

git clone git@bitbucket.org:mituq/muq2.git

git checkout linus/mimcmc-cleanups

Adapt paths to your local archives!
cmake -DCMAKE_CXX_FLAGS="-std=c++0x" -DMUQ_USE_OPENMP=OFF -DMUQ_USE_MPI=ON -DCMAKE_CXX_COMPILER=mpiicpc -DMPI_CXX_COMPILER=mpiicpc -DCMAKE_INSTALL_PREFIX=$PWD/install -DPARCER_EXTERNAL_SOURCE=/hppfs/work/pr83no/ge68wax4/MUQ/muq_dependencies/parcer.zip -DEIGEN_EXTERNAL_SOURCE=/hppfs/work/pr83no/ge68wax4/MUQ/muq_dependencies/eigen-3.3.7.tar.bz2 -DNLOPT_EXTERNAL_SOURCE=/hppfs/work/pr83no/ge68wax4/MUQ/muq_dependencies/nlopt-2.4.2.tar.gz -DBOOST_EXTERNAL_SOURCE=/hppfs/work/pr83no/ge68wax4/MUQ/muq_dependencies/boost_1_63_0.tar.gz -DNANOFLANN_EXTERNAL_SOURCE=/hppfs/work/pr83no/ge68wax4/MUQ/muq_dependencies/nanoflann.zip -DSTANMATH_EXTERNAL_SOURCE=/hppfs/work/pr83no/ge68wax4/MUQ/muq_dependencies/stanmath.zip ..

make -j40 install

#MUQ+ExaHyPE:

git clone git@github.com:annereinarz/MultilevelMCMCExaHyPE.git

CC=mpicc CXX=mpicxx cmake -DCMAKE_PREFIX_PATH=/dss/dsshome1/0C/ge68wax4/MUQ/muq2/build ..

