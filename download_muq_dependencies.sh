#!/bin/bash

mkdir muq_dependencies
cd muq_dependencies

wget https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.bz2
wget http://downloads.sourceforge.net/project/boost/boost/1.63.0/boost_1_63_0.tar.gz
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.19/src/CMake-hdf5-1.8.19.tar.gz
wget https://github.com/jlblancoc/nanoflann/archive/master.zip
wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz
wget https://bitbucket.org/mituq/parcer/get/master.zip
wget https://github.com/stan-dev/math/archive/release/v2.18.0.zip

cd ..
tar -czvf muq_dependencies.tar.gz muq_dependencies

