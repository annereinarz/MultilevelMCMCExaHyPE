  export COMPILER_CFLAGS=-DEXAHYPE_LATE_TAKEOVER
  export EXAHYPE_LATE_TAKEOVER

   ../../../Toolkit/toolkit.sh ../SWE_MC_ADERDG.exahype2
   make clean
   make -j8 link_muq
