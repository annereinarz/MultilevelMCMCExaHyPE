// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ProbeWriter_p3.h"

SWE::ProbeWriter_p3::ProbeWriter_p3(SWE::MySWESolver_p3& solver) {
  // @TODO Please insert your code here.
}

SWE::ProbeWriter_p3::~ProbeWriter_p3() {
}

void SWE::ProbeWriter_p3::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void SWE::ProbeWriter_p3::finishPlotting() {
  // @TODO Please insert your code here.
}

void SWE::ProbeWriter_p3::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 5;
  for (int i=0; i<writtenUnknowns-1; i++){ 
    outputQuantities[i] = Q[i];
  }
  outputQuantities[4] = 0.0;
  if(Q[3] < 0.0)
      outputQuantities[4] = Q[3]+Q[0];
}
