// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ConservedWriter.h"

SWE::ConservedWriter::ConservedWriter(SWE::MySWESolver& solver) {
  // @TODO Please insert your code here.
}

SWE::ConservedWriter::~ConservedWriter() {
}

void SWE::ConservedWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void SWE::ConservedWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void SWE::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 5;
  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i] = Q[i];
  }
}