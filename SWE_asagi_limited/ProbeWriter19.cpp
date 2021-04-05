// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ProbeWriter19.h"
#include "muq_globals.h"

SWE::ProbeWriter19::ProbeWriter19(SWE::MySWESolver& solver) {
  // @TODO Please insert your code here.
}

SWE::ProbeWriter19::~ProbeWriter19() {
}

void SWE::ProbeWriter19::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void SWE::ProbeWriter19::finishPlotting() {
  // @TODO Please insert your code here.
}

void SWE::ProbeWriter19::mapQuantities(
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

  //std::vector<std::vector<double>> probe_point = {{ 545.735266126, 62.7164740303 },
  //						     { 1050.67821,   798.352124}};
  if(outputQuantities[4] > muq::solution[1+2*1]){
	  muq::solution[0+2*1] = timeStamp; 
	  muq::solution[1+2*1] = outputQuantities[4];
	  //std::cout <<"Probe" << 0 << " has time " << muq::solution[0+2*0]/60 << " and height " << muq::solution[1+2*0]*1000 << std::endl;
  }
}
