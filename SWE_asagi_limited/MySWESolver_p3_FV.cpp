#include "MySWESolver_p3_FV.h"
#include "MySWESolver_p3_FV_Variables.h"

#include "kernels/KernelUtils.h"
#include "muq_globals.h"

using namespace kernels;

double grav_p3;
double epsilon_p3;
bool arrived_fv_p3;

tarch::logging::Log SWE::MySWESolver_p3_FV::_log( "SWE::MySWESolver_p3_FV" );

void SWE::MySWESolver_p3_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
    if (constants.isValueValidDouble( "grav_p3" )) {
        grav_p3 = constants.getValueAsDouble("grav_p3")/1000.0;
    }
    if (constants.isValueValidDouble( "epsilon_p3" )) {
        epsilon_p3 = constants.getValueAsDouble( "epsilon_p3" )/100.0;
    }
    arrived_fv_p3 = false;
}



void SWE::MySWESolver_p3_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
	// Dimensions             = 2
	// Number of variables    = 4 + #parameters

	if (tarch::la::equals(t,0.0)) {
		muq::initialData->getInitialData(x, Q);
	}
	else{
		if(Q[0] < epsilon_p3){
			Q[1] = 0;
			Q[2] = 0;      
		}
	}
}

void SWE::MySWESolver_p3_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters

  ReadOnlyVariables vars(Q);
  Variables eigs(lambda);

  const double c = std::sqrt(grav_p3 * vars.h());
  double u_n = Q[dIndex + 1] * vars.h() * std::sqrt(2)/std::sqrt(std::pow(vars.h(), 4) + std::pow(std::max(vars.h(), epsilon_p3), 4));

  eigs.h() = u_n + c;
  eigs.hu() = u_n - c;
  eigs.hv() = u_n;
  eigs.b() = 0.0;

  if(tarch::la::equals(u_n,0.0) && tarch::la::equals(c,0.0)){
    eigs.h() = std::sqrt(grav_p3 * epsilon_p3);
  }

}

void SWE::MySWESolver_p3_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int d,
    const double* const stateInside,
    double* const stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 4 + #parameters


  //OL Code
//  stateOutside[0] = stateInside[0];
//  stateOutside[1] = 0.0;
//  stateOutside[2] = 0.0;
//  stateOutside[3] = 0.0;


  //normal Code
  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];

  //Wall
  stateOutside[d + 1] = -stateInside[d + 1];
}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void SWE::MySWESolver_p3_FV::flux(const double* const Q,double** const F) {
  // Dimensions                        = 2
  // Number of variables + parameters  = 4 + 0
 
  ReadOnlyVariables vars(Q);

  double* f = F[0];
  double* g = F[1];

  double u_n = vars.hu() * vars.h() *std::sqrt(2)/std::sqrt(std::pow(vars.h(), 4) + std::pow(std::max(vars.h(), epsilon_p3), 4));
  double v_n = vars.hv() * vars.h() *std::sqrt(2)/std::sqrt(std::pow(vars.h(), 4) + std::pow(std::max(vars.h(), epsilon_p3), 4));

  f[0] = vars.h() * u_n;
  f[1] = vars.h() * u_n * u_n; // 0.5 * grav_p3 * vars.h() * vars.h();
  f[2] = vars.h() * u_n * v_n;
  f[3] = 0.0;

  g[0] = vars.h() * v_n;
  g[1] = vars.h() * u_n * v_n;
  g[2] = vars.h() * v_n * v_n; // 0.5 * grav_p3 * vars.h() * vars.h();
  g[3] = 0.0;
 
}
double SWE::MySWESolver_p3_FV::riemannSolver(double* const fL, double* const fR, const double* const qL, const double* const qR, const double* QL, const double* QR, const double* cellSize, int direction){
    double LL[NumberOfVariables] = {0.0};
    double LR[NumberOfVariables] = {0.0};

    eigenvalues(qL, direction, LL);
    eigenvalues(qR, direction, LR);

    double smax = 0.0;
    for (int i = 0; i < NumberOfVariables; i++) {
        const double abs_sL_i = std::abs(LL[i]);
        smax = std::max( abs_sL_i, smax );
    }
    for (int i = 0; i < NumberOfVariables; i++) {
        const double abs_sR_i = std::abs(LR[i]);
        smax = std::max( abs_sR_i, smax );
    }

    double FL2[DIMENSIONS][NumberOfVariables] = {0.0};
    double FR2[DIMENSIONS][NumberOfVariables] = {0.0};
    double* FL[DIMENSIONS]={FL2[0], FL2[1]};
    double* FR[DIMENSIONS]={FR2[0], FR2[1]};
    flux(qL, FL);
    flux(qR, FR);

    double flux[NumberOfVariables] = {0.0};

    flux[0] = 0.5 * (FL[direction][0] + FR[direction][0]) - 0.5 * smax * (qR[0] + qR[3] - qL[0] - qL[3]);
    flux[1] = 0.5 * (FL[direction][1] + FR[direction][1]) - 0.5 * smax * (qR[1] - qL[1]);
    flux[2] = 0.5 * (FL[direction][2] + FR[direction][2]) - 0.5 * smax * (qR[2] - qL[2]);
    flux[3] = 0.5 * (FL[direction][3] + FR[direction][3]);

	if (FL[direction][3] != 0.0 || FR[direction][3] != 0.0)
	std::cout << FL[direction][3] << ", " <<  FR[direction][3] << std::endl;

    double hRoe = 0.5*(qL[0] + qR[0]);

    double bm = std::max(qL[3], qR[3]);
    double Deta = std::max(qR[0]+qR[3] - bm, 0.0) - std::max(qL[0]+qL[3] - bm, 0.0);

    double djump[NumberOfVariables] = {0.0};

    djump[direction + 1] = 0.5*grav_p3*hRoe*Deta;


    flux[0] = 0.5 * (FL[direction][0] + FR[direction][0]) - 0.5*smax*Deta;
    for (int i = 0; i < NumberOfVariables; i++){
        fL[i] = flux[i] + djump[i];
        fR[i] = flux[i] - djump[i];
    }

    return smax;
}

void SWE::MySWESolver_p3_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  //do nothing: Should never be called
  std::cout << "Called fv ncp" << std::endl;
  std::terminate();
}

