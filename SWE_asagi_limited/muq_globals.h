#ifndef MUQ_EXTERN
#define MUQ_EXTERN

/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
#include <vector>
#include "InitialData.h"

namespace muq{

    extern int subcommunicator_rank;
    extern double mesh_width;
    extern std::vector<double> param;
    extern std::vector<double> solution;
    extern InitialData* initialData_nobath;
    extern InitialData* initialData;

    /*void gatherProbes(int numProbes){
	    //Ensure that the solutions are distributed to all nodes
	    MPI_Comm comm = tarch::parallel::Node::getInstance().getCommunicator();
	    int numRanks = tarch::parallel::Node::getInstance().getNumberOfNodes();
	    for(int idx = 0; idx < numProbes; idx++){
		    double solutions[numRanks];
		    MPI_Allgather(&muq::solution[idx],1,MPI_DOUBLE,&solutions,1,MPI_DOUBLE,comm);
		    double compare = solutions[0];
		    for(int rank = 1; rank < numRanks; rank++) //find correct solution for each probe point
			    if(std::abs(solutions[rank]) > std::abs(compare))
				    muq::solution[idx] = solutions[rank];
	    }	
	    MPI_Barrier(comm);
    }*/



}


#endif
