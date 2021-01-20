#ifndef MUQ_INIT
#define MUQ_INIT

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


namespace muq{

    int init(int argc, char** argv) ;
    std::vector<double> run_exahype(std::vector<double> param, int sr, int level=0);
    int finalize();
    bool setCommunicator(MPI_Comm comm);

}


#endif
