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

#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/logging/LogFilterFileReader.h"
#include "tarch/parallel/Node.h"

#include "peano/peano.h"
#include "peano/parallel/JoinDataBufferPool.h"

#include "exahype/main.h"
#include "exahype/parser/Parser.h"
#include "exahype/Vertex.h"
#include "exahype/runners/Runner.h"

#include "muq_globals.h"
#include "kernels/KernelCalls.h"

#include "kernels/GaussLegendreBasis.h"

#include <vector>
#include <cstdlib> // getenv, exit
#include <iostream>
#include <cstdio>

#include <exahype/main.h>

tarch::logging::Log _log("exahype");


namespace muq{

	exahype::parser::Parser parser;
	std::vector<std::string> cmdlineargs;
	bool programExitCode;
	
	//initialise global variables
	InitialData* initialData;
	InitialData* initialDataFV;
	
	double grav = 9.81;
	double epsilon_DG = 1e-2;
        double epsilon = 1e-2;
	int level = 0;
	int subcommunicator_rank = 0;
	int sample_number = 0;
	std::vector<double> param = {0.6,0.6};
	std::vector<double> solution = {};

	int finalize(){
		peano::shutdownParallelEnvironment();
		peano::shutdownSharedMemoryEnvironment();
		peano::releaseCachedData();
		return programExitCode;
	}

	bool initParallelEnvironment(int* argc, char*** argv){
		//
		//   Setup environment
		//   =================
		//

		int parallelSetup = peano::initParallelEnvironment(argc, argv);
		if (parallelSetup != 0) {
#ifdef Parallel
			// Please do not use the logging if MPI doesn't work properly.
			std::cerr << "mpi initialisation wasn't successful. Application shut down"
				<< std::endl;
#else
			_log.error("main()",
					"mpi initialisation wasn't successful. Application shut down");
#endif
			return false;
		}
		int sharedMemorySetup = peano::initSharedMemoryEnvironment();
		if (sharedMemorySetup != 0) {
			logError("main()",
					"shared memory initialisation wasn't successful. Application shut "
					"down");
			return false;
		}
		return true;
	}

	bool setCommunicator(MPI_Comm communicator){
		tarch::parallel::Node::getInstance().setCommunicator(communicator);
		return true;
	}

	int init(int argc, char** argv) {
		//
		//   Parse config file
		// =====================
		//
		std::string progname = argv[0];

		if (argc < 2) {
			logError("main()", "Usage: " << progname << " --help");
			return -1;
		}

		// cmdlineargs contains all argv expect the progname.
		std::vector<std::string> cmdlineargs_(argv + 1, argv + argc);
		cmdlineargs = cmdlineargs_;
		std::string firstarg = cmdlineargs[0];

		bool showHelp    = firstarg == "-h" || firstarg == "--help";
		bool showVersion = firstarg == "-v" || firstarg == "--version";
		bool runTests    = firstarg == "-t" || firstarg == "--tests";
		bool runPingPong = firstarg == "-p" || firstarg == "--pingpong";
		bool showCompiledSpecfile = firstarg == "--show-specfile";
		bool runCompiledSpecfile  = firstarg == "--built-in-specfile";

		//
		//   Early standalone options
		//   ========================
		//

		if(showHelp) {
			//help(progname);
			return EXIT_SUCCESS;
		}

		if(showVersion) {
			//std::cout << version(progname);
			return EXIT_SUCCESS;
		}

		if(showCompiledSpecfile) {
			// Unfortunately, we cannot avoid here to get the output dirtied by the
			// tarch::parallel::Node<static>::reserveFreeTag() log outputs.
			// The only alternative to get the clean specfile would be to dump it to
			// a file.

			// if this line does not compile for you, rebuild and rerun the toolkit.
			std::cout << std::string(kernels::compiledSpecfile());
			return EXIT_SUCCESS;
		}

		//
		//   Setup environment
		//   =================
		//
		peano::fillLookupTables();

		if (runTests) {
			//
			//   Run tests
			// =============
			// Our unit tests do cover the generic ADER-DG kernels. The generic kernels do
			// parallelise. As a consequence, they connect to the autotuning feature.
			// Autotuning however is not set up yet, so this will fail. We therefore
			// disable the unit tests in shared memory mode.
			//

			tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().run();
			int testExitCode = tarch::tests::TestCaseRegistry::getInstance()
				.getTestCaseCollection()
				.getNumberOfErrors();

			if (testExitCode != 0) {
				logError("main()", "unit tests failed. Quit.");
				return -2;
			}
			else {
				logInfo("main()", "all unit tests completed successfully.");
				return EXIT_SUCCESS;
			}
		}

		//
		//   Parse specification file
		// =====================================
		//
		std::stringstream specfile;
		std::string specFileName;
		if(runCompiledSpecfile) {
			specFileName = "builtin";
			specfile.str(std::string(kernels::compiledSpecfile()));
		} else {
			specFileName = firstarg;
			specfile.str(kernels::readSpecificationFileToJSON(specFileName));
		}
		parser.readFile(specfile, specFileName);

		if (!parser.isValid()) {
			logError("main()", "invalid config file. Quit");
			return -2;
		}

		//
		//   Init solver registries
		// =====================================
		//

		//
		//   Configure the logging
		// =========================
		//
		tarch::logging::CommandLineLogger::getInstance().clearFilterList();
#if defined(Parallel) || defined(PerformanceAnalysis)
		tarch::logging::CommandLineLogger::getInstance().setLogFormat(
				" ",    // columnSeparator
				true,   // logTimeStamp
				false,  // logTimeStampHumanReadable
				true,   // logMachineName
				true,   // logMessageType
				true,   // logTrace
				parser.getLogFileName() );
#elif defined(Asserts) || defined(Debug)
		tarch::logging::CommandLineLogger::getInstance().setLogFormat(
				" ",    // columnSeparator
				true,   // logTimeStamp
				false,  // logTimeStampHumanReadable
				false,  // logMachineName
				true,   // logMessageType
				true,   // logTrace
				parser.getLogFileName() );
#else
		tarch::logging::CommandLineLogger::getInstance().setLogFormat(
				" ",    // columnSeparator
				true,   // logTimeStamp
				false,  // logTimeStampHumanReadable
				false,  // logMachineName
				true,   // logMessageType
				false,   // logTrace
				parser.getLogFileName() );
#endif

		tarch::logging::CommandLineLogger::getInstance().clearFilterList();
		tarch::logging::LogFilterFileReader::parsePlainTextFile( "exahype.log-filter" );

		muq::param.resize(2);
		muq::solution.resize(4);
		return 0;
	}

	std::vector<double> run_exahype(std::vector<double> param_, int globalcommunicator_rank_, int level_=0){
		muq::subcommunicator_rank = globalcommunicator_rank_;
		param=param_;  //store parameters
		level = level_;
		for(int i=0; i<muq::solution.size();i++)
			muq::solution[i] = -1234.5;
		if(level == 0){
			initialData = new InitialData(15,"data_gmt.yaml");
			initialDataFV = new InitialData(15,"data_gmt.yaml");
		}
		else{
			initialData = new InitialData(14,"data_gmt.yaml");
			initialDataFV = new InitialData(14,"data_gmt.yaml");

		}
		kernels::registerSolvers(parser,0);
		exahype::runners::Runner runner(parser, cmdlineargs);
		int programExitCode = runner.run();
		if (programExitCode == 0) {
#ifdef Parallel
			if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
				logInfo("main()", "Peano terminates successfully");
			}
#else
			logInfo("main()", "Peano terminates successfully");
#endif
		} else {
			logInfo("main()", "quit with error code " << programExitCode);
		}
		delete initialDataFV;
		delete initialData;
		kernels::finalise(0);
		sample_number++;
		return solution;
	}

}
