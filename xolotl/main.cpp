/**
 * Main.c, currently only able to load clusters
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <Reactant.h>
#include <PSIClusterNetworkLoader.h>
#include <PetscSolver.h>
#include <mpi.h>
#include "xolotlCore/io/MPIUtils.h"
#include "xolotlPerf/HandlerRegistryFactory.h"
#include "xolotlCore/commandline/XolotlOptions.h"


using namespace std;
using std::shared_ptr;

//! This operation prints the start message
void printStartMessage() {
	cout << "Starting Xolotl Plasma-Surface Interactions Simulator" << endl;
	// TODO! Print copyright message
	// TODO! Print date and time
}


//! Main program
int main(int argc, char **argv) {

	// Local Declarations
	shared_ptr<std::istream> networkStream;
	std::shared_ptr<PSIClusterNetworkLoader> networkLoader;
	int rank;

	// Check the command line arguments.
    // Skip the executable name before parsing.
    argc -= 1;  // one for the executable name
    argv += 1;  // one for the executable name
    XolotlOptions xopts;
    int nOptsUsed = xopts.parseCommandLine( argc, argv );
    if( !xopts.shouldRun() )
    {
        return xopts.getExitCode();
    }
    argc -= nOptsUsed;
    argv += nOptsUsed;

	// Extract the argument for the file name
    std::string networkFilename = xopts.getNetworkFilename();
    assert( !networkFilename.empty() );

	try {
        // Set up our performance data infrastructure
        if( !xolotlPerf::initialize( xopts.useStandardHandlers() ) )
        {
            std::cerr << "Unable to initialize requested performance data infrastructure.  Aborting" << std::endl;
            return EXIT_FAILURE;
        }
        std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry = xolotlPerf::getHandlerRegistry();
        std::shared_ptr<xolotlPerf::ITimer> totalTimer = handlerRegistry->getTimer( "total" );
        totalTimer->start();

		// Setup the solver
        std::shared_ptr<xolotlPerf::ITimer> solverInitTimer = handlerRegistry->getTimer( "initSolver" );
        solverInitTimer->start();
		xolotlSolver::PetscSolver solver;
		solver.setCommandLineOptions(argc, argv);
		solver.initialize();
        solverInitTimer->stop();

        std::shared_ptr<xolotlPerf::ITimer> networkLoadTimer =  handlerRegistry->getTimer( "loadNetwork" );
        networkLoadTimer->start();

		// Get the MPI rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		// Setup the master
		if (rank == 0) {
			// Say hello
			printStartMessage();
			// Set the input stream on the master
			networkStream = make_shared<std::ifstream>(networkFilename);
		}

		// Broadcast the stream to all worker tasks
		networkLoader = std::make_shared<
					PSIClusterNetworkLoader>();
		networkStream = xolotlCore::MPIUtils::broadcastStream(networkStream, 0,
				MPI_COMM_WORLD );

		// Create a network loader and set the stream on every MPI task
		networkLoader->setInputstream(networkStream);
		// Give the network loader to PETSc as input
		solver.setNetworkLoader(networkLoader);
        networkLoadTimer->stop();

		// Launch the PetscSolver
        std::shared_ptr<xolotlPerf::ITimer> solverTimer = handlerRegistry->getTimer( "solve" );
        solverTimer->start();
		solver.solve();
        solverTimer->stop();

		solver.finalize();
        totalTimer->stop();

        // Report the performance data about the run we just completed
        // TODO implement our own performance data output mechanism
        // rather than relying on GPTL's output.

	} catch (std::string & error) {
		std::cout << error << std::endl;
		std::cout << "Aborting." << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
