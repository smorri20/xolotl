#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <W111TrapMutationHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <Options.h>
#include <DummyHandlerRegistry.h>
#include <DummyAdvectionHandler.h>
#include <mpi.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the W111TrapMutationHandler.
 */
BOOST_AUTO_TEST_SUITE(W111TrapMutationHandler_testSuite)

/**
 * Method checking the initialization and the compute modified trap-mutation methods.
 */
BOOST_AUTO_TEST_CASE(checkModifiedTrapMutation) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 0 0 10 6" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	argv = new char*[2];
	std::string parameterFile = "param.txt";
	argv[0] = new char[parameterFile.length() + 1];
	strcpy(argv[0], parameterFile.c_str());
	argv[1] = 0; // null-terminate the array
	opts.readParams(argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Create the network
	auto network = loader.generate(opts);
	// Get its size
	const int dof = network->getDOF();

	// Suppose we have a grid with 13 grip points and distance of
	// 0.1 nm between grid points
	int nGrid = 14;
	// Initialize the rates
	network->addGridPoints(nGrid);
	std::vector<double> grid;
	for (int l = 0; l < nGrid; l++) {
		grid.push_back((double) l * 0.1);
		network->setTemperature(1200.0, l);
	}
	// Set the surface position
	int surfacePos = 0;

	// Create the modified trap-mutation handler
	W111TrapMutationHandler trapMutationHandler;

	// Create the advection handlers needed to initialize the trap mutation handler
	std::vector<xolotlCore::IAdvectionHandler *> advectionHandlers;
	advectionHandlers.push_back(new DummyAdvectionHandler());

	// Initialize it
	trapMutationHandler.initialize(*network, grid);
	trapMutationHandler.initializeIndex1D(surfacePos, *network,
			advectionHandlers, grid);

	// The arrays of concentration
	double concentration[nGrid * dof];
	double newConcentration[nGrid * dof];

	// Initialize their values
	for (int i = 0; i < nGrid * dof; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the seventh grid point
	double *concOffset = conc + 6 * dof;
	double *updatedConcOffset = updatedConc + 6 * dof;

	// Putting the concentrations in the network so that the rate for
	// desorption is computed correctly
	network->updateConcentrationsFromArray(concOffset);

	// Compute the modified trap mutation at the seventh grid point
	trapMutationHandler.computeTrapMutation(*network, concOffset,
			updatedConcOffset, 6, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 5.0377e+26, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(updatedConcOffset[6], -5.0377e+26, 0.01); // He
	BOOST_REQUIRE_CLOSE(updatedConcOffset[15], 5.0377e+26, 0.01); // Create HeV

	// Get the offset for the twelfth grid point
	concOffset = conc + 11 * dof;
	updatedConcOffset = updatedConc + 11 * dof;

	// Putting the concentrations in the network so that the rate for
	// desorption is computed correctly
	network->updateConcentrationsFromArray(concOffset);

	// Compute the modified trap mutation at the twelfth grid point
	trapMutationHandler.computeTrapMutation(*network, concOffset,
			updatedConcOffset, 11, 0);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 3.1353e+22, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], 0.0, 0.01); // He2
	BOOST_REQUIRE_CLOSE(updatedConcOffset[16], 0.0, 0.01); // Doesn't create He2V
	BOOST_REQUIRE_CLOSE(updatedConcOffset[12], -7.8497e+21, 0.01); // He7
	BOOST_REQUIRE_CLOSE(updatedConcOffset[31], 7.8497e+21, 0.01); // He7V2

	// Initialize the indices and values to set in the Jacobian
	int nHelium = network->getAll(ReactantType::He).size();
	int indices[3 * nHelium];
	double val[3 * nHelium];
	// Get the pointer on them for the compute modified trap-mutation method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the modified trap-mutation at the grid point 11
	int nMutating = trapMutationHandler.computePartialsForTrapMutation(*network,
			valPointer, indicesPointer, 11, 0);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(nMutating, 5);
	BOOST_REQUIRE_EQUAL(indices[0], 8); // He3
	BOOST_REQUIRE_EQUAL(indices[1], 17); // He3V
	BOOST_REQUIRE_EQUAL(indices[2], 0); // I
	BOOST_REQUIRE_EQUAL(indices[3], 9); // He4
	BOOST_REQUIRE_EQUAL(indices[4], 18); // He4V
	BOOST_REQUIRE_EQUAL(indices[5], 0); // I

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 6.575931697e+14, 0.01);

	// Change the temperature of the network
	network->setTemperature(500.0);

	// Reinitialize the handler
	trapMutationHandler.initialize(*network, grid);
	// Update the bursting rate
	trapMutationHandler.updateTrapMutationRate(*network);

	// Compute the partial derivatives for the bursting a the grid point 11
	nMutating = trapMutationHandler.computePartialsForTrapMutation(*network,
			valPointer, indicesPointer, 11, 0);

	// Check values
	BOOST_REQUIRE_EQUAL(nMutating, 5);
	BOOST_REQUIRE_CLOSE(val[0], -5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], -5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 5.536237e+14, 0.01);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
