#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <GBAdvectionHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the GBAdvectionHandler.
 */
BOOST_AUTO_TEST_SUITE(GBAdvectionHandler_testSuite)

/**
 * Method checking the initialization and the compute advection methods.
 */
BOOST_AUTO_TEST_CASE(checkAdvection) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();
	// Get its size
	const int size = network->getAll()->size();

	// Create the advection handler and initialize it
	GBAdvectionHandler advectionHandler;
	advectionHandler.initialize(network);

	// Check the total number of advecting clusters
	BOOST_REQUIRE_EQUAL(advectionHandler.getNumberOfAdvecting(), 6);

	// Set the size parameters
	double hx = 1.0;
	double hy = 0.5;
	double hz = 2.0;

	// The arrays of concentration
	double concentration[9*size];
	double newConcentration[9*size];

	// Initialize their values
	for (int i = 0; i < 9*size; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000K to initialize the diffusion coefficients
	auto reactants = network->getAll();
	for (int i = 0; i < size; i++) {
		auto cluster = (PSICluster *) reactants->at(i);
		cluster->setTemperature(1000.0);
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 9 grid points are laid-out as follow:
	// 6 | 7 | 8
	// 3 | 4 | 5
	// 0 | 1 | 2
	double *concOffset = conc + 4 * size;
	double *updatedConcOffset = updatedConc + 4 * size;

	// Fill the concVector with the pointer to the middle, left, right, bottom, and top grid points
	double **concVector = new double*[5];
	concVector[0] = concOffset; // middle
	concVector[1] = conc + 3 * size; // left
	concVector[2] = conc + 5 * size; // right
	concVector[3] = conc + 1 * size; // bottom
	concVector[4] = conc + 7 * size; // top

	// Set the grid position
	std::vector<double> gridPosition = { hx, hy, 0 };
	double h[3] = {hx, hy, hz};

	// Compute the advection at this grid point
	advectionHandler.computeAdvection(network, h, gridPosition,
			concVector, updatedConcOffset);

	// Check the new values of updatedConcOffset
//	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], -1.54376e+11, 0.01);
//	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], -1.49058e+11, 0.01);
//	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], -1.89362e+11, 0.01);
//	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], -6.65778e+10, 0.01);
//	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], -2.66416e+11, 0.01);
//	BOOST_REQUIRE_CLOSE(updatedConcOffset[5], -9.09579e+09, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[6], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[7], 0.0, 0.01); // Does not advect
	BOOST_REQUIRE_CLOSE(updatedConcOffset[8], 0.0, 0.01); // Does not advect

	// Initialize the rows, columns, and values to set in the Jacobian
	int nAdvec = advectionHandler.getNumberOfAdvecting();
	int indices[nAdvec];
	double val[2*nAdvec];
	// Get the pointer on them for the compute advection method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	advectionHandler.computePartialsForAdvection(network, h, valPointer,
			indicesPointer, gridPosition);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0);
	BOOST_REQUIRE_EQUAL(indices[1], 1);
	BOOST_REQUIRE_EQUAL(indices[2], 2);
	BOOST_REQUIRE_EQUAL(indices[3], 3);
	BOOST_REQUIRE_EQUAL(indices[4], 4);
	BOOST_REQUIRE_EQUAL(indices[5], 5);

	// Check values
//	BOOST_REQUIRE_CLOSE(val[0], -815207266.0, 0.01);
//	BOOST_REQUIRE_CLOSE(val[1], 50950454.0, 0.01);
//	BOOST_REQUIRE_CLOSE(val[2], -700039625.0, 0.01);
//	BOOST_REQUIRE_CLOSE(val[3], 43752477.0, 0.01);

	// Get the stencil
	auto stencil = advectionHandler.getStencilForAdvection(gridPosition);

	// Check the value of the stencil
	BOOST_REQUIRE_EQUAL(stencil[0], 0);
	BOOST_REQUIRE_EQUAL(stencil[1], -1); //y
	BOOST_REQUIRE_EQUAL(stencil[2], 0);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
