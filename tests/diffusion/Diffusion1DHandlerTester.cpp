#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <Diffusion1DHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the Diffusion1DHandler.
 */
BOOST_AUTO_TEST_SUITE(Diffusion1DHandler_testSuite)

/**
 * Method checking the initialization of the off-diagonal part of the Jacobian,
 * and the compute diffusion methods.
 */
BOOST_AUTO_TEST_CASE(checkDiffusion) {
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

	// Create the diffusion handler
	Diffusion1DHandler diffusionHandler;

	// Create ofill
	int mat[size*size];
	int *ofill = &mat[0];

	// Initialize it
	diffusionHandler.initializeOFill(network, ofill);

	// All the clusters diffuse except the last one
	BOOST_REQUIRE_EQUAL(ofill[0], 1);
	BOOST_REQUIRE_EQUAL(ofill[6], 1);
	BOOST_REQUIRE_EQUAL(ofill[12], 1);
	BOOST_REQUIRE_EQUAL(ofill[18], 1);

	// Check the total number of diffusing clusters
	BOOST_REQUIRE_EQUAL(diffusionHandler.getNumberOfDiffusing(), 4);

	// The size parameter
	double s = 1.0;

	// The arrays of concentration
	double concentration[3*size];
	double newConcentration[3*size];

	// Initialize their values
	for (int i = 0; i < 3*size; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Set the temperature to 1000 K to initialize the diffusion coefficients
	auto reactants = network->getAll();
	for (int i = 0; i < size; i++) {
		auto cluster = (PSICluster *) reactants->at(i);
		cluster->setTemperature(1000.0);
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	// Supposing the 3 grid points are laid-out as follow:
	// 0 | 1 | 2
	double *concOffset = conc + size;
	double *updatedConcOffset = updatedConc + size;

	// Fill the concVector with the pointer to the middle, left, and right grid points
	double **concVector = new double*[3];
	concVector[0] = concOffset; // middle
	concVector[1] = conc; // left
	concVector[2] = conc + 2 * size; // right

	// Compute the diffusion at this grid point
	diffusionHandler.computeDiffusion(network, s, concVector,
			updatedConcOffset);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[0], 3917913916907.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[1], 320772236824.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2], 157095680823.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3], 25265634.0, 0.01);
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4], 0.0, 0.01); // Does not diffuse

	// Initialize the indices and values to set in the Jacobian
	int nDiff = diffusionHandler.getNumberOfDiffusing();
	int indices[nDiff];
	double val[3*nDiff];
	// Get the pointer on them for the compute diffusion method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the diffusion a the grid point 1
	diffusionHandler.computePartialsForDiffusion(network, s, valPointer,
			indicesPointer);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(indices[0], 0);
	BOOST_REQUIRE_EQUAL(indices[1], 1);
	BOOST_REQUIRE_EQUAL(indices[2], 2);
	BOOST_REQUIRE_EQUAL(indices[3], 3);

	// Check some values
	BOOST_REQUIRE_CLOSE(val[1], 78358278338.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[4], 6415444736.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[5], 6415444736.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[6], -6283827232.0, 0.01);
	BOOST_REQUIRE_CLOSE(val[9], -1010625.0, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
