#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <HDF5Utils.h>
#include <PSIClusterReactionNetwork.h>
#include <DummyHandlerRegistry.h>
#include <PSIClusterNetworkLoader.h>
#include <XolotlConfig.h>
#include <mpi.h>
#include <memory>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the HDF5Utils
 */
BOOST_AUTO_TEST_SUITE(HDF5Utils_testSuite)

/**
 * Method checking the writing and reading of the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkOI) {
	// Initialize MPI for HDF5
	int argc;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	PSIClusterNetworkLoader loader =
			PSIClusterNetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/tungsten.txt");
	string filename = sourceDir + pathToFile;
	// Create the network stream
	shared_ptr<istream> networkStream;
	networkStream = make_shared<ifstream>(filename);
	// Read the buffer of the stream
	auto bufferSS = make_shared<stringstream>();
	(*bufferSS) << networkStream->rdbuf();
	// Give the network stream to the network loader
	loader.setInputstream(bufferSS);

	// Load the network
	auto network = loader.load();

	// Get the size of the network
	int networkSize = network->size();
	// Set the time step number
	int timeStep = 0;
	// Initialize the HDF5 file
	HDF5Utils::initializeFile(timeStep, networkSize, 1);

	// Set the physical dimension of the grid and the refinement
	int dimension = 5;
	int refinement = 0;
	// Set the time information
	double currentTime = 0.0001;
	double currentTimeStep = 0.000001;
	// Write the header in the HDF5 file
	HDF5Utils::fillHeader(dimension, refinement, currentTime, currentTimeStep);

	// Write the network in the HDF5 file
	HDF5Utils::fillNetwork(network);

	// Create an array of concentration for one grid point
	double concentrations[networkSize];
	// Fill it
	for (int i = 0; i < networkSize; i++) {
		concentrations[i] = (double) i / 5.0;
	}
	double * conc = &concentrations[0];
	// Set the position at this grid point
	int gridPoint = 0;
	double position = 1.5;
	// Write the concentrations in the HDF5 file
	HDF5Utils::fillConcentrations(conc, gridPoint, position);

	// Finalize the HDF5 file
	xolotlCore::HDF5Utils::finalizeFile();

	// Read the header of the written file
	int dim = 0;
	double t = 0.0, dt = 0.0;
	HDF5Utils::readHeader("xolotlStop_0.h5", dim, t, dt);
	// Check the obtained values
	BOOST_REQUIRE_EQUAL(dim, dimension);
	BOOST_REQUIRE_EQUAL(t, currentTime);
	BOOST_REQUIRE_EQUAL(dt, currentTimeStep);

	// Read the network of the written file
	auto networkVector = HDF5Utils::readNetwork("xolotlStop_0.h5");
	// Get all the reactants
	auto reactants = network->getAll();
	// Check the network vector
	for (int i = 0; i < networkSize; i++) {
		// Get the i-th reactant in the network
		shared_ptr<PSICluster> reactant = static_pointer_cast<PSICluster>(reactants->at(i));
		int id = reactant->getId() - 1;
		// Get the corresponding line from the HDF5 file
		auto line = networkVector.at(id);

		// Check the composition
		auto composition = reactant->getComposition();
		BOOST_REQUIRE_EQUAL((int) line[0], composition["He"]);
		BOOST_REQUIRE_EQUAL((int) line[1], composition["V"]);
		BOOST_REQUIRE_EQUAL((int) line[2], composition["I"]);

		// Check the binding energies
		auto bindingEnergies = reactant->getBindingEnergies();
		BOOST_REQUIRE_EQUAL(line[3], bindingEnergies.at(0)); // Helium binding energy
		BOOST_REQUIRE_EQUAL(line[4], bindingEnergies.at(1)); // Vacancy binding energy
		BOOST_REQUIRE_EQUAL(line[5], bindingEnergies.at(2)); // Interstitial binding energy

		// Check the migration energy
		double migrationEnergy = reactant->getMigrationEnergy();
		BOOST_REQUIRE_EQUAL(line[6], migrationEnergy);

		// Check the diffusion factor
		double diffusionFactor = reactant->getDiffusionFactor();
		BOOST_REQUIRE_EQUAL(line[7], diffusionFactor);
	}

	// Read the concentrations at the given grid point
	double newConcentrations[networkSize];
	double * newConc = &newConcentrations[0];
	HDF5Utils::readGridPoint("xolotlStop_0.h5", networkSize, gridPoint, newConc);
	// Check them
	for (int i = 0; i < networkSize; i++) {
		BOOST_REQUIRE_EQUAL(newConcentrations[i], concentrations[i]);
	}

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
