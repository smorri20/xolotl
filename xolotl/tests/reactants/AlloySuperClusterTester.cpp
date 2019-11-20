#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <AlloyCluster.h>
#include <AlloySuperCluster.h>
#include <AlloyClusterNetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <Constants.h>
#include <Options.h>
#include <fstream>
#include <iostream>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the AlloySuperCluster.
 */
BOOST_AUTO_TEST_SUITE (AlloySuperCluster_testSuite)

/**
 * This operation checks the ability of the AlloySuperCluster to describe
 * its connectivity to other clusters.
 */
BOOST_AUTO_TEST_CASE(checkConnectivity) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=200 0 0 5 4" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array
	// Initialize MPI for HDF5
	MPI_Init(&argc, &argv);

	// Read the options
	Options opts;
	opts.readParams(argc, argv);

	// Create the loader
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setMin(51);
	loader.setWidth(5);

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto& reactant = network->getAll(ReactantType::VoidSuper).begin()->second;

	// Check the type name
	BOOST_REQUIRE(ReactantType::VoidSuper == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for it
	int connectivityExpected[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the AlloySuperCluster get*Flux methods.
 */
BOOST_AUTO_TEST_CASE(checkFluxCalculations) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=200 0 0 5 4" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array

	// Read the options
	Options opts;
	opts.readParams(argc, argv);

	// Create the loader
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setMin(51);
	loader.setWidth(5);

	// Generate the network from the options
	auto network = loader.generate(opts);
	// Add a grid point for the rates
	network->addGridPoints(1);
	// Recompute Ids and network size and redefine the connectivities
	network->reinitializeConnectivities();

	// Get the super cluster
	auto& cluster = network->getAll(ReactantType::VoidSuper).begin()->second;

	// Get one that it combines with (V1)
	auto secondCluster = (AlloyCluster *) network->get(Species::V, 1);
	// Set the temperature and concentration
	network->setTemperature(1000.0, 0);
	cluster->setConcentration(0.5);
	secondCluster->setConcentration(0.5);

	// The flux can pretty much be anything except "not a number" (nan).
	double flux = cluster->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(115720.767655, flux, 0.000001);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the AlloySuperCluster get*PartialDerivatives methods.
 */
BOOST_AUTO_TEST_CASE(checkPartialDerivatives) {
	// Local Declarations
	// The vector of partial derivatives to compare with
	double knownPartials[] = {231441, 117091, 78270.5, 58598.9, 46665.2, 2.49036e+09,
		1.25993e+09, 8.42211e+08, 6.30539e+08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		6.24138e+07, 5.03899e+07, 4.22269e+07, 3.62979e+07, 3.17813e+07, 2.82163e+07,
		2.53239e+07, 2.29252e+07, 2.09001e+07, 1.91648e+07, 1.7659e+07, 1.63382e+07,
		1.51689e+07, 1.41252e+07, 1.3187e+07, 1.23381e+07, 1.15657e+07, 1.08593e+07,
		1.02103e+07, 9.61142e+06, 9.05672e+06, 8.54113e+06, 8.06034e+06, 7.61065e+06,
		7.1889e+06, 6.79233e+06, 6.41462e+06, 6.05753e+06, 5.71984e+06, 5.39987e+06,
		5.09613e+06, 4.80731e+06, 4.53222e+06, 4.26981e+06, 4.01913e+06, 3.77934e+06,
		3.54965e+06, 3.32938e+06, 3.11789e+06, 2.9146e+06, 2.71899e+06, 2.53058e+06,
		0, 0, 0, 0, 3.79335e-06, -0.0122539, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 3.79335e-06, 0.0122539, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=200 0 0 5 4" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array

	// Read the options
	Options opts;
	opts.readParams(argc, argv);

	// Create the loader
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setMin(51);
	loader.setWidth(5);

	// Generate the network from the options
	auto network = loader.generate(opts);
	// Add a grid point for the rates
	network->addGridPoints(1);

	// Set the temperature in the network
	double temperature = 1000.0;
	network->setTemperature(temperature, 0);
	// Redefine the connectivities
	network->reinitializeConnectivities();

	// Check the reaction connectivity of the super cluster
	auto& cluster = network->getAll(ReactantType::VoidSuper).begin()->second;
	// Set the cluster concentration
	cluster->setConcentration(0.5);

	// Get the vector of partial derivatives
	auto partials = cluster->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 224U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.001);
	}

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	return;
}

/**
 * This operation checks the reaction radius for AlloySuperCluster.
 */
BOOST_AUTO_TEST_CASE(checkReactionRadius) {
	// Create the parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=200 0 0 5 4" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array

	// Read the options
	Options opts;
	opts.readParams(argc, argv);

	// Create the loader
	AlloyClusterNetworkLoader loader = AlloyClusterNetworkLoader(
			std::make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Set grouping parameters
	loader.setMin(51);
	loader.setWidth(5);

	// Generate the network from the options
	auto network = loader.generate(opts);

	// Check the reaction radius of the super cluster
	auto& cluster = network->getAll(ReactantType::VoidSuper).begin()->second;
	BOOST_REQUIRE_CLOSE(0.775269, cluster->getReactionRadius(), 0.001);

	// Remove the created file
	std::string tempFile = "param.txt";
	std::remove(tempFile.c_str());

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
