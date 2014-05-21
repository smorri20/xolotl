#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <HDF5Utils.h>
#include <PSIClusterReactionNetwork.h>
#include <DummyHandlerRegistry.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <memory>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the HDF5NetworkLoader.
 */
BOOST_AUTO_TEST_SUITE(HDF5NetworkLoader_testSuite)

/**
 * Method checking the loading of the network from the HDF5 file.
 */
BOOST_AUTO_TEST_CASE(checkLoad) {
	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/reactants/testfiles/tungsten_small.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = loader.load();

	// Get the size of the network
	int networkSize = network->size();
	// Check the value
	BOOST_REQUIRE_EQUAL(networkSize, 100);
}

BOOST_AUTO_TEST_SUITE_END()
