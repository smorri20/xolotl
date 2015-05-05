#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "FeFitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the FeFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (FeFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 21; l++) {
		grid.push_back((double) l * 1.05);
	}
	// Specify the surface position
	int surfacePos = 0;

	// Create the Fe flux handler
    auto testFitFlux = make_shared<FeFitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime, surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 0.116361, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[5], 0.098120, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[15], 0.0, 0.01);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
