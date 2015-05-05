#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "W310FitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W310FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l * 1.25);
	}
	// Specify the surface position
	int surfacePos = 0;

	// Create the W310 flux handler
    auto testFitFlux = make_shared<W310FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(surfacePos, grid);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto testFluxVec = testFitFlux->getIncidentFluxVec(currTime, surfacePos);

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testFluxVec[1], 0.520850, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[2], 0.213396, 0.01);
	BOOST_REQUIRE_CLOSE(testFluxVec[3], 0.065754, 0.01);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
