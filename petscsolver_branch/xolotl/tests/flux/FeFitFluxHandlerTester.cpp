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

	// Create the Fe flux handler
    auto testFitFlux = make_shared<FeFitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(grid);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto fluxVector = testFitFlux->getIncidentFluxVec(currTime);

	// Check some values
	BOOST_REQUIRE_CLOSE(fluxVector[1], 0.116361, 0.1);
	BOOST_REQUIRE_CLOSE(fluxVector[10], 0.0467781, 0.1);

}

BOOST_AUTO_TEST_SUITE_END()
