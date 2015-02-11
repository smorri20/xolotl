#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "W221FitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W221FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l * 1.25);
	}

	// Create the W221 flux handler
    auto testFitFlux = make_shared<W221FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(grid);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto fluxVector = testFitFlux->getIncidentFluxVec(currTime);

	// Check a value
	BOOST_REQUIRE_CLOSE(fluxVector[1], 0.431739, 0.01);
}

BOOST_AUTO_TEST_SUITE_END()
