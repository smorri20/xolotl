#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include "W100FitFluxHandler.h"

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (W100FitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkgetIncidentFlux) {
	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l * 1.25);
	}

	// Create the W100 flux handler
    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(grid);

	// Create a time
	double currTime = 1.0;

	// Get the flux vector
	auto fluxVector = testFitFlux->getIncidentFluxVec(currTime);

	// Check a value
	BOOST_REQUIRE_CLOSE(fluxVector[1], 0.476819, 0.01);
}

BOOST_AUTO_TEST_CASE(checkHeFluence) {
	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l * 1.25);
	}

	// Create the W100 flux handler
    auto testFitFlux = make_shared<W100FitFluxHandler>();
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(grid);

	// Check that the fluence is 0.0 at the beginning
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), 0.0);

	// Increment the helium fluence
	testFitFlux->incrementHeFluence(1.0e-8);
	// Check that the fluence is not 0.0 anymore
	BOOST_REQUIRE_EQUAL(testFitFlux->getHeFluence(), 1.0e-8);
}

BOOST_AUTO_TEST_CASE(checkHeFlux) {
	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 5; l++) {
		grid.push_back((double) l * 1.25);
	}

	// Create the W100 flux handler
    auto testFitFlux = make_shared<W100FitFluxHandler>();

    // Set the factor to change the Helium flux
    testFitFlux->setHeFlux(2.5);
    // Initialize the flux handler
    testFitFlux->initializeFluxHandler(grid);

    BOOST_REQUIRE_EQUAL(testFitFlux->getHeFlux(), 2.5);
}


BOOST_AUTO_TEST_SUITE_END()
