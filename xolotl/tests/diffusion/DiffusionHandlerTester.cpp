#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DiffusionHandler.h>

using namespace std;
using namespace xolotlCore;
/**
 * This suite is responsible for testing the DiffusionHandler.
 */
BOOST_AUTO_TEST_SUITE(DiffusionHandler_testSuite)

BOOST_AUTO_TEST_CASE(firstTest) {
	BOOST_TEST_MESSAGE("DiffusionHandlerTester Message: "
			<< "The tests are not implemented yet because the diffusion handler "
			<< "is going to change in the near future.");
}

BOOST_AUTO_TEST_SUITE_END()
