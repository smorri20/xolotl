#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <GPTLHardwareCounter.h>
#include <iostream>
#include <math.h>
#include <papi.h>
#include <vector>
#include <string>

using namespace std;
using namespace xolotlPerf;

const std::vector<HardwareQuantities> test_hardwareQuantities =
	{L1_CACHE_MISS,L2_CACHE_MISS,L3_CACHE_MISS,BRANCH_MISPRED,TOTAL_CYCLES,TOTAL_INSTRUC,FLPT_INSTRUC};

/**
 * This suite is responsible for testing the HardwareCounter.
 */
BOOST_AUTO_TEST_SUITE (HardwareCounter_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	GPTLHardwareCounter tester("test",test_hardwareQuantities);

	//Output the version of PAPI that is being used
	BOOST_TEST_MESSAGE("\n" << "PAPI_VERSION = " << PAPI_VERSION_MAJOR(PAPI_VERSION) << "."
			  << PAPI_VERSION_MINOR(PAPI_VERSION) << "." << PAPI_VERSION_REVISION(PAPI_VERSION) << "\n");

	BOOST_REQUIRE_EQUAL("test", tester.getName());
}

BOOST_AUTO_TEST_CASE(check_getHardwareQuantities) {

	GPTLHardwareCounter tester("test",test_hardwareQuantities);

	BOOST_TEST_MESSAGE("\n" << "GPTLHardwareCounter Message: \n" << "test_hardwareQuantities = ");
	for(unsigned i = 0; i < test_hardwareQuantities.size(); i++){

		BOOST_TEST_MESSAGE(" " << tester.getHardwareQuantities()[i] << " ");
	}

#if READY
	for(unsigned i = 0; i < test_hardwareQuantities.size(); i++){

        // this attempts to compare a string against an enum value
		BOOST_REQUIRE_EQUAL(test_hardwareQuantities[i], tester.getHardwareQuantities()[i]);
	}
#endif // READY
}


BOOST_AUTO_TEST_SUITE_END()
