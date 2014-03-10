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

//const std::vector<HardwareQuantities> test_hardwareQuantities ={L1_CACHE_MISS,L2_CACHE_MISS,
//		L3_CACHE_MISS,BRANCH_MISPREDICTIONS,TOTAL_CYCLES,TOTAL_INSTRUCTIONS,FLPT_INSTRUCTIONS};

const std::vector<HardwareQuantities> test_hardwareQuantities ={L1_CACHE_MISS,
		TOTAL_CYCLES,TOTAL_INSTRUCTIONS,FLPT_INSTRUCTIONS};

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

	BOOST_TEST_MESSAGE("\n" << "GPTLHardwareCounter Message: \n"
							<< "test_hardwareQuantities = ");
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
}

BOOST_AUTO_TEST_CASE(check_getValues) {

	GPTLsetoption(PAPI_L1_TCM, 1);
//	GPTLsetoption(PAPI_L2_TCM, 1);
//	GPTLsetoption(PAPI_L3_TCM, 1);
//	GPTLsetoption(PAPI_BR_MSP, 1);
	GPTLsetoption(PAPI_TOT_CYC, 1);
	GPTLsetoption(PAPI_TOT_INS, 1);
	GPTLsetoption(PAPI_FP_INS, 1);
	GPTLsetoption(PAPI_FP_OPS, 1);
	GPTLinitialize();

	GPTLHardwareCounter tester("test",test_hardwareQuantities);
	GPTLstart("test");
	std::string testName = tester.getName();

	double a;
	for(unsigned i = 0; i < 4; i++)
		a = (double)i/10;

	BOOST_TEST_MESSAGE("\n" << "GPTLHardwareCounter Message: \n"
							<< "test_hardwareQuantities = ");

	GPTLstop("test");

	double val; // = 0.0;
	int gret = GPTLget_eventvalue("test", "PAPI_FP_OPS", -1, &val);
	if (gret < 0)
	{
		BOOST_TEST_MESSAGE("\n" << "GPTLget_eventvalue failed: " << gret);
	}
	else
		BOOST_TEST_MESSAGE("\n" << "PAPI_FP_OPS = " << val);

	for(unsigned i = 0; i < tester.getValues().size(); i++){

		BOOST_TEST_MESSAGE(" " << tester.getHardwareQuantities()[i]
		                       << " = " << tester.getValues()[i] << " ");
	}

//	GPTLpr(0);

	GPTLfinalize();
}

BOOST_AUTO_TEST_SUITE_END()
