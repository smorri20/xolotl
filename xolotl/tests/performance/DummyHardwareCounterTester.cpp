#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
//#include <HardwareCounter.h>
//#include <HardwareQuantities.h>
#include <DummyHardwareCounter.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

//std::vector<HardwareQuantities> test_hquan = {L1_CACHE_MISS,L2_CACHE_MISS,L3_CACHE_MISS,BRANCH_MISPRED,TOTAL_CYCLES,TOTAL_INSTRUC,FLPT_INSTRUC};

/**
 * This suite is responsible for testing the DummyHardwareCounter.
 */
BOOST_AUTO_TEST_SUITE (DummyHardwareCounter_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

//	DummyHardwareCounter tester("test",test_hquan);
	DummyHardwareCounter tester("test");

	BOOST_REQUIRE_EQUAL("", tester.getName());
}

//BOOST_AUTO_TEST_CASE(checkInitialValue) {
//
//	DummyHardwareCounter tester("test",test_hquan);;
//
//	BOOST_REQUIRE_EQUAL(0, tester.getValues());
//
//}

//BOOST_AUTO_TEST_CASE(checkCounting) {
//
//	DummyHardwareCounter tester("test",test_hquan);;
//
//	long count = 3;
//
//	for(int i = 0; i < 3; i++){
//		tester.increment();
//	}
//
//	BOOST_REQUIRE_EQUAL(0, tester.getValues());
//
//}


BOOST_AUTO_TEST_SUITE_END()

