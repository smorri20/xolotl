#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DummyHardwareCounter.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the DummyHardwareCounter.
 */
BOOST_AUTO_TEST_SUITE (DummyHardwareCounter_testSuite)

const std::vector<HardwareQuantities> test_hquan = {L1_CACHE_MISS,L2_CACHE_MISS,L3_CACHE_MISS,BRANCH_MISPRED,TOTAL_CYCLES,TOTAL_INSTRUC,FLPT_INSTRUC};

BOOST_AUTO_TEST_CASE(checkName) {

	DummyHardwareCounter tester("test",test_hquan);

	BOOST_REQUIRE_EQUAL("", tester.getName());
}

BOOST_AUTO_TEST_CASE(checkInitialValue) {

	DummyHardwareCounter tester("test",test_hquan);

	for (unsigned i = 0; i < tester.getValues().size(); i++)
	    BOOST_REQUIRE_EQUAL(0, tester.getValues().at(i));

}

BOOST_AUTO_TEST_CASE(checkCounting) {

	DummyHardwareCounter tester("test",test_hquan);;

	for (unsigned i = 0; i < tester.getValues().size(); i++)
	{
		tester.increment();
	    BOOST_REQUIRE_EQUAL(0, tester.getValues().at(i));
	}

}


BOOST_AUTO_TEST_SUITE_END()

