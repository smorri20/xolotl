#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DummyHardwareCounter.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

const std::vector<HardwareQuantities> test_hquan = {L1_CACHE_MISS,L2_CACHE_MISS,L3_CACHE_MISS,BRANCH_MISPRED,TOTAL_CYCLES,TOTAL_INSTRUC,FLPT_INSTRUC};

/**
 * This suite is responsible for testing the DummyHardwareCounter.
 */
BOOST_AUTO_TEST_SUITE (DummyHardwareCounter_testSuite)

//const std::vector<HardwareQuantities> test_hquan = {L1_CACHE_MISS,L2_CACHE_MISS,L3_CACHE_MISS,BRANCH_MISPRED,TOTAL_CYCLES,TOTAL_INSTRUC,FLPT_INSTRUC};

BOOST_AUTO_TEST_CASE(checkName) {

	DummyHardwareCounter tester("test",test_hquan);

	BOOST_REQUIRE_EQUAL("test", tester.getName());
}

BOOST_AUTO_TEST_SUITE_END()

