#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <string>
#include <limits>
#include <typeinfo>
#include <boost/test/included/unit_test.hpp>
#include "xolotlMemUsage/dummy/DummyMemUsageSampler.h"

using namespace std;
using namespace xolotlMemUsage;

/**
 * This suite is responsible for testing the DummyMemUsageSampler.
 */
BOOST_AUTO_TEST_SUITE (DummyMemUsageSampler_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	DummyMemUsageSampler sampler("test");

	BOOST_REQUIRE_EQUAL("unused", sampler.getName());
}


BOOST_AUTO_TEST_CASE(checkSampling) {

	DummyMemUsageSampler sampler("test");

	sampler.start();
	sleep(3);
	sampler.stop();

    // These are dummies, so the value can be empty.
    BOOST_REQUIRE_EQUAL(typeid(*(sampler.getValue().get())).name(), typeid(IMemUsageSampler::MemUsageData).name());
}

BOOST_AUTO_TEST_SUITE_END()
