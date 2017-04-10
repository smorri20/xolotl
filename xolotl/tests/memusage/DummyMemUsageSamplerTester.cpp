#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <string>
#include <limits>
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


void
CheckMemUsagePageStatsInitialValue(std::string memberName, 
                                    const MemUsagePageStats& stats)
{
    BOOST_REQUIRE_EQUAL(std::numeric_limits<uint64_t>::max(), stats.min);
    BOOST_REQUIRE_EQUAL(std::numeric_limits<uint64_t>::min(), stats.max);
    BOOST_REQUIRE(std::isnan(stats.avg));
    BOOST_REQUIRE(std::isnan(stats.stdev));
}


BOOST_AUTO_TEST_CASE(checkInitialValue) {

	DummyMemUsageSampler sampler("test");

    CheckMemUsagePageStatsInitialValue("vmSize", sampler.getValue().vmSize);
    CheckMemUsagePageStatsInitialValue("vmRSS", sampler.getValue().vmRSS);
    CheckMemUsagePageStatsInitialValue("rss", sampler.getValue().rss);
    CheckMemUsagePageStatsInitialValue("text", sampler.getValue().text);
    CheckMemUsagePageStatsInitialValue("dataAndStack", sampler.getValue().dataAndStack);
}

BOOST_AUTO_TEST_CASE(checkSampling) {

	DummyMemUsageSampler sampler("test");

	sampler.start();
	sleep(3);
	sampler.stop();

    // These are dummies, so their values after being started and stopped
    // should be the same as their initial values.
    CheckMemUsagePageStatsInitialValue("vmSize", sampler.getValue().vmSize);
    CheckMemUsagePageStatsInitialValue("vmRSS", sampler.getValue().vmRSS);
    CheckMemUsagePageStatsInitialValue("rss", sampler.getValue().rss);
    CheckMemUsagePageStatsInitialValue("text", sampler.getValue().text);
    CheckMemUsagePageStatsInitialValue("dataAndStack", sampler.getValue().dataAndStack);
}

BOOST_AUTO_TEST_SUITE_END()

