#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <limits>
#include <string>
#include <cstdlib>


#include "xolotlMemUsage/summaryproc/MemUsageSampler.h"

using namespace std;
using namespace xolotlMemUsage;

/**
 * This suite is responsible for testing the [PerProc]MemUsageSampler.
 */
BOOST_AUTO_TEST_SUITE (MemUsageSampler_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

    std::string samplerName = "testSampler";

	MemUsageSampler sampler(samplerName);

	BOOST_TEST_MESSAGE(
			"\n" << "MemUsageSampler Message: \n" << "sampler.GetName() " << sampler.GetName() << "\n");

	// Check name of the sampler object.
	BOOST_REQUIRE_EQUAL(samplerName, sampler.GetName());
}


void
CheckMemUsagePageStatsInitialValue(std::string memberName, 
                                    const MemUsagePageStats& stats)
{
    BOOST_REQUIRE_EQUAL(std::numeric_limits<uint64_t>::max(), stats.min);
    BOOST_REQUIRE_EQUAL(std::numeric_limits<uint64_t>::min(), stats.max);
    BOOST_REQUIRE(std::isnan(stats.average));
    BOOST_REQUIRE(std::isnan(stats.stdev));
}



BOOST_AUTO_TEST_CASE(checkInitialValue) {

	MemUsageSampler sampler("test");

    BOOST_TEST_MESSAGE("\nChecking MemUsageSampler initial values.\n");

	// Check initial values.
    auto testVal = std::dynamic_pointer_cast<MemUsageStats>(sampler.getValue());
    BOOST_REQUIRE(testVal);
    CheckMemUsagePageStatsInitialValue("vmSize", testVal->vmSize);
    CheckMemUsagePageStatsInitialValue("vmRSS", testVal->vmRSS);
    CheckMemUsagePageStatsInitialValue("rss", testVal->rss);
    CheckMemUsagePageStatsInitialValue("text", testVal->text);
    CheckMemUsagePageStatsInitialValue("dataAndStack", testVal->dataAndStack);
}

inline
void
CheckMinMaxAvg(std::string metric, const MemUsagePageStats& stats)
{
    BOOST_TEST_MESSAGE("\nChecking page metric statistic relationships for " << metric << "\n");
    BOOST_REQUIRE(stats.min <= stats.max);
    BOOST_REQUIRE(stats.min <= stats.average);
    BOOST_REQUIRE(stats.average <= stats.max);
}


BOOST_AUTO_TEST_CASE(checkSampling) {

    std::string samplerName = "testSampler";

	MemUsageSampler sampler(samplerName);


    // Collect some samples.
	BOOST_TEST_MESSAGE("\nSampling memory usage functionality.\n");
    sampler.start();
    // What to do here?
    sleep(5);
    sampler.stop();
    BOOST_TEST_MESSAGE("done.");

    // Verify statistics collected, as much as possible.
    auto& runningSampleData = sampler.GetRunningSampleData();
    BOOST_REQUIRE_CLOSE(5.0, (double)(runningSampleData.nSamples), 20.0);

    auto testVal = std::dynamic_pointer_cast<MemUsageStats>(sampler.getValue());
    BOOST_REQUIRE(testVal);

    CheckMinMaxAvg("vmSize", testVal->vmSize);
    CheckMinMaxAvg("vmRSS", testVal->vmRSS);
    CheckMinMaxAvg("rss", testVal->rss);
    CheckMinMaxAvg("text", testVal->text);
    CheckMinMaxAvg("dataAndStack", testVal->dataAndStack);

    BOOST_REQUIRE(testVal->vmSize.average > (testVal->text.average + testVal->dataAndStack.average));
}

BOOST_AUTO_TEST_SUITE_END()

