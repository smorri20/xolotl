#define BOOST_TEST_MODULE Regression

#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <unistd.h>
#include <boost/test/included/unit_test.hpp>
#include "xolotlMemUsage/xolotlMemUsage.h"
#include "xolotlMemUsage/MemUsageObjStatistics.h"

namespace xmem = xolotlMemUsage;

// our coordinates in the MPI world
int cwRank = -1;
int cwSize = -1;

/**
 * Test suite for HandlerRegistry classes (mainly StdHandlerRegistry).
 */
BOOST_AUTO_TEST_SUITE (StdHandlerRegistry_testSuite)

struct MPIFixture {
	MPIFixture(void) {
		MPI_Init(&boost::unit_test::framework::master_test_suite().argc,
				&boost::unit_test::framework::master_test_suite().argv);

		MPI_Comm_rank(MPI_COMM_WORLD, &cwRank);
		MPI_Comm_size(MPI_COMM_WORLD, &cwSize);
	}

	~MPIFixture(void) {
		MPI_Finalize();
	}
};

#if BOOST_VERSION >= 105900
// In Boost 1.59, the semicolon at the end of the definition of BOOST_GLOBAL_FIXTURE is removed
BOOST_GLOBAL_FIXTURE(MPIFixture);
#else
// With earlier Boost versions, naively adding a semicolon to our code will generate compiler
// warnings about redundant semicolons
BOOST_GLOBAL_FIXTURE (MPIFixture)
#endif

BOOST_AUTO_TEST_CASE(createDummyHandlerReg) {
	unsigned int nGoodInits = 0;

	try {
		xmem::initialize(xmem::IHandlerRegistry::dummy);
		nGoodInits++;

		std::shared_ptr<xmem::IHandlerRegistry> reg =
				xmem::getHandlerRegistry();
		if (reg) {
			nGoodInits++;
		}

		BOOST_TEST_MESSAGE("Dummy handler registry created successfully.");
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"DummyHandlerRegistry creation failed: " << e.what());
	}

	BOOST_REQUIRE_EQUAL(nGoodInits, 2U);
}

BOOST_AUTO_TEST_CASE(createStdHandlerReg) {
	unsigned int nGoodInits = 0;

	try {
		xmem::initialize(xmem::IHandlerRegistry::std);
		nGoodInits++;

		std::shared_ptr<xmem::IHandlerRegistry> reg =
				xmem::getHandlerRegistry();
		if (reg) {
			nGoodInits++;
		}

		BOOST_TEST_MESSAGE("Standard handler registry created successfully.");
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE("StdHandlerRegistry creation failed: " << e.what());
	}

	BOOST_REQUIRE_EQUAL(nGoodInits, 2U);
}

BOOST_AUTO_TEST_CASE(aggregateStats) {
	try {
		xmem::initialize(xmem::IHandlerRegistry::std);
		std::shared_ptr<xmem::IHandlerRegistry> reg =
				xmem::getHandlerRegistry();

        const std::string samplerName = "testSampler";

        std::shared_ptr<xmem::IMemUsageSampler> sampler = reg->getMemUsageSampler(samplerName);

		if (!sampler) {
			throw std::runtime_error("Failed to create MemUsageSampler");
		}

        // Simulate sampling memory for some event.
        // TODO How can we do something that we can verify?
		BOOST_TEST_MESSAGE("Simulating memory usage...");
        sampler->start();
        // What to do here?
        sleep(3);
        sampler->stop();
		const unsigned int nTimedSeconds = 5;
		BOOST_TEST_MESSAGE("done.");

		// compute statistics about the program's memory usage.
		auto memStats = reg->collectStatistics();

		// Verify the statistics collected.
		// Only rank 0 does the verification.
		if (cwRank == 0) {

			// First check times.  Should be very close to the nTimedSeconds
			// with little spread.
			BOOST_REQUIRE_EQUAL(memStats.memStats.size(), 1U);
			xmem::MemUsageObjStatistics<xmem::IMemUsageSampler::ValType>& memStatsObj =
					memStats.memStats.begin()->second;

			BOOST_TEST_MESSAGE("mem usage sampler name: " << memStats.memStats.begin()->first);
			BOOST_REQUIRE_EQUAL(memStats.memStats.begin()->first, samplerName);

            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.vmSize.min), (double)(sampler->getValue().vmSize.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.vmSize.max), (double)(sampler->getValue().vmSize.max), 3.0);
            BOOST_REQUIRE_CLOSE(memStats.memStats[samplerName].stats.vmSize.avg, sampler->getValue().vmSize.avg, 3.0);

            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.vmRSS.min), (double)(sampler->getValue().vmRSS.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.vmRSS.max), (double)(sampler->getValue().vmRSS.max), 3.0);
            BOOST_REQUIRE_CLOSE(memStats.memStats[samplerName].stats.vmRSS.avg, sampler->getValue().vmRSS.avg, 3.0);

            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.rss.min), (double)(sampler->getValue().rss.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.rss.max), (double)(sampler->getValue().rss.max), 3.0);
            BOOST_REQUIRE_CLOSE(memStats.memStats[samplerName].stats.rss.avg, sampler->getValue().rss.avg, 3.0);

            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.text.min), (double)(sampler->getValue().text.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.text.max), (double)(sampler->getValue().text.max), 3.0);
            BOOST_REQUIRE_CLOSE(memStats.memStats[samplerName].stats.text.avg, sampler->getValue().text.avg, 3.0);

            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.dataAndStack.min), (double)(sampler->getValue().dataAndStack.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(memStats.memStats[samplerName].stats.dataAndStack.max), (double)(sampler->getValue().dataAndStack.max), 3.0);
            BOOST_REQUIRE_CLOSE(memStats.memStats[samplerName].stats.dataAndStack.avg, sampler->getValue().dataAndStack.avg, 3.0);
		}
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"Test of aggregating memory usage stats failed: " << e.what());
	}
}

BOOST_AUTO_TEST_SUITE_END()
