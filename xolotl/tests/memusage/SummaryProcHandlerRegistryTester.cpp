#define BOOST_TEST_MODULE Regression

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <typeinfo>
#include <unistd.h>
#include <boost/test/included/unit_test.hpp>
#include "xolotlMemUsage/xolotlMemUsage.h"
#include "xolotlMemUsage/summaryproc/SummaryProcHandlerRegistry.h"

namespace xmem = xolotlMemUsage;

// our coordinates in the MPI world
int cwRank = -1;
int cwSize = -1;

/**
 * Test suite for SummaryProcHandlerRegistry.
 */
BOOST_AUTO_TEST_SUITE (MemUsageHandlerRegistry_testSuite)

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


BOOST_AUTO_TEST_CASE(createHandlerReg) {
	unsigned int nGoodInits = 0;

	try {
		xmem::initialize(xmem::IHandlerRegistry::summaryProc);
		nGoodInits++;

		std::shared_ptr<xmem::IHandlerRegistry> reg = xmem::getHandlerRegistry();
		if (reg) {
			nGoodInits++;
		}

        auto snreg = std::dynamic_pointer_cast<xmem::SummaryProcHandlerRegistry>(reg);
        if(snreg) {
            nGoodInits++;
        }

		BOOST_TEST_MESSAGE("Summary per-proc handler registry created successfully.");
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE("Summary per-proc HandlerRegistry creation failed: " << e.what());
	}

    BOOST_REQUIRE_EQUAL(nGoodInits, 3U);
}

BOOST_AUTO_TEST_CASE(accessStats) {
	try {
		xmem::initialize(xmem::IHandlerRegistry::summaryProc);
        auto reg = xmem::getHandlerRegistry();
        BOOST_REQUIRE(reg);

        const std::string samplerName = "testSampler";

        std::shared_ptr<xmem::IMemUsageSampler> sampler = reg->getMemUsageSampler(samplerName);

		if (!sampler) {
			throw std::runtime_error("Failed to create MemUsageSampler");
		}

        // Ensure that we get the right type for statistics.
        auto gd = reg->collectData();
        BOOST_REQUIRE_EQUAL(typeid(*(gd.get())).name(), typeid(xmem::SummaryProcHandlerRegistry::GlobalMemUsageStats).name());
        
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"Test of summary memory usage registry handler failed: " << e.what());
	}
}

BOOST_AUTO_TEST_CASE(aggregateStats) {
	try {
		xmem::initialize(xmem::IHandlerRegistry::summaryProc);
        auto reg = xmem::getHandlerRegistry();

        const std::string samplerName = "testSampler";
        std::shared_ptr<xmem::IMemUsageSampler> sampler = reg->getMemUsageSampler(samplerName);

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
        auto rawStats = reg->collectData();

		// Verify the statistics collected.
        auto stats = std::dynamic_pointer_cast<xmem::SummaryProcHandlerRegistry::GlobalMemUsageStats>(rawStats);
        BOOST_REQUIRE(stats);
        
        // Only rank 0 should have valid data.
        if(cwRank == 0) {

            // First check times.  Should be very close to the nTimedSeconds
            // with little spread.
            BOOST_REQUIRE_EQUAL(stats->memStats.size(), 1U);

            xmem::MemUsageStats const& mapMemStats = stats->memStats[samplerName].stats;
            auto samplerMemStats = std::dynamic_pointer_cast<xmem::MemUsageStats>(sampler->getValue());
            BOOST_REQUIRE(samplerMemStats);


            BOOST_TEST_MESSAGE("mem usage sampler name: " << stats->memStats.begin()->first);
            BOOST_REQUIRE_EQUAL(stats->memStats.begin()->first, samplerName);

#if defined(HAVE_STATM)
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.vmSize.min), (double)(samplerMemStats->vmSize.min), 10.0);
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.vmSize.max), (double)(samplerMemStats->vmSize.max), 3.0);
            BOOST_REQUIRE_CLOSE(mapMemStats.vmSize.average, samplerMemStats->vmSize.average, 3.0);

            BOOST_REQUIRE_CLOSE((double)(mapMemStats.vmRSS.min), (double)(samplerMemStats->vmRSS.min), 10.0);
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.vmRSS.max), (double)(samplerMemStats->vmRSS.max), 3.0);
            BOOST_REQUIRE_CLOSE(mapMemStats.vmRSS.average, samplerMemStats->vmRSS.average, 3.0);

#else
#error "Configuration error: thought we had a per-proc data source, but none is actually available."
#endif // defined(HAVE_SYSINFO)

        }
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"Test of aggregating memory usage stats failed: " << e.what());
	}
}

BOOST_AUTO_TEST_SUITE_END()
