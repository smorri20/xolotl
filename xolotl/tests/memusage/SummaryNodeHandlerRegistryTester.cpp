#define BOOST_TEST_MODULE Regression

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <typeinfo>
#include <unistd.h>
#include <boost/test/included/unit_test.hpp>
#include "xolotlMemUsage/xolotlMemUsage.h"
#include "xolotlMemUsage/summarynode/SummaryNodeHandlerRegistry.h"

namespace xmem = xolotlMemUsage;

// our coordinates in the MPI world
int cwRank = -1;
int cwSize = -1;

/**
 * Test suite for SummaryNodeHandlerRegistry.
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
		xmem::initialize(xmem::IHandlerRegistry::summaryNode);
		nGoodInits++;

		std::shared_ptr<xmem::IHandlerRegistry> reg = xmem::getHandlerRegistry();
		if (reg) {
			nGoodInits++;
		}

        auto snreg = std::dynamic_pointer_cast<xmem::SummaryNodeHandlerRegistry>(reg);
        if(snreg) {
            nGoodInits++;
        }

		BOOST_TEST_MESSAGE("Summary per-node handler registry created successfully.");
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE("Summary per-node HandlerRegistry creation failed: " << e.what());
	}

    // One rank on each node should get the actual registry.
    // We know rank 0 will be one, but can't say anything about
    // which others will without knowing how many processes were
    // created on each node.
    // TODO this only works if running the test on a single node.
    // How can we tell who are supposed to have a "real" handler,
    // without duplicating the entire logic that it uses?
    if(cwRank == 0) {
	    BOOST_REQUIRE_EQUAL(nGoodInits, 3U);
    }
#if READY    
    else {
	    BOOST_REQUIRE_EQUAL(nGoodInits, 2U);
    }
#endif // READY
}

BOOST_AUTO_TEST_CASE(accessStats) {
	try {
		xmem::initialize(xmem::IHandlerRegistry::summaryNode);
        auto reg = xmem::getHandlerRegistry();
        BOOST_REQUIRE(reg);

        const std::string samplerName = "testSampler";

        std::shared_ptr<xmem::IMemUsageSampler> sampler = reg->getMemUsageSampler(samplerName);

		if (!sampler) {
			throw std::runtime_error("Failed to create MemUsageSampler");
		}

        // Ensure that we get the right type for statistics.
        auto gd = reg->collectData();
        if(cwRank == 0) {
            BOOST_REQUIRE_EQUAL(typeid(*(gd.get())).name(), typeid(xmem::SummaryNodeHandlerRegistry::GlobalMemUsageStats).name());
        }
        
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"Test of summary memory usage registry handler failed: " << e.what());
	}
}

BOOST_AUTO_TEST_CASE(aggregateStats) {
	try {
		xmem::initialize(xmem::IHandlerRegistry::summaryNode);
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
		// Only rank 0 should have valid data.
		if (cwRank == 0) {

            auto stats = std::dynamic_pointer_cast<xmem::SummaryNodeHandlerRegistry::GlobalMemUsageStats>(rawStats);
            
			// First check times.  Should be very close to the nTimedSeconds
			// with little spread.
			BOOST_REQUIRE_EQUAL(stats->memStats.size(), 1U);

            xmem::NodeMemUsageStats const& mapMemStats = stats->memStats[samplerName].stats;
            auto samplerMemStats = std::dynamic_pointer_cast<xmem::NodeMemUsageStats>(sampler->getValue());
            BOOST_REQUIRE(samplerMemStats);


			BOOST_TEST_MESSAGE("mem usage sampler name: " << stats->memStats.begin()->first);
			BOOST_REQUIRE_EQUAL(stats->memStats.begin()->first, samplerName);

#if defined(HAVE_SYSINFO)
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.totalram.min), (double)(samplerMemStats->totalram.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.totalram.max), (double)(samplerMemStats->totalram.max), 3.0);
            BOOST_REQUIRE_CLOSE(mapMemStats.totalram.average, samplerMemStats->totalram.average, 3.0);

            BOOST_REQUIRE_CLOSE((double)(mapMemStats.freeram.min), (double)(samplerMemStats->freeram.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.freeram.max), (double)(samplerMemStats->freeram.max), 3.0);
            BOOST_REQUIRE_CLOSE(mapMemStats.freeram.average, samplerMemStats->freeram.average, 3.0);

#elif defined(HAVE_MACH_HOST_STATISTICS)
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.free_count.min), (double)(samplerMemStats->free_count.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.free_count.max), (double)(samplerMemStats->free_count.max), 3.0);
            BOOST_REQUIRE_CLOSE(mapMemStats.free_count.average, samplerMemStats->free_count.average, 3.0);

            BOOST_REQUIRE_CLOSE((double)(mapMemStats.active_count.min), (double)(samplerMemStats->active_count.min), 3.0);
            BOOST_REQUIRE_CLOSE((double)(mapMemStats.active_count.max), (double)(samplerMemStats->active_count.max), 3.0);
            BOOST_REQUIRE_CLOSE(mapMemStats.active_count.average, samplerMemStats->active_count.average, 3.0);

#else
#error "Configuration error: thought we had a per-node data source, but none is actually available."
#endif // defined(HAVE_SYSINFO)

		}
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"Test of aggregating memory usage stats failed: " << e.what());
	}
}

BOOST_AUTO_TEST_SUITE_END()
