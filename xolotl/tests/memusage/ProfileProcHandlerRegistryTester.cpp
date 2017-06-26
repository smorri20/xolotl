#define BOOST_TEST_MODULE Regression

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <typeinfo>
#include <system_error>
#include <unistd.h>
#include <boost/test/included/unit_test.hpp>

#include "xolotlMemUsage/xolotlMemUsage.h"
#include "xolotlMemUsage/profileproc/ProfileProcHandlerRegistry.h"

#include "tests/memusage/memUsageTestConfig.h"

#if defined(HAVE_STD_FILESYSTEM)
namespace fs = std::filesystem;
#elif defined(HAVE_BOOST_FILESYSTEM)
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#endif // defined(HAVE_STD_FILESYSTEM)

namespace xmem = xolotlMemUsage;

// our coordinates in the MPI world
int cwRank = -1;
int cwSize = -1;

/**
 * Test suite for ProfileProcHandlerRegistry.
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
		xmem::initialize(xmem::IHandlerRegistry::profileProc);
		nGoodInits++;

		std::shared_ptr<xmem::IHandlerRegistry> reg = xmem::getHandlerRegistry();
		if (reg) {
			nGoodInits++;
		}

        auto snreg = std::dynamic_pointer_cast<xmem::ProfileProcHandlerRegistry>(reg);
        if(snreg) {
            nGoodInits++;
        }

		BOOST_TEST_MESSAGE("Profile per-proc handler registry created successfully.");
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE("Profile per-proc HandlerRegistry creation failed: " << e.what());
	}

    BOOST_REQUIRE_EQUAL(nGoodInits, 3U);
}

BOOST_AUTO_TEST_CASE(accessStats) {
	try {
		xmem::initialize(xmem::IHandlerRegistry::profileProc);
        auto reg = xmem::getHandlerRegistry();
        BOOST_REQUIRE(reg);

        const std::string samplerName = "testSampler";

        std::shared_ptr<xmem::IMemUsageSampler> sampler = reg->getMemUsageSampler(samplerName);

		if (!sampler) {
			throw std::runtime_error("Failed to create MemUsageSampler");
		}

        // Ensure that we get the right type for statistics.
        auto gd = reg->collectData();
        BOOST_REQUIRE_EQUAL(typeid(*(gd.get())).name(), typeid(xmem::ProfileProcHandlerRegistry::GlobalMemUsageData).name());
        
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"Test of profile memory usage registry handler failed: " << e.what());
	}
}

#if defined(HAVE_STD_FILESYSTEM) || defined(HAVE_BOOST_FILESYSTEM)
BOOST_AUTO_TEST_CASE(genProfiles) {

	try {
		xmem::initialize(xmem::IHandlerRegistry::profileProc, 
                            std::chrono::duration<uint64_t, std::milli>(100),
                            "bogus.%r");
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

        // Generate the profiles.
        auto gd = reg->collectData();

        // Figure out the file name of our profile.
        std::ostringstream profilePathStr;
        uint32_t maxWidth = std::rint(std::ceil(std::log10(cwSize - 1)));
        profilePathStr << "bogus." << std::setw(maxWidth) << std::setfill('0') << cwRank;
        fs::path profileFilename(profilePathStr.str());
        fs::path profilePath = fs::current_path() / profileFilename;

        // Verify that we generated a profile.
        BOOST_CHECK(fs::exists(profilePath));
        BOOST_CHECK(fs::is_regular_file(profilePath));
        BOOST_CHECK(not fs::is_empty(profilePath));

        // Everyone tries to remove their own file.
        if(fs::exists(profilePath)) {
            fs::remove(profilePath);
        }

	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"Test of generating memory usage profiles failed: " << e.what());
	}
}
#endif // defined(HAVE_STD_FILESYSTEM) || defined(HAVE_BOOST_FILESYSTEM)

BOOST_AUTO_TEST_SUITE_END()
