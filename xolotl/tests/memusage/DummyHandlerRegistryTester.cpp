#define BOOST_TEST_MODULE Regression

#include "mpi.h"
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <cmath>
#include <boost/test/included/unit_test.hpp>

#include "xolotlMemUsage/xolotlMemUsage.h"
#include "xolotlMemUsage/dummy/DummyHandlerRegistry.h"

namespace xmem = xolotlMemUsage;

// our coordinates in the MPI world
int cwRank = -1;
int cwSize = -1;

/**
 * Test suite for dummy handler registry.
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

        // Check for correct type.
        BOOST_REQUIRE_EQUAL(typeid(*(reg.get())).name(), typeid(xmem::DummyHandlerRegistry).name());

		BOOST_TEST_MESSAGE("Dummy handler registry created successfully.");
	} catch (const std::exception& e) {
		BOOST_TEST_MESSAGE(
				"DummyHandlerRegistry creation failed: " << e.what());
	}

	BOOST_REQUIRE_EQUAL(nGoodInits, 2U);
}

BOOST_AUTO_TEST_CASE(checkStubDummyFuncs) {

    xmem::initialize(xmem::IHandlerRegistry::dummy);
    auto reg = xmem::getHandlerRegistry();

    // Verify that dummy registry handler functions are stubs.
    std::ostringstream ostr;
    reg->reportData(ostr, reg->collectData());
    BOOST_REQUIRE_EQUAL(ostr.tellp(), 0U);
}

BOOST_AUTO_TEST_SUITE_END()
