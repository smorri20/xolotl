#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSICluster.h>
#include "SimpleReactionNetwork.h"
#include <HeVCluster.h>
#include <SuperCluster.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <algorithm>
#include <math.h>
#include <limits>

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry = std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the HeVCluster.
 */
BOOST_AUTO_TEST_SUITE(SuperCluster_testSuite)


BOOST_AUTO_TEST_CASE(getSpeciesSize) {
	SuperCluster cluster(4.5, 5.0, 2, 0.1, 200.0, registry);

	// Get the composition back
	auto composition = cluster.getComposition();

	// Check the composition is the created one
	BOOST_REQUIRE_EQUAL(composition["He"], 9);
	BOOST_REQUIRE_EQUAL(composition["V"], 10);
	BOOST_REQUIRE_EQUAL(composition["I"], 0);

	// Check if it is a mixed cluster
	BOOST_REQUIRE_EQUAL(cluster.isMixed(), true);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
