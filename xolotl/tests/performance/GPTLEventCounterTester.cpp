#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <EventCounter.h>
#include <GPTLEventCounter.h>
#include "papi.h"
#include <string>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the GPTLEventCounter.
 */
BOOST_AUTO_TEST_SUITE (GPTLEventCounter_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	GPTLinitialize();

	GPTLEventCounter tester("test");
//	EventCounter * testme = &tester;

	BOOST_REQUIRE_EQUAL("test", tester.getName());
}

BOOST_AUTO_TEST_CASE(checkInitialValue) {

	GPTLinitialize();
	GPTLEventCounter tester("test");
//	EventCounter * testme = &tester;

	BOOST_REQUIRE_EQUAL(0, tester.getValue());

}

BOOST_AUTO_TEST_CASE(checkCounting) {

	GPTLinitialize();
	GPTLEventCounter tester("test");
//	EventCounter * testme = &tester;

	for(int i = 0; i < 3; i++){
		tester.increment();
	}

	BOOST_REQUIRE_EQUAL(3, tester.getValue());

}


BOOST_AUTO_TEST_SUITE_END()





