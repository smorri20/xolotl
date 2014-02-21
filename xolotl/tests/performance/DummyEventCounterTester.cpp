#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <EventCounter.h>
#include <DummyEventCounter.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the DummyEventCounter.
 */
BOOST_AUTO_TEST_SUITE (DummyEventCounter_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	DummyEventCounter tester("test");
//	EventCounter * testme = &tester;

	BOOST_REQUIRE_EQUAL("", tester.getName());
//}
//
//BOOST_AUTO_TEST_CASE(checkInitialValue) {
//
//	DummyEventCounter tester("test");
//	EventCounter * testme = &tester;
//
//	BOOST_REQUIRE_EQUAL(0, testme->getValue());
//
//}
//
//BOOST_AUTO_TEST_CASE(checkCounting) {
//
//	DummyEventCounter tester("test");
//	EventCounter * testme = &tester;
//
//	long count = 3;
//
//	for(int i = 0; i < 3; i++){
//		testme->increment();
//	}
//
//	BOOST_REQUIRE_EQUAL(0, testme->getValue());
//
}


BOOST_AUTO_TEST_SUITE_END()

