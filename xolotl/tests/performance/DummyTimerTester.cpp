#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <Timer.h>
#include <DummyTimer.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the DummyTimer.
 */
BOOST_AUTO_TEST_SUITE (DummyTimer_testSuite)

BOOST_AUTO_TEST_CASE(checkName) {

	DummyTimer tester("test");
	Timer * testme = &tester;

	BOOST_REQUIRE_EQUAL("", testme->getName());
}

BOOST_AUTO_TEST_CASE(checkInitialValue) {

	DummyTimer tester("test");
	Timer * testme = &tester;

	BOOST_REQUIRE_EQUAL(0, testme->getValue());

}

BOOST_AUTO_TEST_CASE(checkTiming) {

	DummyTimer tester("test");
	Timer * testme = &tester;

	testme->start();
	sleep(3);
	testme->stop();

	BOOST_REQUIRE_EQUAL(0, testme->getValue());
	BOOST_REQUIRE_EQUAL(0, testme->getUnits());

}


BOOST_AUTO_TEST_SUITE_END()




