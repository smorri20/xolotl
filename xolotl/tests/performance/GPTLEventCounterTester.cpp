#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <EventCounter.h>
#include <GPTLEventCounter.h>
#include <stdlib.h>
#include <stdio.h>
#include <papi.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the GPTLEventCounter.
 */
BOOST_AUTO_TEST_SUITE (GPTLEventCounter_testSuite)

//BOOST_AUTO_TEST_CASE(checkName) {
//
//	GPTLinitialize();
//
//	GPTLEventCounter tester("test");
//
//	std::cout << "\n" << "GPTLEventCounter Message: \n" << "tester.getName() " << tester.getName() << "\n"
//					  << std::endl;
//
//	// Require that the name of this GPTLEventCounter is "test"
//	BOOST_REQUIRE_EQUAL("test", tester.getName());
//}
//
//BOOST_AUTO_TEST_CASE(checkInitialValue) {
//
////	GPTLinitialize();
//	GPTLEventCounter tester("test");
//
//	std::cout << "\n" << "GPTLEventCounter Message: \n" << "tester.getValue() " << tester.getValue() << "\n" << std::endl;
//
//	// Require that the value of this GPTLEventCounter is 0
//	BOOST_REQUIRE_EQUAL(0, tester.getValue());
//
//}

BOOST_AUTO_TEST_CASE(checkCounting) {

	GPTLinitialize();
	GPTLEventCounter tester("test");

	std::cout << "\n" << "PAPI_VERSION = " << PAPI_VERSION_MAJOR(PAPI_VERSION) << "."
			  << PAPI_VERSION_MINOR(PAPI_VERSION) << "." << PAPI_VERSION_REVISION(PAPI_VERSION) << std::endl;

	std::cout << "\n" << "GPTLEventCounter Message: \n" << "tester.getName() = " << tester.getName() << "\n"
					  << std::endl;

	// Require that the name of this GPTLEventCounter is "test"
	BOOST_REQUIRE_EQUAL("test", tester.getName());

	std::cout << "\n" << "GPTLEventCounter Message: \n" << "tester.getValue() = " << tester.getValue() << "\n" << std::endl;

	// Require that the value of this GPTLEventCounter is 0
	BOOST_REQUIRE_EQUAL(0, tester.getValue());

	for(int i = 0; i < 3; i++){

		//increment the GPTLEventCounter
		tester.increment();
	}

	std::cout << "\n" << "GPTLEventCounter Message: \n" << "tester.getValue() = " << tester.getValue() << "\n" << std::endl;

	// Require that the value of this GPTLEventCounter is 3
	BOOST_REQUIRE_EQUAL(3, tester.getValue());

}


BOOST_AUTO_TEST_SUITE_END()





