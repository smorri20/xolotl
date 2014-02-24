#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <Timer.h>
#include <GPTLTimer.h>
#include <string>

using namespace std;
using namespace xolotlPerf;

/**
 * This suite is responsible for testing the GPTLTimer.
 */
BOOST_AUTO_TEST_SUITE (GPTLTimer_testSuite)

//BOOST_AUTO_TEST_CASE(checkName) {
//	GPTLinitialize();
//
//	GPTLTimer tester("test");
//	Timer * testme = &tester;
//
//	BOOST_REQUIRE_EQUAL("test", testme->getName());
//}
//
//BOOST_AUTO_TEST_CASE(checkInitialValue) {
//
//	GPTLTimer tester("test");
//	Timer * testme = &tester;
//
//	BOOST_REQUIRE_EQUAL(0.0, testme->getValue());
//
//}

BOOST_AUTO_TEST_CASE(checkTiming) {

	GPTLinitialize();

	GPTLTimer testme("test");
//	Timer * testme = &tester;

//	Timer testme("test")= new GPTLTimer;

	BOOST_REQUIRE_EQUAL(0.0, testme.getValue());

	BOOST_TEST_MESSAGE("GPTLTimer Message: \n"
			  << "testme.getName() " << testme.getName() << "\n"
			  << "testme.getValue() " << testme.getValue() << "\n"
			  << "\n");

	std::cout << "\n" << "GPTLTimer Message: \n" << "testme.getName() " << testme.getName() << "\n"
			  << "testme.getValue() " << testme.getValue() << "\n" << std::endl;

	testme.start();
	sleep(1);
	testme.stop();

	std::cout << "\n" << "GPTLTimer Message: \n" << "testme.getName() " << testme.getName() << "\n"
			  << "testme.getValue() " << testme.getValue() << "\n" << std::endl;

	BOOST_REQUIRE_EQUAL("test", testme.getName());

	BOOST_REQUIRE_CLOSE(1, testme.getValue(),0.0001);
//	BOOST_REQUIRE_EQUAL(0, testme->getUnits());

}


BOOST_AUTO_TEST_SUITE_END()








//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MODULE Regression
//
//#include <boost/test/included/unit_test.hpp>
//#include "Timer.h"
//#include <string>
//#include <chrono>
//#include <ctime>
//
//using namespace std;
//using namespace xolotlPerf;
//
//class TimerTester : public Timer {
//protected:
//
//	/**
//	 * The default constructor in private because TestEventCounters
//	 * must always be initialized with a name.
//	 */
////	TimerTester():value(0),Timer("") {}
//	std::chrono::high_resolution_clock::time_point t_start, t_stop;
////	long units;
//public:
//	TimerTester(const std::string name) : Timer(name) {
////		units = 0;
//	}
//	~TimerTester() {};
//
//	void start() {
//		t_start = std::chrono::high_resolution_clock::now();
//	}
//
//	void stop() {
//		t_stop = std::chrono::high_resolution_clock::now();
//	}
//
//	long getValue() const {
//		return std::chrono::duration_cast<std::chrono::seconds>(t_stop-t_start).count();
//	}
////	long getUnits() {return units;}
//
//};
//
///**
// * This suite is responsible for testing the Timer.
// */
//BOOST_AUTO_TEST_SUITE (Timer_testSuite)
//
//BOOST_AUTO_TEST_CASE(checkName) {
//
//	TimerTester tester("test");
//	Timer * testme = &tester;
//
//	BOOST_REQUIRE_EQUAL("test", testme->getName());
//}
//
//BOOST_AUTO_TEST_CASE(checkInitialValue) {
//
//	TimerTester tester("test");
//	Timer * testme = &tester;
//
//	BOOST_REQUIRE_EQUAL(0, testme->getValue());
//
//}
//
//BOOST_AUTO_TEST_CASE(checkTiming) {
//
//	//Local Declarations
//	TimerTester tester("test");
//	Timer * testme = &tester;
//	long sleeper = 5;
//
////	clock_t t = clock();
////	sleep(5);
////	t = clock() - t;
////	float val = ((float)t)/CLOCKS_PER_SEC;
//
//	testme->start();
//	sleep(5);
//	testme->stop();
//
////	BOOST_TEST_MESSAGE("TimerTester Message: \n" << "Wallclock time passed " << std::chrono::duration_cast<std::chrono::seconds>(testme->start() - testme->stop()).count() );
//
//	BOOST_REQUIRE_EQUAL(sleeper, testme->getValue());
//}
//
//BOOST_AUTO_TEST_SUITE_END()
