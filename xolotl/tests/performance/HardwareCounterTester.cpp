//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MODULE Regression
//
//#include <boost/test/included/unit_test.hpp>
//#include <HardwareCounter.h>
//#include <HardwareQuantities.h>
//#include <string>
//#include <chrono>
//#include <ctime>
//
//using namespace std;
//using namespace xolotlPerf;
//
//class HardwareCounterTester : public HardwareCounter {
//private:
//	std::vector<long> values;
//public:
//	HardwareCounterTester(const std::string& name, const std::vector<HardwareQuantities>& quantities) :
//		HardwareCounter(name, quantities) {values.assign((quantities.size()),0);}
//
//	~HardwareCounterTester() {};
//
////	std::vector<long> getValues() const {
//////		for(std::vector<long>::iterator it = values.begin(); it !=values.end(); ++it)
//////			return *it;
////		return values;
////	}
//	long getValues() const {
////		for(std::vector<long>::iterator it = values.begin(); it !=values.end(); ++it)
////			return *it;
//		return values.at(FLPT_INSTRUC);
//	}
//
//	void increment() {
////		for(std::vector<long>::iterator it = values.begin(); it !=values.end(); ++it)
////			++(*it);
//		++values.at(FLPT_INSTRUC);
//	}
//
//};
//
///**
// * This suite is responsible for testing the HardwareCounter.
// */
//BOOST_AUTO_TEST_SUITE (HardwareCounter_testSuite)
//
//BOOST_AUTO_TEST_CASE(checkName) {
//
//	std::vector<HardwareQuantities> hwquant(7);
//	HardwareCounterTester tester("test", hwquant);
//	HardwareCounter * testme = &tester;
//
//	BOOST_REQUIRE_EQUAL("test", testme->getName());
//}
//
//BOOST_AUTO_TEST_CASE(checkInitialValue) {
//
//	std::vector<HardwareQuantities> hwquant;
//	HardwareCounterTester tester("test", hwquant);
//	HardwareCounter * testme = &tester;
//
//	BOOST_TEST_MESSAGE("testme->getValues is: "
//				<< testme->getValues());
//
////	BOOST_REQUIRE_EQUAL(0, testme->getValues());
//
//}
//
//BOOST_AUTO_TEST_CASE(checkCounting) {
//
//	std::vector<HardwareQuantities> hwquant;
//	HardwareCounterTester tester("test", hwquant);
//	HardwareCounter * testme = &tester;
//	long count = 3;
//
//	for(int i = 0; i < 3; i++){
//		testme->increment();
//	}
//
//	BOOST_REQUIRE_EQUAL(count, testme->getValues());
//
//}
//
//
//BOOST_AUTO_TEST_SUITE_END()
