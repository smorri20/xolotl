#include "DummyHardwareCounter.h"

using namespace xolotlPerf;

//DummyHardwareCounter::DummyHardwareCounter(std::string counterName,
//		const std::vector<HardwareQuantities> &counterQuantities) :
//		IHardwareCounter(counterName, counterQuantities),
//		quantities(0), values(0) {
//
//}


//DummyHardwareCounter::~DummyHardwareCounter() {
//
//}

//std::vector<int> DummyHardwareCounter::getValues() const {
//
//	return values;
//}
std::vector<double> DummyHardwareCounter::getValues(void) const {
    return std::vector<double>();
}

//const std::string DummyHardwareCounter::getName() const {
//
//	return "";
//}

std::vector<std::string> DummyHardwareCounter::getHardwareQuantities() const {

	return std::vector<std::string>();
}
