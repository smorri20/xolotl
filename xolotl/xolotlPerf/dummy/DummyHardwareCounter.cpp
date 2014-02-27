#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include "DummyHardwareCounter.h"

using namespace xolotlPerf;

DummyHardwareCounter::DummyHardwareCounter(std::string counterName, const std::vector<HardwareQuantities> &counterQuantities) :
		IHardwareCounter(counterName, counterQuantities) {

//	name = counterName;

}


DummyHardwareCounter::~DummyHardwareCounter() {

}

//std::vector<int> DummyHardwareCounter::getValues() const {
//
//	return values;
//}

const std::string DummyHardwareCounter::getName() const {

	return "";
}

std::vector<HardwareQuantities> DummyHardwareCounter::getHardwareQuantities() const {

	return quantities;
}
