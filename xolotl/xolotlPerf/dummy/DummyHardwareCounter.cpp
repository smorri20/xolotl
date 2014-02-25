#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include "DummyHardwareCounter.h"

using namespace xolotlPerf;

DummyHardwareCounter::DummyHardwareCounter(std::string aname, const std::vector<HardwareQuantities> &hquantities) :
		HardwareCounter(aname, hquantities) {

	name = aname;
	values = std::vector<int>(hquantities.size(), 0);

}


DummyHardwareCounter::~DummyHardwareCounter() {

}

//std::shared_ptr<HardwareCounter> DummyHardwareCounter::clone() {
//	std::shared_ptr<HardwareCounter> eventCounter(new DummyHardwareCounter(*this));
//	return eventCounter;
//}

std::vector<int> DummyHardwareCounter::getValues() const {

	return values;
}

const std::string DummyHardwareCounter::getName() const {

	return "";
}

void DummyHardwareCounter::increment(){

}
