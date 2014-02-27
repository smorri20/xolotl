#include "DummyEventCounter.h"

using namespace xolotlPerf;


DummyEventCounter::DummyEventCounter(std::string eventCounterName) : IEventCounter(eventCounterName) {

	name = eventCounterName;
	value = 0;
}


DummyEventCounter::~DummyEventCounter() {

}

int DummyEventCounter::getValue() {
	return 0;
}

const std::string DummyEventCounter::getName() const {
	return "";
}

void DummyEventCounter::increment(){

}
