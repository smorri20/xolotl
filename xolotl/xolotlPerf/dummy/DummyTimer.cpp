#include "DummyTimer.h"

using namespace xolotlPerf;

DummyTimer::DummyTimer(std::string timerName) : ITimer(timerName) {

}

DummyTimer::~DummyTimer() {

}

void DummyTimer::start() {

}

void DummyTimer::stop() {

}

const std::string DummyTimer::getName() const {
	return "";
}

double DummyTimer::getValue() {
	return 0;
}

long DummyTimer::getUnits() const {
	return 0;
}
