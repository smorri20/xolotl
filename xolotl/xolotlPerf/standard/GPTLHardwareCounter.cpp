#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <papi.h>
#include "GPTLHardwareCounter.h"

using namespace xolotlPerf;

// Create the static map of hardware quantities
//std::unordered_map<HardwareQuantities, int> GPTLHardwareCounter::hardwareQuantitiesMap;

GPTLHardwareCounter::GPTLHardwareCounter(std::string counterName,
		const std::vector<HardwareQuantities> &counterQuantities) :
		IHardwareCounter(counterName, counterQuantities) {

	name = counterName;
	quantities = counterQuantities;

	// Set up the hardware quantities map
//	hardwareQuantitiesMap = { {L1_CACHE_MISS, PAPI_L1_TCM},
//			{L2_CACHE_MISS, PAPI_L2_TCM}, {L3_CACHE_MISS, PAPI_L3_TCM},
//			{BRANCH_MISPRED, PAPI_BR_MSP}, {TOTAL_CYCLES, PAPI_TOT_CYC},
//			{TOTAL_INSTRUC, PAPI_TOT_INS}, {FLPT_INSTRUC, PAPI_FP_INS} };

}

GPTLHardwareCounter::~GPTLHardwareCounter() {

}

//std::vector<int> GPTLHardwareCounter::getValues() const {

	// The following documentation was taken directly from gptl.c
	/*
	 ** GPTLget_eventvalue: return PAPI-based event value for a timer. All values will be
	 ** returned as doubles, even if the event is not derived.
	 **
	 ** Input args:
	 ** const char *timername: timer name
	 ** const char *eventname: event name (must be currently enabled)
	 ** int t: thread number (if < 0, the request is for the current thread)
	 **
	 ** Output args:
	 ** double *value: current value of the event for this timer
	 */
	//	int gret = GPTLget_eventvalue( name.c_str(), -1, &papival );

//	return values;
//}

const std::string GPTLHardwareCounter::getName() const {

	return name;
}

std::vector<HardwareQuantities> GPTLHardwareCounter::getHardwareQuantities() const {

	// TO DO:  return the quantities as vector of strings

	return quantities;
}


