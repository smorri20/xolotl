#ifndef DUMMYEVENTCOUNTER_H
#define DUMMYEVENTCOUNTER_H

#include <string>
#include <memory>
#include "IEventCounter.h"

using namespace std;

namespace xolotlPerf{

// The DummyEventCounter class is instantiated by the DummyHandlerRegistry
// class and realizes the DummyEventCounter interface.
class DummyEventCounter : public IEventCounter
{

private:

	/**
	 * The name of this IEventCounter.
	 */
	std::string name;

	/**
	 * The value of this IEventCounter.
	 */
	int value;

	/**
	 * The default constructor is declared private since all EventCounters
	 *  must be initialized with a name.
	 */
	DummyEventCounter():IEventCounter(""), name("private"), value(0) {}


public:

	/**
	 * DummyEventCounter constructor that takes the argument name
	 *
	 * @param eventCounterName The DummyEventCounter's name
	 */
	DummyEventCounter(std::string eventCounterName);

	/**
	 * The destructor
	 */
	~DummyEventCounter();

	/**
	 * This operation returns the value of the DummyEventCounter,
	 * the frequency of the specified event.
	 */
	int getValue();

	/**
	 * This operation returns the name of the DummyEventCounter.
	 *
	 * @return The name of this DummyEventCounter
	 */
	const std::string getName() const;

	/**
	 * This operation increments the DummyEventCounter.
	 */
	void increment();


};  //end class DummyEventCounter

}  //end namespace xolotlPerf

#endif
