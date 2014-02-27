#ifndef IEVENTCOUNTER_H
#define IEVENTCOUNTER_H

#include <string>

using namespace std;

namespace xolotlPerf {

//Realizations of this interface are responsible for the collection
//of event performance counter data.
class IEventCounter {

//private:

	/**
	 * The default constructor is declared private since all event counters
	 * must be initialized with a name.
	 */
//	IEventCounter();

public:

	/**
	 * IEventCounter constructor that takes the argument name
	 *
	 * @param eventCounterName The IEventCounter's name
	 */
	IEventCounter(std::string eventCounterName) {}

	/**
	 * The destructor
	 */
	virtual ~IEventCounter(){}

	/**
	 * This operation returns the value of the IEventCounter, the frequency
	 * of the specified event.
	 */
	virtual int getValue() = 0;

	/**
	 * This operation returns the name of the IEventCounter.
	 *
	 * @return The name of this IEventCounter
	 */
	virtual const std::string getName() const = 0;

	/**
	 * This operation increments the IEventCounter.
	 */
	virtual void increment() = 0;

};
//end class IEventCounter

}//end namespace xolotlPerf

#endif
