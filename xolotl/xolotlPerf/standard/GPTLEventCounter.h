#ifndef GPTLEVENTCOUNTER_H
#define GPTLEVENTCOUNTER_H
#include "gptl.h"
#include "papi.h"
#include "IEventCounter.h"

namespace xolotlPerf{

/**
 * The GPTLEventCounter class is instantiated by the StandardHandlerRegistry
 * class and realizes the IEventCounter interface to access event performance
 * counter data found via the General Purpose Timing Library (GPTL).
 */
class GPTLEventCounter : public IEventCounter
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
	GPTLEventCounter():IEventCounter(""), name("private"), value(0) { }


public:

	/**
	 * GPTLEventCounter constructor that takes the argument name
	 *
	 * @param eventCounterName The GPTLEventCounter's name
	 */
	GPTLEventCounter(std::string eventCounterName);

	/**
	 * The destructor
	 */
	~GPTLEventCounter();

	/**
	 * This operation returns the value of the GPTLEventCounter,
	 * the frequency of the specified event.
	 */
	int getValue();

	/**
	 * This operation returns the name of the GPTLEventCounter.
	 *
	 * @return The name of this GPTLEventCounter
	 */
	const std::string getName() const;

	/**
	 * This operation increments the GPTLEventCounter.
	 */
	void increment();

};  //end class GPTLEventCounter


}  //end namespace xolotlPerf

#endif
