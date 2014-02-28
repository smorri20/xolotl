#ifndef DUMMYEVENTCOUNTER_H
#define DUMMYEVENTCOUNTER_H

#include <string>
#include "Identifiable.h"
#include "IEventCounter.h"

using namespace std;

namespace xolotlPerf{

// The DummyEventCounter class is instantiated by the DummyHandlerRegistry
// class and realizes the DummyEventCounter interface.
class DummyEventCounter : public IEventCounter, public xolotlCore::Identifiable
{

private:

	/**
	 * The value of this IEventCounter.
	 */
	unsigned long value;

	/**
	 * The default constructor is declared private since all EventCounters
	 *  must be initialized with a name.
	 */
    DummyEventCounter(void)
      : xolotlCore::Identifiable("unused"),
        value( 0 )
    { }


public:

	/**
	 * DummyEventCounter constructor that takes the argument name
	 *
	 * @param name The DummyEventCounter's name
	 */
	DummyEventCounter(std::string name)
      : xolotlCore::Identifiable( name ),
        value( 0 )
    { }

	/**
	 * The destructor
	 */
	virtual ~DummyEventCounter() { }

	/**
	 * This operation returns the value of the DummyEventCounter,
	 * the frequency of the specified event.
	 */
	virtual unsigned long getValue() const  { return value; }


	/**
	 * This operation increments the DummyEventCounter.
	 */
	virtual void increment()    { ++value; }


};  //end class DummyEventCounter

}  //end namespace xolotlPerf

#endif
