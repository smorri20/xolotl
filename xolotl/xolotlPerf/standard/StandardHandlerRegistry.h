#ifndef STANDARDHANDLERREGISTRY_H
#define STANDARDHANDLERREGISTRY_H

#include <iostream>
#include "IHandlerRegistry.h"
#include "GPTLTimer.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLTimer
#include "GPTLHardwareCounter.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLHardwareCounter
#include "EventCounter.h"


namespace xolotlPerf {


// Factory for creating timers and hardware counter objects
// objects that use GPTL for collecting performance data,
// and event counters that interoperate with the timers and hardware counter
// objects.
//
// This is named "Standard" because it is expected to be the standard
// performance data collection mechanism used by XOLOTL.
// TODO Perhaps "Default" would be a better name?
class StandardHandlerRegistry : public IHandlerRegistry
{
private:
    static std::shared_ptr<StandardHandlerRegistry> theRegistry;

public:
    // Construct a StandardHandlerRegistry.
    StandardHandlerRegistry( void );

    // StandardHandlerRegistry is a Singleton.  (I.e., there is
    // only one instance in the process.)
    static std::shared_ptr<StandardHandlerRegistry> getRegistry( void );

    // Clean up a StandardHandlerRegistry.
    virtual ~StandardHandlerRegistry( void );

    // Obtain a Timer by name.
    virtual std::shared_ptr<ITimer> getTimer(std::string name);

    // Obtain an EventCounter by name.
    virtual std::shared_ptr<IEventCounter> getEventCounter(std::string name);

    // Obtain a HardwareCounter object by name and by the 
    // counter data it collects.
	virtual std::shared_ptr<IHardwareCounter> getHardwareCounter(
			std::string name,
			std::vector<HardwareQuantities> quantities);

	// This operation outputs the information gathered to the given 
    // output stream.
	virtual void dump(std::ostream& os) const;

};  //end class StandardHandlerRegistry

}//end namespace xolotlPerf

#endif
