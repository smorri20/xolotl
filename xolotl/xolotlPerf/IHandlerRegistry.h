#ifndef IHANDLERREGISTRY_H
#define IHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include "ITimer.h" //Dependency Generated Source:IHandlerRegistry Target:ITimer
#include "IEventCounter.h" //Dependency Generated Source:IHandlerRegistry Target:IEventCounter
#include "IHardwareCounter.h" //Dependency Generated Source:IHandlerRegistry Target:IHardwareCounter
#include "HardwareQuantities.h"

namespace xolotlPerf {

/**
 * Realizations of this interface are responsible for the collection of performance data.
 */
class IHandlerRegistry {

public:

	// The destructor
	virtual ~IHandlerRegistry(){}

	// This operation returns the ITimer specified by the parameter.
	virtual std::shared_ptr<ITimer> getTimer(std::string name) const = 0;

	// This operation returns the IEventCounter specified by the parameter.
	virtual std::shared_ptr<IEventCounter> getEventCounter(
			std::string name) const = 0;

	// This operation returns the specified IHardwareCounter.
	virtual std::shared_ptr<IHardwareCounter> getHardwareCounter(
			std::string name,
			std::vector<HardwareQuantities> quantities) const = 0;

	// This operation outputs the information gathered.
	virtual void dump(std::ostream out) const = 0;

};
//end class IHandlerRegistry

}//end namespace xolotlPerff

#endif
