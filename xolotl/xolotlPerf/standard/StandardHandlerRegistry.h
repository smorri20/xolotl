#ifndef STANDARDHANDLERREGISTRY_H
#define STANDARDHANDLERREGISTRY_H

#include <string>
#include <vector>
#include <memory>
#include "IHandlerRegistry.h"
#include "GPTLTimer.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLTimer
#include "GPTLEventCounter.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLEventCounter
#include "GPTLHardwareCounter.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLHardwareCounter

namespace xolotlPerf {

/*The StandardHandlerRegistry class realizes the interface IHandlerRegistry
 * to acquire the performance data found by implementing the performance interfaces
 * ITimer, IEventCounter, and IHardwareCounter.
 */
class StandardHandlerRegistry : public IHandlerRegistry
{

//    public:
//
//        StandardHandlerRegistry(StandardHandlerRegistry & arg);
////        StandardHandlerRegistry();
//
//        ~StandardHandlerRegistry();
//
//        // This operation returns the ITimer specified by the parameter.
//        std::shared_ptr<ITimer> getTimer(const std::string name) const;
//
//        // This operation returns the IEventCounter specified by the parameter.
//        std::shared_ptr<IEventCounter> getEventCounter(const std::string name) const;
//
//        // This operation returns the specified IHardwareCounter.
//        std::shared_ptr<IHardwareCounter> getHardwareCounter(const std::string name,
//        		const std::shared_ptr<std::vector<HardwareQuantities> > quantities) const;
////        std::shared_ptr<IHardwareCounter> getHardwareCounter(const std::string name,
////        			const std::vector<HardwareQuantities> quantities) const;
//
//        // This operation returns a list of values of the, initially specified, PAPI
//        // preset quantities monitored by the IHardwareCounter.
//        const std::shared_ptr<std::vector<HardwareQuantities> > getHardwareQuantities() const;
//
//        //        const std::vector<HardwareQuantities> & getHardwareQuantities() const;
//
//        // This operation outputs the information gathered.
//        virtual void dump(std::ostream out) const;

};  //end class StandardHandlerRegistry

}//end namespace xolotlPerf

#endif
