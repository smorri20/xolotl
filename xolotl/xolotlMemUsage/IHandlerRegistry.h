#ifndef XMEMUSAGE_IHANDLERREGISTRY_H
#define XMEMUSAGE_IHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include "IMemUsageSampler.h"
#include "xolotlMemUsage/common/AsyncSamplingThread.h"


namespace xolotlMemUsage {

/**
 * Factory for building memory usage data collection objects.
 */
class IHandlerRegistry {

public:

    /**
     * Globally aggregated memory usage data for all
     * data types we know how to collect.
     * This is a base class for derived classes to extend.
     */
    struct GlobalData {
        
        virtual ~GlobalData(void) { }   // required to make GlobalData polymorphic
    };


	/// Possible types of memory usage handler registries.
	enum RegistryType {
		dummy,      //< Use stub classes that do not collect any data.
		std,        //< Collect memory usage statistics.
        profile     //< Collect memory usage profiles.
	};


    /**
     * Type of values used to specify asynchronous sampling interval.
     */
    using SamplingInterval = AsyncSamplingThreadBase::ClockType::duration;


	/**
	 * The destructor
	 */
	virtual ~IHandlerRegistry() {
	}


    /**
     * Obtain a memory usage sampler with the given name.
     */
    virtual std::shared_ptr<IMemUsageSampler> getMemUsageSampler(
            const std::string& name) = 0;


	/**
	 * Collect statistics about any memory usage data collected by
	 * processes of the program.
	 */
    virtual std::shared_ptr<GlobalData> collectData(void) const = 0;

	/**
	 * Report memory usage data statistics to the given stream.
	 *
	 * @param os Stream on which to output statistics.
     * @param stats The memory usage statistics to display.
	 *
	 */
	virtual void reportData(std::ostream& os,
                                std::shared_ptr<GlobalData> data) const = 0;
};

} //end namespace xolotlMemUsage

#endif // XMEMUSAGE_IHANDLERREGISTRY_H
