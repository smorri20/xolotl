#ifndef XMEMUSAGE_IHANDLERREGISTRY_H
#define XMEMUSAGE_IHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include "IMemUsageSampler.h"
#include "MemUsageObjStatistics.h"


namespace xolotlMemUsage {

/**
 * Factory for building memory usage data collection objects.
 */
class IHandlerRegistry {

public:

    /**
     * Globally aggregated memory usage statistics for all
     * data types we know how to collect.
     */
    struct GlobalMemUsageStats {
        MemUsageObjStatsMap<IMemUsageSampler::GlobalStatsType> memStats;
    };


	/// Possible types of memory usage handler registries.
	enum RegistryType {
		dummy,     //< Use stub classes that do not collect any data
		std        //< Use the default API.
	};


    /**
     * Type of values used to specify asynchronous sampling interval.
     */
    typedef std::chrono::duration<double> SamplingInterval;


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
    virtual GlobalMemUsageStats collectStatistics(void) const = 0;

	/**
	 * Report memory usage data statistics to the given stream.
	 *
	 * @param os Stream on which to output statistics.
     * @param stats The memory usage statistics to display.
	 *
	 */
	virtual void reportStatistics(std::ostream& os, 
                                const GlobalMemUsageStats& stats) const = 0;
};

} //end namespace xolotlMemUsage

#endif // XMEMUSAGE_IHANDLERREGISTRY_H
