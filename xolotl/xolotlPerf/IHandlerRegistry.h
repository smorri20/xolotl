#ifndef IHANDLERREGISTRY_H
#define IHANDLERREGISTRY_H

// Includes
#include <string>
#include <vector>
#include <memory>
#include "ITimer.h"
#include "IEventCounter.h"
#include "IHardwareCounter.h"
#include "IMemSamplingRegion.h"
#include "PerfObjStatistics.h"


namespace xolotlPerf {

/**
 * Factory for building performance data collection objects, such
 * as timers and counters.
 */
class IHandlerRegistry {

public:

    /**
     * Globally aggregated performance statistics for all
     * performance data types we know how to collect.
     */
    struct GlobalPerfStats {
        PerfObjStatsMap<ITimer::GlobalStatsType> timerStats;
        PerfObjStatsMap<IEventCounter::GlobalStatsType> counterStats;
        PerfObjStatsMap<IHardwareCounter::GlobalStatsType> hwCounterStats;
        PerfObjStatsMap<IMemSamplingRegion::GlobalStatsType> memStats;
    };


	/// Possible types of performance handler registries.
	enum RegistryType {
		dummy,     //< Use stub classes that do not collect any performance data
		std,        //< Use the best available API.
		os,         //< Use operating system/runtime API.
		papi,       //< Use PAPI to collect performance data.
	};

	/**
	 * The destructor
	 */
	virtual ~IHandlerRegistry() {
	}

	/**
	 * This operation returns the ITimer specified by the parameter.
	 */
	virtual std::shared_ptr<ITimer> getTimer(const std::string& name) = 0;

	/**
	 * This operation returns the IEventCounter specified by the parameter.
	 */
	virtual std::shared_ptr<IEventCounter> getEventCounter(
			const std::string& name) = 0;

	/**
	 * This operation returns the specified IHardwareCounter.
	 */
	virtual std::shared_ptr<IHardwareCounter> getHardwareCounter(
			const std::string& name,
			const IHardwareCounter::SpecType& ctrSpec) = 0;

    /**
     * Obtain a memory sampling region.
     */
    virtual std::shared_ptr<IMemSamplingRegion> getMemSamplingRegion(
            const std::string& name) = 0;


	/**
	 * Collect statistics about any performance data collected by
	 * processes of the program.
	 */
    virtual GlobalPerfStats collectStatistics(void) const = 0;

	/**
	 * Report performance data statistics to the given stream.
	 *
	 * @param os Stream on which to output statistics.
     * @param stats The performance statistics to display.
	 *
	 */
	virtual void reportStatistics(std::ostream& os, 
                                const GlobalPerfStats& stats) const = 0;
};

} //end namespace xolotlPerf

#endif
