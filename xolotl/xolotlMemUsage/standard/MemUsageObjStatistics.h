#ifndef MEM_USAGE_OBJ_STATISTICS_H
#define MEM_USAGE_OBJ_STATISTICS_H

#include <string>
#include <map>
#include <iostream>
#include "xolotlMemUsage/standard/MemUsageStats.h"

namespace xolotlMemUsage {


struct MemUsageObjStatisticsBase {
    uint32_t processCount; //< number of processes that collected data for this metric

	/**
	 * Construct an object with default values.
	 */
    MemUsageObjStatisticsBase(void) :
        processCount(0) {
        // Nothing else to do.
    }

    MemUsageObjStatisticsBase(const MemUsageObjStatisticsBase& other) :
        processCount(other.processCount) {
        // Nothing else to do.
    }

    MemUsageObjStatisticsBase& operator=(const MemUsageObjStatisticsBase& other) {
        if(&other != this) {
            processCount = other.processCount;
        }
        return *this;
    }
    
	/**
	 * Output our name and statistics to the given output stream.
	 * @param os The output stream on which we will write our statistics.
	 */
	void outputTo(std::ostream& os) const {
		os << "    " << "process_count: " << processCount << '\n';
    }
};

/**
 * Global (aggregated) statistics for a metric we collected
 * during a program run.
 */
template<class T>
struct MemUsageObjStatistics : public MemUsageObjStatisticsBase {

	T min; //< min value across all processes that collected data for the metric
	T max; //< max value across all processes that collected data for the metric
	double average; //< average value across all processes that collected data for the metric
	double stdev; //< standard deviation of values across all processes that collected data for the metric

	/**
	 * Construct a MemUsageObjStatistics struct with default values.
	 */
	MemUsageObjStatistics(void) :
        min(std::numeric_limits<T>::max()),
        max(std::numeric_limits<T>::min()),
        average(std::numeric_limits<T>::quiet_NaN()),
        stdev(std::numeric_limits<T>::quiet_NaN()) {
        // Nothing else to do.
	}

	/**
	 * Construct a MemUsageObjStatistics struct as a copy of another.
	 * @param obj The object to be copied.
	 */
	MemUsageObjStatistics(const MemUsageObjStatistics& obj) :
        MemUsageObjStatisticsBase(obj),
        min(obj.min),
        max(obj.max),
        average(obj.average),
        stdev(obj.stdev) {
        // Nothing else to do.
	}

	/**
	 * Replace my own values to be a copy of another MemUsageObjStatistics.
	 * @param obj The object to be copied.
	 */
	MemUsageObjStatistics& operator=(const MemUsageObjStatistics& obj) {
		if (&obj != this) {
            MemUsageObjStatisticsBase::operator=(obj);
			min = obj.min;
			max = obj.max;
			average = obj.average;
			stdev = obj.stdev;
		}
		return *this;
	}

	/**
	 * Output our name and statistics to the given output stream.
	 * @param os The output stream on which we will write our statistics.
	 */
	void outputTo(std::ostream& os) const {
        MemUsageObjStatisticsBase::outputTo(os);
		os << "    " << "min: " << min << '\n'
            << "    " << "max: " << max << '\n' 
            << "    " << "average: " << average << '\n' 
            << "    " << "stdev: " << stdev << '\n';
	}
};


// Specialization for memory statistics.
// The default definition is suitable for scalars, but not for 
// structures/classes.
template<>
struct MemUsageObjStatistics<MemUsageStats> : public MemUsageObjStatisticsBase {

    MemUsageStats stats; //< statistics across all processes


	/**
	 * Construct a MemUsageObjStatistics struct with default values.
	 * @param _name The metric name.
	 */
	MemUsageObjStatistics(void) {
        // Nothing else to do.
	}

	/**
	 * Construct a MemUsageObjStatistics struct as a copy of another.
	 * @param obj The object to be copied.
	 */
	MemUsageObjStatistics(const MemUsageObjStatistics& other) :
        MemUsageObjStatisticsBase(other),
        stats(other.stats) {
        // Nothing else to do.
	}

	/**
	 * Replace my own values to be a copy of another MemUsageObjStatistics.
	 * @param obj The object to be copied.
	 */
	MemUsageObjStatistics& operator=(const MemUsageObjStatistics& other) {
		if (&other != this) {
            MemUsageObjStatisticsBase::operator=(other);
            stats = other.stats;
		}
		return *this;
	}

	/**
	 * Output our name and statistics to the given output stream.
	 * @param os The output stream on which we will write our statistics.
	 */
	void outputTo(std::ostream& os) const {
        MemUsageObjStatisticsBase::outputTo(os);
        stats.OutputTo(os, "    ");
	}
};


/* C++11 supports template typedefs.   Earlier C++ standards did not,
 * but we are requiring C++11 for the shared_ptr support so we know
 * that we are using a C++11 compiler.
 */
template<typename T>
using MemUsageObjStatsMap = std::map<std::string, T>;

} // namespace xolotlMemUsage

#endif // MEM_USAGE_OBJ_STATISTICS_H
