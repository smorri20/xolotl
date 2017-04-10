#ifndef PERFOBJSTATISTICS_H
#define PERFOBJSTATISTICS_H

#include <string>
#include <map>
#include <iostream>

namespace xolotlPerf {


struct PerfObjStatisticsBase {
    uint32_t processCount; //< number of processes that collected data for this metric

	/**
	 * Construct an object with default values.
	 */
    PerfObjStatisticsBase(void) :
        processCount(0) {
        // Nothing else to do.
    }

    PerfObjStatisticsBase(const PerfObjStatisticsBase& other) :
        processCount(other.processCount) {
        // Nothing else to do.
    }

    PerfObjStatisticsBase& operator=(const PerfObjStatisticsBase& other) {
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
 * Global (aggregated) statistics for a performance metric we collected
 * during a program run.
 */
template<class T>
struct PerfObjStatistics : public PerfObjStatisticsBase {

	T min; //< min value across all processes that collected data for the metric
	T max; //< max value across all processes that collected data for the metric
	double average; //< average value across all processes that collected data for the metric
	double stdev; //< standard deviation of values across all processes that collected data for the metric

	/**
	 * Construct a PerfObjStatistics struct with default values.
	 */
	PerfObjStatistics(void) :
        min(std::numeric_limits<T>::max()),
        max(std::numeric_limits<T>::min()),
        average(std::numeric_limits<T>::quiet_NaN()),
        stdev(std::numeric_limits<T>::quiet_NaN()) {
        // Nothing else to do.
	}

	/**
	 * Construct a PerfObjStatistics struct as a copy of another.
	 * @param obj The object to be copied.
	 */
	PerfObjStatistics(const PerfObjStatistics& obj) :
        PerfObjStatisticsBase(obj),
        min(obj.min),
        max(obj.max),
        average(obj.average),
        stdev(obj.stdev) {
        // Nothing else to do.
	}

	/**
	 * Replace my own values to be a copy of another PerfObjStatistics.
	 * @param obj The object to be copied.
	 */
	PerfObjStatistics& operator=(const PerfObjStatistics& obj) {
		if (&obj != this) {
            PerfObjStatisticsBase::operator=(obj);
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
        PerfObjStatisticsBase::outputTo(os);
		os << "    " << "min: " << min << '\n'
            << "    " << "max: " << max << '\n' 
            << "    " << "average: " << average << '\n' 
            << "    " << "stdev: " << stdev << '\n';
	}
};


/* C++11 supports template typedefs.   Earlier C++ standards did not,
 * but we are requiring C++11 for the shared_ptr support so we know
 * that we are using a C++11 compiler.
 */
template<typename T>
using PerfObjStatsMap = std::map<std::string, T>;

} // namespace xolotlPerf

#endif // PERFOBJSTATISTICS_H
