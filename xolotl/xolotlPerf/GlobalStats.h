#ifndef XOLOTL_PERF_GLOBAL_STATS_H
#define XOLOTL_PERF_GLOBAL_STATS_H

namespace xolotlPerf {


/**
 * Type of global statisics, suitable for a scalar.
 */
template<typename T>
struct GlobalStats {

    T min;
    T max;
    double average;
    double stdev;
};


}; // namespace xolotlPerf

#endif // XOLOTL_PERF_GLOBAL_STATS_H
