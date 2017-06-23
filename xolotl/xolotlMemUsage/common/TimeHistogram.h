#ifndef XMEMUSAGE_TIME_HISTOGRAM_H
#define XMEMUSAGE_TIME_HISTOGRAM_H


#include <chrono>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <cmath>



namespace xolotlMemUsage {


template<typename SType,
            typename RType,
            typename MType,
            typename ClockType = std::chrono::system_clock>
class TimeHistogram
{
public:
    using TSType = typename ClockType::time_point;
    using DType = typename ClockType::duration;

private:
    struct BinData
    {
        RType runningValue;
        uint64_t nSamples;

        BinData(void)
        {
            // We do not set our initial values directly
            // in the constructor so that a user
            // can provide their own initial values,
            // and so that we don't force them to 
            // implement a particular type of constructor.
            Reset();
        }

        void Reset(void);

        void SetToCombinationOf(typename std::vector<BinData>::const_iterator first,
                                typename std::vector<BinData>::const_iterator last);

        uint64_t GetNumSamples(void) const  { return nSamples; }
        const RType& GetRunningData(void) const    { return runningValue; }
        MType GetMetricValue(void) const;

        void HandleSample(const SType& s);
    };

    unsigned int nBins;
    DType binWidth;
    TSType startTimestamp;
    std::vector<BinData> bins;


    void Fold(void);
    unsigned int FindBin(const TSType& timestamp);

public:
    TimeHistogram(unsigned int _nBins,
                    DType _binWidth,
                    TSType _startTimestamp);

    void HandleSample(const TSType& timestamp,
                        const SType& sample);

    const std::vector<BinData>& GetBins(void) const { return bins; }

    TSType GetStartTimestamp(void) const { return startTimestamp; }
    DType GetBinWidth(void) const   { return binWidth; }
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_TIME_HISTOGRAM_H
