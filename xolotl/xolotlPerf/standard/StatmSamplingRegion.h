#ifndef MEM_SAMPLING_SUPPORT_H
#define MEM_SAMPLING_SUPPORT_H

#include <sstream>
#include <limits>
#include <unistd.h>
#include "xolotlPerf/standard/AsyncSamplingRegion.h"
#include "xolotlPerf/MemStats.h"

namespace xolotlPerf {

struct RunningPageStats
{
    uint64_t total;
    uint64_t min;
    uint64_t max;

    RunningPageStats(void)
      : total(0),
        min(std::numeric_limits<uint64_t>::max()),
        max(std::numeric_limits<uint64_t>::min())
    {
        // Nothing else to do.
    }

    void HandleNewSample(uint64_t val)
    {
        total += val;
        if(val < min)
        {
            min = val;
        }
        if(val > max)
        {
            max = val;
        }
    }

    MemPageStats GetPageStats(uint64_t nSamples) const
    {
        return MemPageStats(min, max, ((double)total) / nSamples);
    }
};

struct StatmData
{
    uint64_t nSamples;

    RunningPageStats vmSize;
    RunningPageStats vmRSS;
    RunningPageStats rss;
    RunningPageStats text;
    RunningPageStats dataAndStack;

    StatmData(void)
      : nSamples(0)
    {
        // Nothing else to do.
    }

    void HandleNewSample(uint64_t _vmSize,
                        uint64_t _vmRSS,
                        uint64_t _rss,
                        uint64_t _text,
                        uint64_t _dataAndStack)
    {
        vmSize.HandleNewSample(_vmSize);
        vmRSS.HandleNewSample(_vmRSS);
        rss.HandleNewSample(_rss);
        text.HandleNewSample(_text);
        dataAndStack.HandleNewSample(_dataAndStack);
        nSamples++;
    }

    MemStats GetCurrentStats(void) const;
};

struct StatmSupportData
{
    std::string statmFilePath;

    StatmSupportData(void)
      : statmFilePath([]{std::ostringstream s; s << "/proc/" << getpid() << "/statm"; return s.str();}())
    {
        // Nothing else to do.
    }
};

typedef AsyncSamplingRegion<StatmData, StatmSupportData> StatmSamplingRegion;

} // end namespace xolotlPerf

#endif // MEM_SAMPLING_SUPPORT_H
