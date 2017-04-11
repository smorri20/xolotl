#include <fstream>
#include "xolotlMemUsage/standard/StatmSampler.h"

namespace xolotlMemUsage {

template<>
StatmSampler::SamplingThreadType::SamplingIntervalType StatmSampler::samplingInterval = std::chrono::duration<double>(1);


template<>
void
StatmSampler::SamplingThreadType::Sample(void)
{
    // Collect the sample.
    std::ifstream ifs(supportData.statmFilePath);
    uint64_t vmSize;
    uint64_t vmRSS;
    uint64_t rss;
    uint64_t text;
    uint64_t libUnused;
    uint64_t dataAndStack;
    uint64_t dirtyUnused;

    ifs >> vmSize
        >> vmRSS
        >> rss
        >> text
        >> libUnused
        >> dataAndStack
        >> dirtyUnused;

    // Ensure we are the only ones updating our data structures.
    std::lock_guard<std::mutex> ourLock(mtx);

    // Add the sample to the data for all active samplers.
    for(auto currSampler : activeSamplers)
    {
        currSampler->GetRunningSampleData().HandleNewSample(vmSize,
                                                            vmRSS,
                                                            rss,
                                                            text,
                                                            dataAndStack);
    }
}


std::shared_ptr<MemUsageStats>
StatmData::GetCurrentStats(void) const
{
    return std::make_shared<MemUsageStats>(vmSize.GetPageStats(nSamples),
                                vmRSS.GetPageStats(nSamples),
                                rss.GetPageStats(nSamples),
                                text.GetPageStats(nSamples),
                                dataAndStack.GetPageStats(nSamples));
}


} // end namespace xolotlMemUsage

