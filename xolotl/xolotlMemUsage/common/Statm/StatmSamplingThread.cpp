#include <fstream>
#include "xolotlMemUsage/common/Statm/StatmSamplingThread.h"

namespace xolotlMemUsage {

template<>
std::tuple<AsyncSamplingThreadBase::ClockType::time_point, Statm::Sample>
StatmSamplingThread::CollectSample(const Statm::SupportData& supportData) const
{
    // Collect the sample.
    std::ifstream ifs(supportData.statmFilePath);
    uint64_t unused;
    Statm::Sample sample;
    ifs >> sample.vmSize
        >> sample.vmRSS
        >> sample.rss
        >> sample.text
        >> unused
        >> sample.dataAndStack
        >> unused;

    return std::make_tuple(AsyncSamplingThreadBase::ClockType::now(), sample);
}

} // namespace xolotlMemUsage

