#include <sys/sysinfo.h>
#include "xolotlMemUsage/common/SysInfo/SysInfoSamplingThread.h"

namespace xolotlMemUsage {

template<>
std::tuple<AsyncSamplingThreadBase::ClockType::time_point, SysInfo::Sample>
SysInfoSamplingThread::CollectSample(const SysInfo::SupportData& supportData) const
{
    // Collect the sample.
    struct sysinfo si;
    sysinfo(&si);

    SysInfo::Sample sample(si);
    return std::make_tuple(AsyncSamplingThreadBase::ClockType::now(), sample);
}

} // namespace xolotlMemUsage

