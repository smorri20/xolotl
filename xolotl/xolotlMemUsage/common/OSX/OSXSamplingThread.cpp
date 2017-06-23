#include <mach/host_info.h>
#include <mach/mach_host.h>
#include <mach/task_info.h>
#include <mach/task.h>
#include "xolotlMemUsage/common/OSX/OSXSamplingThread.h"

namespace xolotlMemUsage {

template<>
std::tuple<AsyncSamplingThreadBase::ClockType::time_point, OSX::Sample>
OSXSamplingThread::CollectSample(const OSX::SupportData& supportData) const
{
    // Collect the sample.
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
    vm_statistics64_data_t vmstats;
    auto callret = host_statistics64(mach_host_self(),
                                    HOST_VM_INFO64,
                                    (host_info_t)&vmstats,
                                    &count);
    if(callret != KERN_SUCCESS)
    {
        throw std::runtime_error("Unable to get memory usage info.");
    }

    OSX::Sample sample(vmstats);
    return std::make_tuple(AsyncSamplingThreadBase::ClockType::now(), sample);
}

} // namespace xolotlMemUsage

