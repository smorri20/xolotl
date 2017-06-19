#ifndef XMEMUSAGE_SYSINFO_SAMPLING_THREAD_H
#define XMEMUSAGE_SYSINFO_SAMPLING_THREAD_H

#include "xolotlMemUsage/common/AsyncSamplingThread.h"
#include "xolotlMemUsage/common/SysInfo/SysInfoSamplerBase.h"

namespace xolotlMemUsage {

using SysInfoSamplingThread = AsyncSamplingThread<SysInfo::Sample, SysInfo::SupportData>;

} // xolotlMemUsage

#endif // XMEMUSAGE_SYSINFO_SAMPLING_THREAD_H

