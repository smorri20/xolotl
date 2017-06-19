#ifndef XMEMUSAGE_STATM_SAMPLING_THREAD_H
#define XMEMUSAGE_STATM_SAMPLING_THREAD_H

#include "xolotlMemUsage/common/AsyncSamplingThread.h"
#include "xolotlMemUsage/common/Statm/StatmSamplerBase.h"

namespace xolotlMemUsage {

using StatmSamplingThread = AsyncSamplingThread<Statm::Sample, Statm::SupportData>;

} // xolotlMemUsage

#endif // XMEMUSAGE_STATM_SAMPLING_THREAD_H

