#ifndef XMEMUSAGE_OSX_SAMPLING_THREAD_H
#define XMEMUSAGE_OSX_SAMPLING_THREAD_H

#include "xolotlMemUsage/common/AsyncSamplingThread.h"
#include "xolotlMemUsage/common/OSX/OSXSamplerBase.h"

namespace xolotlMemUsage {

using OSXSamplingThread = AsyncSamplingThread<OSX::Sample, OSX::SupportData>;

} // xolotlMemUsage

#endif // XMEMUSAGE_OSX_SAMPLING_THREAD_H

