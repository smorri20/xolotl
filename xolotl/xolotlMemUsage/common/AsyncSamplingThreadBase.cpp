#include "xolotlMemUsage/common/AsyncSamplingThread.h"

namespace xolotlMemUsage {

// Default to sampling every second.
AsyncSamplingThreadBase::ClockType::duration AsyncSamplingThreadBase::samplingInterval = std::chrono::duration<uint64_t, std::milli>(1000);

} // namespace xolotlMemUsage

