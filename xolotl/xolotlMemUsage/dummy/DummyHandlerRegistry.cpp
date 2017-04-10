#include "xolotlMemUsage/dummy/DummyHandlerRegistry.h"

namespace xolotlMemUsage {

/**
 * Obtain a memory usage sampler.
 */
std::shared_ptr<IMemUsageSampler> DummyHandlerRegistry::getMemUsageSampler(
        const std::string& name) {

    return std::make_shared<DummyMemUsageSampler>(name);
}

} // namespace xolotlMemUsage

