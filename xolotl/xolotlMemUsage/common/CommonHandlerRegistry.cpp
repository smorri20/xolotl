#include "xolotlMemUsage/common/CommonHandlerRegistry.h"

namespace xolotlMemUsage {

std::shared_ptr<IMemUsageSampler> CommonHandlerRegistry::getMemUsageSampler(const std::string& name) {

    std::shared_ptr<IMemUsageSampler> ret;

    // Have we already created a sampler with this name?
    auto iter = allSamplers.find(name);
    if(iter != allSamplers.end()) {

        // We have already created a memory usage sampler with this name,
        // so return that one.
        ret = iter->second;
    }
    else {
        // We have not already created a memory usage sampler with this name.
        // Create one, keep track of it, and return it.
        ret = MakeMemUsageSampler(name);
        allSamplers.emplace(name, ret);
    }

    return ret;
}

} // xolotlMemUsage

