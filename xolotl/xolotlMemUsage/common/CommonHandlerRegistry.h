#ifndef XMEMUSAGE_COMMON_HANDLER_REGISTRY_H
#define XMEMUSAGE_COMMON_HANDLER_REGISTRY_H

#include <map>
#include <memory>
#include "xolotlMemUsage/IHandlerRegistry.h"
#include "xolotlMemUsage/IMemUsageSampler.h"


namespace xolotlMemUsage {

/**
 * Base class for for building memory usage data collection objects that
 * collect data (as opposed to low-overhead stubs).
 */
class CommonHandlerRegistry : public IHandlerRegistry {
protected:
    /**
     * Known samplers, keyed by name.
     */
    std::map<std::string, std::shared_ptr<IMemUsageSampler> > allSamplers;

    /**
     * Construct a MemUsageSampler (of the appropriate derived type)
     * with the given name.
     *
     * @param name Name to associate with the sampler.
     */
    virtual std::shared_ptr<IMemUsageSampler> MakeMemUsageSampler(std::string name) = 0;


public:

	/**
	 * Construct a CommonHandlerRegistry.
	 */
	CommonHandlerRegistry(void) = default;

	/**
	 * Destroy the CommonHandlerRegistry.
	 */
	virtual ~CommonHandlerRegistry(void) {
        allSamplers.clear();
    }


    /**
     * Obtain a memory usage sampler.
     *
     * @param name The name of the sampler to look up or create.
     */
    virtual std::shared_ptr<IMemUsageSampler> getMemUsageSampler(const std::string& name);
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_COMMON_HANDLER_REGISTRY_H
