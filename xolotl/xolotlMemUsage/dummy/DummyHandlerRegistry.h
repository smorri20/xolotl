#ifndef XMEMUSAGE_DUMMY_HANDLER_REGISTRY_H
#define XMEMUSAGE_DUMMY_HANDLER_REGISTRY_H

#include "xolotlMemUsage/IHandlerRegistry.h"
#include "xolotlMemUsage/dummy/DummyMemUsageSampler.h"


namespace xolotlMemUsage {

// Facgtory for creating memory usage tracking objects that are
// dummies, i.e., they provide the right interface but don't do
// anything.  They are stubs.  This is so that the client code can be 
// written to use the memory usage data collection infrastructure without 
// regard to whether data collection is active or disabled.
//
class DummyHandlerRegistry: public IHandlerRegistry {
public:
	DummyHandlerRegistry(IHandlerRegistry::SamplingInterval /* samplingInterval */) {
	}

	virtual ~DummyHandlerRegistry(void) {
	}

    /**
     * Obtain a memory usage sampler.
     */
    virtual std::shared_ptr<IMemUsageSampler> getMemUsageSampler(
            const std::string& name);

	/**
	 * Collect statistics about any memory usage data collected by
	 * processes of the program.
	 * This method is a stub.
	 */
	virtual GlobalMemUsageStats collectStatistics(void) const {
        return GlobalMemUsageStats();
    }

	/**
	 * Report memory usage data statistics to the given stream.
	 * This method is a stub in this class.
	 *
	 * @param os Stream on which to output statistics.
     * @param stats Collected statistics.
	 */
	virtual void reportStatistics(std::ostream& os,
                                    const GlobalMemUsageStats& stats) const {

        // Nothing to do.
    }
};

} //end namespace xolotlMemUsage

#endif
