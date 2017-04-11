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
            const std::string& name) {

        return std::make_shared<DummyMemUsageSampler>(name);
    }


	/**
	 * Collect data about any memory usage collected by
	 * any processes of the program.
	 * This method is a stub.
	 */
	virtual std::shared_ptr<GlobalData> collectData(void) const {

        return std::make_shared<GlobalData>();
    }

	/**
	 * Report memory usage data to the given stream.
	 * This method is a stub in this class.
	 *
	 * @param os Stream on which to output data.
     * @param stats Collected data.
	 */
	virtual void reportData(std::ostream& os,
                                    std::shared_ptr<GlobalData> data) const {

        // Nothing to do.
    }
};

} //end namespace xolotlMemUsage

#endif
