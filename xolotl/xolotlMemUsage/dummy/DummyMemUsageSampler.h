#ifndef DUMMY_MEM_USAGE_SAMPLER_H
#define DUMMY_MEM_USAGE_SAMPLER_H

#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlCore/Identifiable.h"


namespace xolotlMemUsage {

/**
 * A "dummy" memory usage sampler.
 * Implements the interface, so that the caller's code doesn't
 * have to change, but does nothing.
 */
class DummyMemUsageSampler : public IMemUsageSampler, xolotlCore::Identifiable {

private:
    /**
     * Construct a DummyMemUsageSampler.
     * Declared private to enforce that MemUsageSamplers must
     * be constructed with a name.
     */
    DummyMemUsageSampler(void)
      : xolotlCore::Identifiable("unused") {
        // Nothing else to do.
    }
    
public:

	/**
     * Construct a DummyMemUsageSampler with the given name.
	 *
	 * @param name The object's name.
	 */
	DummyMemUsageSampler(const std::string& name)
      : xolotlCore::Identifiable("unused") {
	}

	/**
	 * Destroy the MemUsageSampler.
	 */
	virtual ~DummyMemUsageSampler(void) {
        // Nothing to do.
	}

	/**
	 * Start sampling.
	 */
	virtual void start(void) {
        // Nothing to do.
    }

	/**
	 * Stop sampling.
	 */
	virtual void stop(void) {
        // Nothing to do.
    }

	/**
	 * Obtain the sampler's collected metric values.
	 */
	virtual IMemUsageSampler::ValType getValue(void) const {
        return IMemUsageSampler::ValType();
    }
};
//end class DummyMemUsageSampler

}//end namespace xolotlMemUsage

#endif // DUMMY_MEM_USAGE_SAMPLER_H
