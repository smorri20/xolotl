#ifndef DUMMY_MEM_SAMPLING_REGION_H
#define DUMMY_MEM_SAMPLING_REGION_H

#include "xolotlPerf/IMemSamplingRegion.h"
#include "xolotlCore/Identifiable.h"


namespace xolotlPerf {

/**
 * A "dummy" memory sampling region.
 * Implements the interface, so that the caller's code doesn't
 * have to change, but does nothing.
 */
class DummyMemSamplingRegion : public IMemSamplingRegion, xolotlCore::Identifiable {

private:
    /**
     * Construct a DummyMemSamplingRegion.
     * Declared private to enforce that MemSamplingRegions must
     * be constructed with a name.
     */
    DummyMemSamplingRegion(void)
      : xolotlCore::Identifiable("unused") {
        // Nothing else to do.
    }
    
public:

	/**
     * Construct a DummyMemSamplingRegion with the given name.
	 *
	 * @param name The object's name.
	 */
	DummyMemSamplingRegion(const std::string& name)
      : xolotlCore::Identifiable("unused") {
	}

	/**
	 * Destroy the MemSamplingRegion.
	 */
	virtual ~DummyMemSamplingRegion(void) {
        // Nothing to do.
	}

	/**
	 * Start sampling for the region.
	 */
	virtual void start(void) {
        // Nothing to do.
    }

	/**
	 * Stop sampling for the region.
	 */
	virtual void stop(void) {
        // Nothing to do.
    }

	/**
	 * Obtain the region's metric values.
	 */
	virtual IMemSamplingRegion::ValType getValue(void) const {
        return IMemSamplingRegion::ValType();
    }
};
//end class DummyMemSamplingRegion

}//end namespace xolotlPerf

#endif // DUMMY_MEM_SAMPLING_REGION_H
