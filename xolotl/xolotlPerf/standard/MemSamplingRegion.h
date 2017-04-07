#ifndef MEM_SAMPLING_REGION_H
#define MEM_SAMPLING_REGION_H

#include "xolotlCore/Identifiable.h"
#include "xolotlPerf/IMemSamplingRegion.h"
#include "xolotlPerf/standard/StatmSamplingRegion.h"


namespace xolotlPerf {

class MemSamplingRegion : public IMemSamplingRegion,
                            public StatmSamplingRegion,
                            public xolotlCore::Identifiable {

private:
    /**
     * Construct a MemSamplingRegion.
     * Declared private to enforce that MemSamplingRegions must
     * be constructed with a name.
     */
    MemSamplingRegion(void)
      : xolotlCore::Identifiable("unused"),
        StatmSamplingRegion("unused") {
        // Nothing else to do.
    }
    
public:

	/**
     * Construct a MemSamplingRegion with the given name.
	 *
	 * @param name The object's name.
	 */
	MemSamplingRegion(const std::string& name)
      : xolotlCore::Identifiable(name),
        StatmSamplingRegion(name) {

        // Nothing else to do.
    }

	/**
	 * Destroy the MemSamplingRegion.
	 */
	virtual ~MemSamplingRegion(void) {

        // Ensure we've stopped sampling.
        stop();
    }

	/**
	 * Start sampling for this memory region.
	 */
	virtual void start(void) {

        StartSampling();
    }

	/**
	 * Stop sampling for this memory region.
	 */
	virtual void stop(void) {

        StopSampling();
    }

	/**
	 * Obtain the region's metric values.
	 */
	virtual IMemSamplingRegion::ValType getValue(void) const {

        return GetRunningSampleData().GetCurrentStats();
    }
};

} // end namespace xolotlPerf

#endif // MEM_SAMPLING_REGION_H
