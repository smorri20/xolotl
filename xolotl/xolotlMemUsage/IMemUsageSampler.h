#ifndef XMEMUSAGE_IMEMUSAGE_SAMPLER_H
#define XMEMUSAGE_IMEMUSAGE_SAMPLER_H

#include <memory>
#include "xolotlCore/IIdentifiable.h"

namespace xolotlMemUsage {

class IMemUsageSampler : public virtual xolotlCore::IIdentifiable {

public:
    /**
     * Memory usage data that we can collect.
     * This is a base class for derived classes to extend.
     */
    struct MemUsageData {
        virtual ~MemUsageData(void) { } // required to make class polymorphic
    };

    /// Destroy the memory usage sampler.
    virtual ~IMemUsageSampler(void) { }

    /// Start collecting memory usage samples.
    virtual void start(void) = 0;

    /// Stop collecting memory usage samples.
    virtual void stop(void) = 0;

    /// Obtain the sampler's current values.
    virtual std::shared_ptr<MemUsageData> getValue(void) const = 0;
};

}//end namespace xolotlMemUsage

#endif // XMEMUSAGE_IMEMUSAGE_SAMPLER_H
