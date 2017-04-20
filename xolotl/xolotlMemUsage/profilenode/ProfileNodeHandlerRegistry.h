#ifndef XMEMUSAGE_PROFILE_NODE_HANDLER_REGISTRY_H
#define XMEMUSAGE_PROFILE_NODE_HANDLER_REGISTRY_H

#include <map>
#include <memory>
#include "xolotlMemUsage/common/CommonHandlerRegistry.h"
#include "xolotlMemUsage/profilenode/NodeMemUsageProfile.h"
#include "xolotlMemUsage/profilenode/NodeMemUsageProfiler.h"


namespace xolotlMemUsage {

/**
 * Registry for building memory usage profilers.
 */
class ProfileNodeHandlerRegistry : public CommonHandlerRegistry {
public:

    /**
     * Global mem usage profile data for all known profilers.
     */
    struct GlobalMemUsageData : public IHandlerRegistry::GlobalData {
        std::map<std::string, std::tuple<uint32_t, NodeMemUsageProfile> > profiles;
    };


private:
    /// Name of file to which memory usage profile data should be saved.
    std::string oFilename;  

    /// Communicator to use for aggregating results.
    MPI_Comm aggComm;

    /// Rank within aggregating communicator.
    int aggCommRank;
    

    /**
     * Construct a profiler (of the appropriate derived type)
     * with the given name.
     *
     * @param name Name to associate with the sampler.
     */
    virtual std::shared_ptr<IMemUsageSampler> MakeMemUsageSampler(std::string name) {
        return std::make_shared<NodeMemUsageProfiler>(name);
    }

public:
	/**
	 * Construct a Registry for collecting memory usage profiles.
	 */
    ProfileNodeHandlerRegistry(std::string outputFilename, 
                                MPI_Comm _aggComm,
                                int _aggCommRank)
      : oFilename(outputFilename),
        aggComm(_aggComm),
        aggCommRank(_aggCommRank) {
        // Nothing else to do.
    }


	/**
	 * Destroy our Registry.
	 */
	virtual ~ProfileNodeHandlerRegistry(void) {
        // Nothing else to do.
    }

	/**
     * Aggregate memory usage data collected by any process in the program.
     * Note that we currently write profile data to per-process memory
     * usage profile files.  Therefore, we do not aggregate data in this
     * function - we just have each rank output its own profile data.
     * TODO if we eventually implement memory profile aggregation, this 
     * will have to be changed back to something that outputs the 
     * aggregated profile data to the given output stream.
	 */
    virtual std::shared_ptr<IHandlerRegistry::GlobalData> collectData(void) const;


	/**
	 * Report memory usage data to the given stream.
     * Note that because we write per-process memory profile files, and
     * this function is called only in the rank 0 process, it is a no-op for us.
     * TODO if we eventually implement memory profile aggregation, this will
     * have to be changed back to something that outputs the aggregated
     * profile data to the given output stream.
	 *
	 * @param os Stream on which to output the data.
     * @param stats Data to be reported.
	 */
	virtual void reportData(std::ostream& os,
                        std::shared_ptr<IHandlerRegistry::GlobalData> data) const {
        // Nothing to do.  See comment above.
    }
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_PROFILE_NODE_HANDLER_REGISTRY_H
