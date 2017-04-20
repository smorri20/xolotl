#ifndef XMEMUSAGE_SUMMARY_NODE_HANDLER_REGISTRY_H
#define XMEMUSAGE_SUMMARY_NODE_HANDLER_REGISTRY_H

#include <map>
#include <memory>
#include "xolotlMemUsage/common/CommonHandlerRegistry.h"
#include "xolotlMemUsage/summarynode/NodeMemUsageObjStatistics.h"
#include "xolotlMemUsage/summarynode/NodeMemUsageSampler.h"


namespace xolotlMemUsage {

/**
 * Base class for for building memory usage data collection objects that
 * collect data (as opposed to low-overhead stubs).
 */
class SummaryNodeHandlerRegistry : public CommonHandlerRegistry {
public:
    /**
     * Globally aggregated memory usage statistics for all
     * data types we know how to collect.
     */
    struct GlobalMemUsageStats : public IHandlerRegistry::GlobalData {
        NodeMemUsageObjStatsMap<NodeMemUsageObjStatistics<NodeMemUsageStats> > memStats;
    };

private:
	/**
	 * Collect memory usage data from all program processes.
	 * In the process with rank 0, compute statistics for each
	 * memory usage metrics and populate the stats map with those
	 * statistics.
	 *
	 * @param myRank This process' MPI rank.
	 * @param myObjs My map of memory usage data objects of type T.
	 * @param stats A map of statistics structs to be populated with
	 * statistics for all memory usage data objects of type T across all
	 * processes in the program.  Only meaningful in process with MPI rank 0.
	 */
	template<typename T, typename V>
	void AggregateStatistics(int myRank,
			const std::map<std::string, std::shared_ptr<T> >& myObjs,
			std::map<std::string, NodeMemUsageObjStatistics<V> >& stats) const;

	/**
	 * Retrieve the value of the named object from the given collection of
	 * objects.
	 *
	 * @param myObjs A collection of memory usage data collection objects.
	 * @param objName The name of one of the objects in the collection.
	 * @return The value of the named object.
	 */
	template<typename T, typename V>
    std::shared_ptr<V> GetObjValue(
			const std::map<std::string, std::shared_ptr<T> >& myObjs,
			const std::string& objName) const;

	/**
	 * Build a collection of names of the memory usage objects of type T
	 * from within my own process.  Type T is a type we know about.
	 *
	 * @param myObjs A map of memory usage data collection objects.
	 * @param objNames A vector of memory usage metric names defined by
	 * the objects in myObjs.
	 */
	template<typename T>
	void CollectMyObjectNames(
			const std::map<std::string, std::shared_ptr<T> >& myObjs,
			std::vector<std::string>& objNames) const;

	/**
	 * Build a collection of all known memory usage data metric names across
	 * all processes in the program.
	 *
	 * @param myRank My MPI rank.
	 * @param myObjs A collection of the memory usage data collection objects from my own process.
	 * @param stats A map of partially constructed MemUsageObjStatistics, keyed
	 * by memory usage metric name.
	 * There is one MemUsageObjStatistics item in the map for each memory usage 
	 * metric known across all processes of the program.
	 * This map will only be populated within the process with MPI rank 0.
	 */
	template<typename T, typename V>
	void CollectAllObjectNames(int myRank,
			const std::map<std::string, std::shared_ptr<T> >& myObjs,
			std::map<std::string, NodeMemUsageObjStatistics<V> >& stats) const;


    /**
     * Construct a MemUsageSampler (of the appropriate derived type)
     * with the given name.
     *
     * @param name Name to associate with the sampler.
     */
    virtual std::shared_ptr<IMemUsageSampler> MakeMemUsageSampler(std::string name) {
        return std::make_shared<NodeMemUsageSampler>(name);
    }


    /// Communicator to use for aggregating collected data.
    MPI_Comm aggComm;


    /// Our rank within the aggregator communicator.
    int aggCommRank;


public:

	/**
	 * Construct a SummaryNodeHandlerRegistry.
	 */
    SummaryNodeHandlerRegistry(MPI_Comm _aggComm, int _aggCommRank)
      : aggComm(_aggComm),
        aggCommRank(_aggCommRank) {
        // Nothing else to do.
    }

	/**
	 * Destroy a SummaryNodeHandlerRegistry.
	 */
	virtual ~SummaryNodeHandlerRegistry(void) {
        // Nothing else to do.
    }


	/**
	 * Collect data about any memory usage data collected by
	 * processes of the program.
	 */
    virtual std::shared_ptr<IHandlerRegistry::GlobalData> collectData(void) const;

	/**
	 * Report memory usage data to the given stream.
	 *
	 * @param os Stream on which to output data.
     * @param stats The data to be reported.
	 */
	virtual void reportData(std::ostream& os,
                        std::shared_ptr<IHandlerRegistry::GlobalData> data) const;
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_SUMMARY_NODE_HANDLER_REGISTRY_H
