#ifndef HDF5NETWORKLOADER_H_
#define HDF5NETWORKLOADER_H_

//Includes
#include <PSIClusterNetworkLoader.h>

namespace xolotlCore {

/**
 * This class overwrites the load() methods of PSIClusterNetworkLoader
 * for HDF5 files.
 */
class HDF5NetworkLoader: public PSIClusterNetworkLoader {
private:

	/**
	 * Name of the file to load the network from.
	 */
	std::string fileName;

	/**
	 * The helium size at wich the grouping scheme starts
	 */
	int heMin;

	/**
	 * The number of clusters that will be gathered in the grouping scheme
	 */
	int sectionWidth;

	/**
	 * Private nullary constructor.
	 */
	HDF5NetworkLoader() {}

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 */
	HDF5NetworkLoader(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry):
		fileName(""),
		heMin(std::numeric_limits<int>::max()),
		sectionWidth(1) {
		handlerRegistry = registry;

		return;
	}

	/**
	 * The destructor.
	 */
	virtual ~HDF5NetworkLoader() {
	}

	/**
	 * This operation will load the reaction network from the HDF5 file in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 *
	 * @return The reaction network.
	 */
	std::shared_ptr<PSIClusterReactionNetwork> load();

	/**
	 * This operation will apply a sectional grouping method to the network.
	 *
	 * @param The network to be modified.
	 */
	void applySectionalGrouping(std::shared_ptr<PSIClusterReactionNetwork> network);

	/**
	 * Find the super cluster that will replace the HeV cluster given by its composition.
	 *
	 * @param clusterList The list of super clusters.
	 * @param comp The composition of the HeV cluster.
	 * @return The pointer to the super cluster
	 */
	PSICluster * findSuperCluster(std::vector<Reactant *> clusterList,
			std::map<std::string, int> comp);

	/**
	 * This operation will set the name of the file where to take the network from.
	 *
	 * @param name The name of the file
	 */
	void setFilename (const std::string& name) {fileName = name;}

	/**
	 * This operation will get the name of the file where to take the network from.
	 *
	 * @return The name of the file
	 */
	std::string getFilename () const {return fileName;}

	/**
	 * This operation will set the helium size at which the grouping scheme starts.
	 *
	 * @param min The value for the size
	 */
	void setHeMin (int min) {heMin = min;}

	/**
	 * This operation will set the width for the grouping scheme.
	 *
	 * @param w The value of the width
	 */
	void setWidth (int w) {sectionWidth = w;}
};

} /* namespace xolotlCore */

#endif /* HDF5NETWORKLOADER_H_ */
