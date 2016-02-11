#ifndef HDF5NETWORKLOADER_H_
#define HDF5NETWORKLOADER_H_

//Includes
#include <PSIClusterNetworkLoader.h>
#include "SuperCluster.h"

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
	 * The vacancy size at which the grouping scheme starts
	 */
	int vMin;

	/**
	 * The width of the group in the helium direction.
	 */
	int heSectionWidth;

	/**
	 * The width of the group in the vacancy direction.
	 */
	int vSectionWidth;

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
		vMin(std::numeric_limits<int>::max()),
		heSectionWidth(1),
		vSectionWidth(1) {
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
	void setVMin (int min) {vMin = min;}

	/**
	 * This operation will set the helium width for the grouping scheme.
	 *
	 * @param w The value of the width
	 */
	void setHeWidth (int w) {heSectionWidth = w;}

	/**
	 * This operation will set the vacancy width for the grouping scheme.
	 *
	 * @param w The value of the width
	 */
	void setVWidth (int w) {vSectionWidth = w;}
};

} /* namespace xolotlCore */

#endif /* HDF5NETWORKLOADER_H_ */
