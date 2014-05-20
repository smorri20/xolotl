#ifndef HDF5NETWORKLOADER_H_
#define HDF5NETWORKLOADER_H_

//Includes
#include <PSIClusterNetworkLoader.h>

namespace xolotlCore {

/**
 * This class overwrites load() for HDF5 files
 */
class HDF5NetworkLoader: public PSIClusterNetworkLoader {

	/**
	 * Private nullary constructor.
	 */
	HDF5NetworkLoader() {
	};

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 */
	HDF5NetworkLoader(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
		handlerRegistry = registry;
	}

	/**
	 * Destructor
	 */
	virtual ~HDF5NetworkLoader() {
	}

	/**
	 * This operation will load the reaction network from the HDF5 file xolotlStart.h5 in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 * @return The reaction network
	 */
	std::shared_ptr<PSIClusterReactionNetwork> load();
};

} /* namespace xolotlCore */
#endif /* HDF5NETWORKLOADER_H_ */
