#ifndef IREACTIONHANDLERFACTORY_H
#define IREACTIONHANDLERFACTORY_H

#include <Options.h>
#include <INetworkLoader.h>
#include <IReactionNetwork.h>

namespace xolotlFactory {

/**
 * Realizations of this interface are responsible for handling the flux and the advection.
 * they are both dependent on the type of material under study.
 */
class IReactionHandlerFactory {
public:

	/**
	 * The destructor
	 */
	~IReactionHandlerFactory() {
	}

	/**
	 * Initialize the reaction network.
	 *
	 * @param options The options.
	 * @param registry The performance registry.
	 * @param bounds1 The bounds for the grouping.
	 * @param bounds2 The bounds for the grouping.
	 * @param padeVector The vector containing the Pade approximation for each cluster
	 * @param idMap The idea map from the previously built network
	 */
	virtual void initializeReactionNetwork(xolotlCore::Options &options,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry,
			const std::vector<xolotlCore::IReactant::SizeType> & bounds1,
			const std::vector<xolotlCore::IReactant::SizeType> & bounds2,
			std::vector<std::vector<double> > & padeVector,
			std::map<std::string, int> & idMap) = 0;

	/**
	 * Return the network loader.
	 *
	 * @return The network loader.
	 */
	virtual std::shared_ptr<xolotlCore::INetworkLoader> getNetworkLoaderHandler() const = 0;

	/**
	 * Return the network.
	 *
	 * @return The network.
	 */
	virtual xolotlCore::IReactionNetwork& getNetworkHandler() const = 0;

	/**
	 * Function that create the wanted reaction handler factory depending on the given type.
	 *
	 * @param problemType The type of wanted problem (PSI or NE).
	 * @return The reaction factory.
	 */
	static std::shared_ptr<IReactionHandlerFactory> createNetworkFactory(
			const std::string& problemType);

};

} // end namespace xolotlFactory

#endif // IREACTIONHANDLERFACTORY_H
