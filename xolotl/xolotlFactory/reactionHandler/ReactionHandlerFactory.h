#ifndef REACTIONHANDLERFACTORY_H
#define REACTIONHANDLERFACTORY_H

#include <memory>
#include <IReactionNetwork.h>
#include <HDF5NetworkLoader.h>
#include <Options.h>

namespace xolotlFactory {

/**
 * Build the desired type of reaction network.
 *
 * @param networkLoader The network loader
 * @param Options The options
 * @return True if the network was created successfully.
 */
bool initializeReactionHandler(std::shared_ptr<xolotlCore::HDF5NetworkLoader> networkLoader,
		xolotlCore::Options &options);

/**
 * Access the network.
 *
 *  @return The network object.
 */
std::shared_ptr<xolotlCore::IReactionNetwork> getNetworkHandler();

} // end namespace xolotlFactory

#endif /* REACTIONHANDLERFACTORY_H */
