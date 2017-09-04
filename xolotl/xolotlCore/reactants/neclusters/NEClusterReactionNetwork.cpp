#include "NEClusterReactionNetwork.h"
#include "NECluster.h"
#include "NESuperCluster.h"
#include <xolotlPerf.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Constants.h>

using namespace xolotlCore;


NEClusterReactionNetwork::NEClusterReactionNetwork(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		ReactionNetwork( {
                    ReactantType::V,
                    ReactantType::I,
                    ReactantType::Xe,
                    ReactantType::XeV,
                    ReactantType::XeI,
                    ReactantType::NESuper
                },
                ReactantType::NESuper, registry) {

	// Initialize default properties
	dissociationsEnabled = true;

	return;
}

double NEClusterReactionNetwork::calculateDissociationConstant(const DissociationReaction& reaction) const {

	// If the dissociations are not allowed
	if (!dissociationsEnabled)
		return 0.0;

	// Compute the atomic volume
	double atomicVolume = 0.5 * xolotlCore::uraniumDioxydeLatticeConstant
			* xolotlCore::uraniumDioxydeLatticeConstant
			* xolotlCore::uraniumDioxydeLatticeConstant;

	// Get the rate constant from the reverse reaction
	double kPlus = reaction.reverseReaction->kConstant;

	// Calculate and return
	double bindingEnergy = computeBindingEnergy(reaction);
	double k_minus_exp = exp(
			-1.0 * bindingEnergy / (xolotlCore::kBoltzmann * temperature));
	double k_minus = (1.0 / atomicVolume) * kPlus * k_minus_exp;

	return k_minus;
}

void NEClusterReactionNetwork::createReactionConnectivity() {
	// Initial declarations
	int firstSize = 0, secondSize = 0, productSize = 0;

	// Single species clustering (Xe)
	// We know here that only Xe_1 can cluster so we simplify the search
	// Xe_(a-i) + Xe_i --> Xe_a
	firstSize = 1;
	auto singleXeCluster = get(ReactantType::Xe, firstSize);
	// Get all the Xe clusters
	auto const& xeClusters = getAll(ReactantType::Xe);
	// Consider each Xe cluster.
    for (auto const& currMapItem : getAll(ReactantType::Xe)) {

        auto& xeReactant = *(currMapItem.second);

		// Get the size of the second reactant and product
		secondSize = xeReactant.getSize();
		productSize = firstSize + secondSize;
		// Get the product cluster for the reaction
		auto product = get(ReactantType::Xe, productSize);
		// Check that the reaction can occur
		if (product
				&& (singleXeCluster->getDiffusionFactor() > 0.0
						|| xeReactant.getDiffusionFactor() > 0.0)) {
			// Create a production reaction
			auto reaction = std::make_shared<ProductionReaction>(
					singleXeCluster, &xeReactant);
			// Tell the reactants that they are in this reaction
			singleXeCluster->createCombination(reaction);
			xeReactant.createCombination(reaction);
			product->createProduction(reaction);

			// Check if the reverse reaction is allowed
			checkDissociationConnectivity(product, reaction);
		}
	}

	return;
}

void NEClusterReactionNetwork::checkDissociationConnectivity(
		IReactant * emittingReactant,
		std::shared_ptr<ProductionReaction> reaction) {
	// Check if at least one of the potentially emitted cluster is size one
	if (reaction->first->getSize() != 1 && reaction->second->getSize() != 1) {
		// Don't add the reverse reaction
		return;
	}

	// The reaction can occur, create the dissociation
	// Create a dissociation reaction
	auto dissociationReaction = std::make_shared<DissociationReaction>(
			emittingReactant, reaction->first, reaction->second);
	// Set the reverse reaction
	dissociationReaction->reverseReaction = reaction.get();
	// Tell the reactants that their are in this reaction
	reaction->first->createDissociation(dissociationReaction);
	reaction->second->createDissociation(dissociationReaction);
	emittingReactant->createEmission(dissociationReaction);

	return;
}

void NEClusterReactionNetwork::setTemperature(double temp) {
	ReactionNetwork::setTemperature(temp);

	computeRateConstants();

	return;
}

double NEClusterReactionNetwork::getTemperature() const {
	return temperature;
}

IReactant * NEClusterReactionNetwork::get(ReactantType type,
		IReactant::SizeType size) const {
	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	// Only pull the reactant if the name and size are valid
	if ((type == ReactantType::Xe || type == ReactantType::V || type == ReactantType::I) && size >= 1) {
        IReactant::Composition composition;
		composition[toCompIdx(toSpecies(type))] = size;

		// Check if we know about the reactant.
        auto const& currReactantMap = clusterTypeMap.at(type);
        auto iter = currReactantMap.find(composition);
        if (iter != currReactantMap.end()) {
            retReactant = iter->second;
        }
	}

	return retReactant.get();
}

IReactant * NEClusterReactionNetwork::getCompound(ReactantType type,
        const IReactant::Composition& comp) const {

	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	// Only pull the reactant if the name is valid and there are enough sizes
	// to fill the composition.
	if (type == ReactantType::XeV || type == ReactantType::XeI) {

		// Check if the reactant is in the map
        auto const& currReactantMap = clusterTypeMap.at(type);
        auto iter = currReactantMap.find(comp);
        if (iter != currReactantMap.end()) {
            retReactant = iter->second;
        }
	}

	return retReactant.get();
}

IReactant * NEClusterReactionNetwork::getSuper(ReactantType type,
		IReactant::SizeType size) const {
	// Local Declarations
	std::shared_ptr<IReactant> retReactant;

	// Only pull the reactant if the name and size are valid.
	if (type == ReactantType::NESuper && size >= 1) {
        IReactant::Composition composition;
		composition[toCompIdx(Species::Xe)] = size;

		// Check if the reactant is in the map
        auto const& currTypeMap = clusterTypeMap.at(type);
        auto iter = currTypeMap.find(composition);
        if (iter != currTypeMap.end()) {
            retReactant = iter->second;
		}
	}

	return retReactant.get();
}


void NEClusterReactionNetwork::add(std::shared_ptr<IReactant> reactant) {

	// Only add a complete reactant
	if (reactant) {
		// Get the composition
		auto& composition = reactant->getComposition();

        // Get the species sizes
        auto numXe = composition[toCompIdx(Species::Xe)] ;
        auto numV = composition[toCompIdx(Species::V)] ;
        auto numI = composition[toCompIdx(Species::I)] ;

        // Add element if we don't already know about it.
        auto& currTypeMap = clusterTypeMap[reactant->getType()];
        auto iter = currTypeMap.find(composition);
        if (iter == currTypeMap.end()) {

            // Set the id for this cluster
            // (It is networkSize+1 because we haven't added 
            // it to the network yet.)
            reactant->setId(size() + 1);

            // Add to our per-type map.
            currTypeMap.emplace(reactant->getComposition(), reactant);

            // Add the pointer to the list of all clusters
            allReactants.push_back(*reactant);

            // Update max cluster size for this type if necessary.
            // TODO if clusterSize always is same as reactant->getSize,
            // just use reactant->getSize and avoid extracting
            // individual species sizes except in the exception case.
            // (See PSI implementation for an example.)
            IReactant::SizeType clusterSize = numXe + numV + numI;
            assert(clusterSize == reactant->getSize());
            maxClusterSizeMap[reactant->getType()] = std::max(reactant->getSize(), maxClusterSizeMap[reactant->getType()]);

        } else {

            std::stringstream errStream;
            errStream << "NEClusterReactionNetwork Message: "
                    << "Duplicate Reactant (Xe=" << numXe << ",V=" << numV
                    << ",I=" << numI << ") not added!" << std::endl;
            throw errStream.str();
        }
	}

	return;
}

void NEClusterReactionNetwork::addSuper(std::shared_ptr<IReactant> reactant) {

	// Only add a complete reactant
	if (reactant) {
		// Get the composition
		auto& composition = reactant->getComposition();

        // Add the element if we don't already know about it.
        auto& currTypeMap = clusterTypeMap[reactant->getType()];
        auto iter = currTypeMap.find(composition);
        if (iter == currTypeMap.end()) {

            // Set the id for this cluster
            // (It is networkSize+1 because we haven't added it 
            // to the network yet.)
            reactant->setId(size() + 1);

            // Add the pointer to the list of all clusters
            allReactants.push_back(*reactant);

            // Add to our per-type map.
            currTypeMap.emplace(reactant->getComposition(), reactant);

		} else {
            // Get the species sizes
            auto numXe = composition[toCompIdx(Species::Xe)] ;
            auto numV = composition[toCompIdx(Species::V)] ;
            auto numI = composition[toCompIdx(Species::I)] ;

			std::stringstream errStream;
			errStream << "NEClusterReactionNetwork Message: "
					<< "Duplicate Super Reactant (Xe=" << numXe << ",V=" << numV
					<< ",I=" << numI << ") not added!" << std::endl;
			throw errStream.str();
		}
	}

	return;
}

// TODO check if this is identical to the PSI vesrion,
// and move upwards in class hierarchy if so.
void NEClusterReactionNetwork::removeReactants(
		const IReactionNetwork::ReactantVector& doomedReactants) {

	// Build a ReactantMatcher functor for the doomed reactants.
	// Doing this here allows us to construct the canonical composition
	// strings for the doomed reactants once and reuse them.
	// If we used an anonymous functor object in the std::remove_if
	// calls we would build these strings several times in this function.
	ReactionNetwork::ReactantMatcher doomedReactantMatcher(doomedReactants);

	// Remove the doomed reactants from our collection of all known reactants.
	auto ariter = std::remove_if(allReactants.begin(), allReactants.end(),
			doomedReactantMatcher);
	allReactants.erase(ariter, allReactants.end());

	// Remove the doomed reactants from the type-specific cluster vectors.
	// First, determine all cluster types used by clusters in the collection
	// of doomed reactants...
	std::set<ReactantType> typesUsed;
	for (auto reactant : doomedReactants) {
		typesUsed.insert(reactant->getType());
	}

	// ...Next, examine each type's collection of clusters and remove the
	// doomed reactants.
	for (auto currType : typesUsed) {
		auto& clusters = clusterTypeMap[currType];
        for (auto const& currDoomedReactant : doomedReactants) {
            auto iter = clusters.find(currDoomedReactant->getComposition());
            assert(iter != clusters.end());
            clusters.erase(iter);
        }
	}

	return;
}

void NEClusterReactionNetwork::reinitializeNetwork() {

	// Reset the Ids
    // TODO std::for_each visits reactants in-order in C++11, so ids are
    // assigned sequentially.
	int id = 0;
    std::for_each(allReactants.begin(), allReactants.end(),
        [&id](IReactant& currReactant) {
            id++;
            currReactant.setId(id);
            currReactant.setXeMomentumId(id);

            currReactant.optimizeReactions();
        });

	// Get all the super clusters and loop on them
    for (auto const& currMapItem : clusterTypeMap.at(ReactantType::NESuper)) {

        auto& currCluster = static_cast<NESuperCluster&>(*(currMapItem.second));
		id++;
		currCluster.setXeMomentumId(id);
	}

	return;
}

void NEClusterReactionNetwork::reinitializeConnectivities() {
	// Loop on all the reactants to reset their connectivities
    std::for_each(allReactants.begin(), allReactants.end(),
        [](IReactant& currReactant) {
            currReactant.resetConnectivities();
        });

	return;
}

void NEClusterReactionNetwork::updateConcentrationsFromArray(
		double * concentrations) {

	// Set the concentrations of each reactant.
    std::for_each(allReactants.begin(), allReactants.end(),
        [&concentrations](IReactant& currReactant) {
            auto id = currReactant.getId() - 1;
            currReactant.setConcentration(concentrations[id]);
        });

	// Set the moments
    for (auto const& currMapItem : getAll(ReactantType::NESuper)) {
		// Get the superCluster
        auto& cluster = static_cast<NESuperCluster&>(*(currMapItem.second));

		auto id = cluster.getId() - 1;
		cluster.setZerothMomentum(concentrations[id]);
		id = cluster.getXeMomentumId() - 1;
		cluster.setMomentum(concentrations[id]);
	}

	return;
}

void NEClusterReactionNetwork::getDiagonalFill(int *diagFill) {

	// Degrees of freedom is the total number of clusters in the network
	const int dof = getDOF();

	// Get the connectivity for each reactant
    std::for_each(allReactants.begin(), allReactants.end(),
        [&diagFill, &dof, this](const IReactant& reactant) {
            // Get the reactant and its connectivity
            auto const& connectivity = reactant.getConnectivity();
            auto connectivityLength = connectivity.size();
            // Get the reactant id so that the connectivity can be lined up in
            // the proper column
            auto id = reactant.getId() - 1;
            // Create the vector that will be inserted into the dFill map
            std::vector<int> columnIds;
            // Add it to the diagonal fill block
            for (int j = 0; j < connectivityLength; j++) {
                // The id starts at j*connectivity length and is always offset
                // by the id, which denotes the exact column.
                auto index = id * dof + j;
                diagFill[index] = connectivity[j];
                // Add a column id if the connectivity is equal to 1.
                if (connectivity[j] == 1) {
                    columnIds.push_back(j);
                }
            }
            // Update the map
            dFillMap[id] = columnIds;
        });

	// Get the connectivity for each moment
    for (auto const& currMapItem : getAll(ReactantType::NESuper)) {

		// Get the reactant and its connectivity
		auto const& reactant = static_cast<NESuperCluster&>(*(currMapItem.second));

		auto const& connectivity = reactant.getConnectivity();
		auto connectivityLength = connectivity.size();
		// Get the xenon momentum id so that the connectivity can be lined up in
		// the proper column
		auto id = reactant.getXeMomentumId() - 1;

		// Create the vector that will be inserted into the dFill map
		std::vector<int> columnIds;
		// Add it to the diagonal fill block
		for (int j = 0; j < connectivityLength; j++) {
			// The id starts at j*connectivity length and is always offset
			// by the id, which denotes the exact column.
			auto index = (id) * dof + j;
			diagFill[index] = connectivity[j];
			// Add a column id if the connectivity is equal to 1.
			if (connectivity[j] == 1) {
				columnIds.push_back(j);
			}
		}
		// Update the map
		dFillMap[id] = columnIds;
	}

	return;
}

void NEClusterReactionNetwork::computeRateConstants() {
	// Local declarations
	double rate = 0.0;
	// Initialize the value for the biggest production rate
	double biggestProductionRate = 0.0;

	// Loop on all the production reactions
    for (auto& currReactionInfo : productionReactionMap) {
			
        auto& currReaction = currReactionInfo.second;

		// Compute the rate
		rate = calculateReactionRateConstant(*currReaction);
		// Set it in the reaction
		currReaction->kConstant = rate;

		// Check if the rate is the biggest one up to now
		if (rate > biggestProductionRate)
			biggestProductionRate = rate;
	}

	// Loop on all the dissociation reactions
    for (auto& currReactionInfo : dissociationReactionMap) {

        auto& currReaction = currReactionInfo.second;

		// Compute the rate
		rate = calculateDissociationConstant(*currReaction);
		// Set it in the reaction
		currReaction->kConstant = rate;
	}

	// Set the biggest rate
	biggestRate = biggestProductionRate;

	return;
}

void NEClusterReactionNetwork::computeAllFluxes(double *updatedConcOffset) {

	// ----- Compute all of the new fluxes -----
    std::for_each(allReactants.begin(), allReactants.end(),
        [&updatedConcOffset](IReactant& cluster) {

            // Compute the flux
            auto flux = cluster.getTotalFlux();
            // Update the concentration of the cluster
            auto reactantIndex = cluster.getId() - 1;
            updatedConcOffset[reactantIndex] += flux;
        });

	// ---- Moments ----
    for (auto const& currMapItem : getAll(ReactantType::NESuper)) {

        auto& superCluster = static_cast<NESuperCluster&>(*(currMapItem.second));

		// Compute the xenon momentum flux
		auto flux = superCluster.getMomentumFlux();
		// Update the concentration of the cluster
		auto reactantIndex = superCluster.getXeMomentumId() - 1;
		updatedConcOffset[reactantIndex] += flux;
	}

	return;
}

void NEClusterReactionNetwork::computeAllPartials(double *vals, int *indices,
		int *size) {

	// Initial declarations
	const int dof = getDOF();
	std::vector<double> clusterPartials;
	clusterPartials.resize(dof, 0.0);

	// Get the super clusters
	auto const& superClusters = getAll(ReactantType::NESuper);

	// Update the column in the Jacobian that represents each normal reactant
    for (auto const& superMapItem : superClusters) {
        auto const& reactant = *(superMapItem.second);

		// Get the reactant index
		auto reactantIndex = reactant.getId() - 1;

		// Get the partial derivatives
		reactant.getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		auto pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];

			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	// Update the column in the Jacobian that represents the moment for the super clusters
    for (auto const& currMapItem : superClusters) {
		auto const& reactant = static_cast<NESuperCluster&>(*(currMapItem.second));

		// Get the super cluster index
		auto reactantIndex = reactant.getId() - 1;

		// Get the partial derivatives
		reactant.getPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		auto pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		auto pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}

		// Get the helium momentum index
		reactantIndex = reactant.getXeMomentumId() - 1;

		// Get the partial derivatives
		reactant.getMomentPartialDerivatives(clusterPartials);
		// Get the list of column ids from the map
		pdColIdsVector = dFillMap.at(reactantIndex);
		// Number of partial derivatives
		pdColIdsVectorSize = pdColIdsVector.size();
		size[reactantIndex] = pdColIdsVectorSize;

		// Loop over the list of column ids
		for (int j = 0; j < pdColIdsVectorSize; j++) {
			// Set the index
			indices[reactantIndex * dof + j] = pdColIdsVector[j];
			// Get the partial derivative from the array of all of the partials
			vals[reactantIndex * dof + j] = clusterPartials[pdColIdsVector[j]];

			// Reset the cluster partial value to zero. This is much faster
			// than using memset.
			clusterPartials[pdColIdsVector[j]] = 0.0;
		}
	}

	return;
}

