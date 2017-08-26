#include <algorithm>
#include "NECluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;


void NECluster::createProduction(std::shared_ptr<ProductionReaction> reaction,
		int a, int b) {

	// Add a cluster pair for given reaction 
	reactingPairs.emplace_back((NECluster *) reaction->first,
			(NECluster *) reaction->second);
	// Setup the connectivity array
	setReactionConnectivity(reaction->first->getId());
	setReactionConnectivity(reaction->second->getId());

	return;
}

void NECluster::createCombination(std::shared_ptr<ProductionReaction> reaction,
		int a, int b) {
	setReactionConnectivity(id);
	// Look for the other cluster
	IReactant * secondCluster;
	if (reaction->first->getId() == id)
		secondCluster = reaction->second;
	else
		secondCluster = reaction->first;

	// Add the combining cluster to list of clusters that combine with us
	combiningReactants.emplace_back((NECluster *) secondCluster);

	// Setup the connectivity array
	setReactionConnectivity(id);
	setReactionConnectivity(secondCluster->getId());

	return;
}

void NECluster::createDissociation(
		std::shared_ptr<DissociationReaction> reaction, int a, int b) {
	// Look for the other cluster
	IReactant * emittedCluster;
	if (reaction->first->getId() == id)
		emittedCluster = reaction->second;
	else
		emittedCluster = reaction->first;

	// Add a pair where it is important that the
	// dissociating cluster is the first one
	dissociatingPairs.emplace_back((NECluster *) reaction->dissociating,
			(NECluster *) emittedCluster);

	// Setup the connectivity array
	setDissociationConnectivity(reaction->dissociating->getId());

	return;
}

void NECluster::createEmission(std::shared_ptr<DissociationReaction> reaction,
		int a, int b) {

	// Add the pair of emitted clusters.
	emissionPairs.emplace_back((NECluster *) reaction->first,
                    			(NECluster *) reaction->second);

	// Setup the connectivity array to itself
	setReactionConnectivity(id);

	return;
}

void NECluster::optimizeReactions() {
	// Loop on the pairs to add reactions to the network
    std::for_each(reactingPairs.begin(), reactingPairs.end(),
        [this](ClusterPair& currPair) {
            // Create the corresponding production reaction
            auto newReaction = std::make_shared<ProductionReaction>(currPair.first,
                    currPair.second);
            // Add it to the network
            newReaction = network.addProductionReaction(newReaction);
            // Link it to the pair
            currPair.reaction = newReaction;
        });

    std::for_each(combiningReactants.begin(), combiningReactants.end(),
        [this](CombiningCluster& cc) {
            // Create the corresponding production reaction
            auto newReaction = std::make_shared<ProductionReaction>(cc.combining,
                    this);
            // Add it to the network
            newReaction = network.addProductionReaction(newReaction);
            // Link it to the pair
            cc.reaction = newReaction;
        });

    std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
        [this](ClusterPair& currPair) {
            // Create the corresponding dissociation reaction
            auto newReaction = std::make_shared<DissociationReaction>(currPair.first,
                    currPair.second, this);
            // Add it to the network
            newReaction = network.addDissociationReaction(newReaction);
            // Link it to the pair
            currPair.reaction = newReaction;
        });

    std::for_each(emissionPairs.begin(), emissionPairs.end(),
        [this](ClusterPair& currPair) {
            // Create the corresponding dissociation reaction
            auto newReaction = std::make_shared<DissociationReaction>(this,
                    currPair.first, currPair.second);
            // Add it to the network
            newReaction = network.addDissociationReaction(newReaction);
            // Link it to the pair
            currPair.reaction = newReaction;
        });

	return;
}

static std::vector<int> getFullConnectivityVector(std::set<int> connectivitySet,
		int size) {
	// Create a vector of zeroes with size equal to the network size
	std::vector<int> connectivity(size);

	// Set the value of the connectivity array to one for each element that is
	// in the set.
	for (auto it = connectivitySet.begin(); it != connectivitySet.end(); ++it) {
		connectivity[(*it) - 1] = 1;
	}

	return connectivity;
}

std::vector<int> NECluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network.getDOF());
}

std::vector<int> NECluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network.getDOF());
}

void NECluster::resetConnectivities() {
	// Shrink the arrays to save some space. (About 10% or so.)
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(xeMomId);
	setDissociationConnectivity(xeMomId);

	// Apply to each reacting pair.
    std::for_each(reactingPairs.begin(), reactingPairs.end(),
        [this](const ClusterPair& currPair) {
            // The cluster is connecting to both clusters in the pair
            setReactionConnectivity(currPair.first->id);
            setReactionConnectivity(currPair.first->xeMomId);
            setReactionConnectivity(currPair.second->id);
            setReactionConnectivity(currPair.second->xeMomId);
        });

	// Apply to each combining cluster.
    std::for_each(combiningReactants.begin(), combiningReactants.end(),
        [this](const CombiningCluster& cc) {
            // The cluster is connecting to the combining cluster
            setReactionConnectivity(cc.combining->id);
            setReactionConnectivity(cc.combining->xeMomId);
        });

	// Apply to each effective dissociating pair
    std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
        [this](const ClusterPair& currPair) {
            // The cluster is connecting to the dissociating cluster which
            // is the first one by definition
            setDissociationConnectivity(currPair.first->id);
            setDissociationConnectivity(currPair.first->xeMomId);
        });

	// Don't apply to the emission pairs because
	// this cluster is not connected to them

	return;
}

void NECluster::updateFromNetwork() {

	// Clear the flux-related arrays
	reactingPairs.clear();
	combiningReactants.clear();
	dissociatingPairs.clear();
	emissionPairs.clear();

	return;
}

double NECluster::getMomentum() const {
	return 0.0;
}

double NECluster::getTotalFlux() {
	// Get the fluxes
	double prodFlux = getProductionFlux();
	double dissFlux = getDissociationFlux();
	double combFlux = getCombinationFlux();
	double emissFlux = getEmissionFlux();

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double NECluster::getDissociationFlux() const {

    // Sum dissociation flux over all pairs that dissociate to form this one.
    double flux = std::accumulate(dissociatingPairs.begin(), dissociatingPairs.end(),
            0.0,
            [](double running, const ClusterPair& currPair) {
                // Get the dissociating cluster
                auto& dissociatingCluster = currPair.first;
                // Calculate the Dissociation flux
                return running + (currPair.reaction->kConstant * 
                        dissociatingCluster->getConcentration(currPair.firstDistance));
            });

	// Return the flux
	return flux;
}

double NECluster::getEmissionFlux() const {

    // Sum reaction rate constants over all emission pair reactions.
    double flux = std::accumulate(emissionPairs.begin(), emissionPairs.end(),
            0.0,
            [](double running, const ClusterPair& currPair) {
                return running + currPair.reaction->kConstant;
            });

	return flux * concentration;
}

double NECluster::getProductionFlux() const {
	// Local declarations
	double flux = 0.0;
	int nPairs = 0;
	NECluster *firstReactant = nullptr, *secondReactant = nullptr;

	// Set the total number of reacting pairs
	nPairs = reactingPairs.size();
	// Loop over all the reacting pairs
	for (int i = 0; i < nPairs; i++) {
		// Get the two reacting clusters
		firstReactant = reactingPairs[i].first;
		secondReactant = reactingPairs[i].second;
		// Update the flux
		flux += reactingPairs[i].reaction->kConstant
				* firstReactant->getConcentration(
						reactingPairs[i].firstDistance)
				* secondReactant->getConcentration(
						reactingPairs[i].secondDistance);
	}

	// Return the production flux
	return flux;
}

double NECluster::getCombinationFlux() const {

    double flux = std::accumulate(combiningReactants.begin(), combiningReactants.end(),
            0.0,
            [](double running, const CombiningCluster& currPair) {
                // Get the cluster that combines with this one
                auto& combiningCluster = currPair.combining;
                // Calculate Second term of production flux
                return running + 
                    (currPair.reaction->kConstant * 
                        combiningCluster->getConcentration(currPair.distance));

            });

	return flux * concentration;
}

std::vector<double> NECluster::getPartialDerivatives() const {
	// Local Declarations
	std::vector<double> partials(network.getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return partials;
}

void NECluster::getPartialDerivatives(std::vector<double> & partials) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials);
	getCombinationPartialDerivatives(partials);
	getDissociationPartialDerivatives(partials);
	getEmissionPartialDerivatives(partials);

	return;
}

void NECluster::getProductionPartialDerivatives(
		std::vector<double> & partials) const {
	// Initial declarations
	int numReactants = 0, index = 0;
	double value = 0.0;

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	numReactants = reactingPairs.size();
	for (int i = 0; i < numReactants; i++) {
		// Compute the contribution from the first part of the reacting pair
		value = reactingPairs[i].reaction->kConstant
				* reactingPairs[i].second->getConcentration(
						reactingPairs[i].secondDistance);
		index = reactingPairs[i].first->id - 1;
		partials[index] += value;
		index = reactingPairs[i].first->xeMomId - 1;
		partials[index] += value * reactingPairs[i].firstDistance;
		// Compute the contribution from the second part of the reacting pair
		value = reactingPairs[i].reaction->kConstant
				* reactingPairs[i].first->getConcentration(
						reactingPairs[i].firstDistance);
		index = reactingPairs[i].second->id - 1;
		partials[index] += value;
		index = reactingPairs[i].second->xeMomId - 1;
		partials[index] += value * reactingPairs[i].secondDistance;
	}

	return;
}

void NECluster::getCombinationPartialDerivatives(
		std::vector<double> & partials) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
    std::for_each(combiningReactants.begin(), combiningReactants.end(),
        [this,&partials](const CombiningCluster& cc) {

            NECluster* cluster = (NECluster *) cc.combining;
            // Remember that the flux due to combinations is OUTGOING (-=)!
            // Compute the contribution from this cluster
            partials[id - 1] -= cc.reaction->kConstant * 
                cluster->getConcentration(cc.distance);
            // Compute the contribution from the combining cluster
            double value = cc.reaction->kConstant * concentration;

            partials[cluster->id - 1] -= value;
            partials[cluster->xeMomId - 1] -= value * cc.distance;
        });

	return;
}

void NECluster::getDissociationPartialDerivatives(
		std::vector<double> & partials) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
    std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
        [&partials](const ClusterPair& currPair) {
            // Get the dissociating cluster
            NECluster* cluster = currPair.first;
            partials[cluster->id - 1] += currPair.reaction->kConstant;
            partials[cluster->xeMomId - 1] += currPair.reaction->kConstant * 
                currPair.firstDistance;
        });

	return;
}

void NECluster::getEmissionPartialDerivatives(
		std::vector<double> & partials) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)
    double emissionFlux = std::accumulate(emissionPairs.begin(), emissionPairs.end(),
            0.0,
            [](double running, const ClusterPair& currPair) {
                return running + currPair.reaction->kConstant;
            });

    // Recall emission flux is OUTGOING
    partials[id - 1] -= emissionFlux;

	return;
}

void NECluster::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

void NECluster::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

double NECluster::getLeftSideRate() const {

    // Sum reaction rate contributions over all combining clusters.
    double combiningRateTotal = std::accumulate(combiningReactants.begin(), combiningReactants.end(),
            0.0,
            [](double running, const CombiningCluster& currPair) {
                auto const& cluster = (NECluster*) currPair.combining;

                return running + (currPair.reaction->kConstant * 
                                    cluster->concentration);
            });

	// Sum reaction rate constants over all emission pairs.
    double emissionRateTotal = std::accumulate(emissionPairs.begin(), emissionPairs.end(),
            0.0,
            [](double running, const ClusterPair& currPair) {
                return running + currPair.reaction->kConstant;
            });

	return combiningRateTotal + emissionRateTotal;
}

std::vector<int> NECluster::getConnectivity() const {
	int connectivityLength = network.getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The reaction vector has an incorrect length");
	}
	if (dissociationConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The dissociation vector has an incorrect length");
	}

	// Merge the two vectors such that the final vector contains
	// a 1 at a position if either of the connectivity arrays
	// have a 1
	for (int i = 0; i < connectivityLength; i++) {
		// Consider each connectivity array only if its type is enabled
		connectivity[i] = reactionConnVector[i] || dissociationConnVector[i];
	}

	return connectivity;
}
