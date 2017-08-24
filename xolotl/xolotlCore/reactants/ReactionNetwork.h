#ifndef REACTION_NETWORK_H
#define REACTION_NETWORK_H

// Includes
#include "IReactionNetwork.h"
#include <Constants.h>
#include <set>
#include <map>
#include <unordered_map>

namespace xolotlPerf {
class IHandlerRegistry;
class IEventCounter;
}

namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants
 *  (combinations of normal reactants). It also manages a set of properties
 *  that describe both.
 */
class ReactionNetwork: public IReactionNetwork {

protected:

    /**
     * Nice name for map of reactants, keyed by their composition.
     */
    using ReactantMap = std::map<IReactant::Composition, std::shared_ptr<IReactant> >;


	/**
	 * A functor useful for identifying a set of reactants by their
	 * composition from a container, e.g., when removing a collection
	 * of doomed reactants from a vector of reactants.
	 */
	class ReactantMatcher {
	private:
		/**
		 * The canonical composition string representations of the reactants
		 * we are to find.
		 */
		std::set<IReactant::Composition> comps;

	public:
		/**
		 * Build a ReactantMatcher.
		 * @param reactants The collection of reactants we are to recognize.
		 */
		ReactantMatcher(const IReactionNetwork::ReactantVector& reactants) {
			for (auto reactant : reactants) {
				comps.insert(reactant->getComposition());
			}
		}

		/**
		 * Determine if the given reactant is in our set.
		 * @param testReactant The reactant to check.
		 * @return true iff the reactant's composition's canonical string
		 * representation is in our set.
		 */
		bool operator()(const IReactant* testReactant) const {
			auto iter = comps.find(testReactant->getComposition());
			return (iter != comps.end());
		}

		/**
		 * Determine if the given reactant is in our set.
		 * @param testReactant The reactant to check.
		 * @return true iff the reactant's composition's canonical string
		 * representation is in our set.
		 */
		bool operator()(const std::shared_ptr<IReactant> testReactant) const {
			return this->operator()(testReactant.get());
		}
	};

	/**
	 * The performance handler registry that will be used with
	 * this class.
	 */
	std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry;

	/**
	 * Counter for the number of times the network concentration is updated.
	 */
	std::shared_ptr<xolotlPerf::IEventCounter> concUpdateCounter;

    /**
     * All known ProductionReactions in the network, keyed by a
     * representation of the reaction.
     */
	std::map<ProductionReaction::KeyType, std::shared_ptr<ProductionReaction> > productionReactionMap;

	/**
	 * All known dissociation reactions in the network, keyed by a
     * representation of the reaction.
	 */
	std::map<DissociationReaction::KeyType,
			std::shared_ptr<DissociationReaction> > dissociationReactionMap;

	/**
	 * A map for storing the dfill configuration and accelerating the formation of
	 * the Jacobian. Its keys are reactant/cluster ids and its values are integer
	 * vectors of the column ids that are marked as connected for that cluster in
	 * the dfill array.
	 */
	std::unordered_map<int, std::vector<int> > dFillMap;

	/**
	 * The current temperature at which the network's clusters exist.
	 */
	double temperature;

	/**
	 * The size of the network. It is also used to set the id of new reactants
	 * that are added to the network.
	 */
	int networkSize;

#if READY
	/**
	 * The names of the reactants supported by this network.
	 */
    std::vector<Species> names;

	/**
	 * The names of the compound reactants supported by this network.
	 */
    std::vector<Species> compoundNames;
#endif // READY

	/**
	 * The biggest rate for this cluster
	 */
	double biggestRate;

	/**
	 * Are dissociations enabled?
	 */
	bool dissociationsEnabled;

	/**
	 * Number of vacancy clusters in our network.
	 */
	int numVClusters;

	/**
	 * Number of interstitial clusters in our network.
	 */
	int numIClusters;

	/**
	 * Number of super clusters in our network.
	 */
	int numSuperClusters;

	/**
	 * Maximum size of a vacancy cluster.
	 */
    IReactant::SizeType maxVClusterSize;

	/**
	 * Maximum size of an interstitial cluster.
	 */
    IReactant::SizeType maxIClusterSize;

	/**
	 * This vector contains the information on the group bounds in both directions.
	 */
	std::vector<int> boundVector;


	/**
     * All reactants known to the network.
	 */
    std::vector<IReactant*> allReactants;


	/**
	 * This map stores all of the clusters in the network by type.
	 */
	std::map<ReactantType, ReactantVector> clusterTypeMap;


	/**
	 * The map of single-species clusters, indexed by their composition.
	 */
    ReactantMap singleSpeciesMap;

	/**
	 * The map of mixed or compound species clusters, indexed by 
     * their composition.
	 */
    ReactantMap mixedSpeciesMap;

	/**
	 * The map of super species clusters, indexed by their composition.
	 */
    ReactantMap superSpeciesMap;

	/**
	 * Calculate the reaction constant dependent on the
	 * reaction radii and the diffusion coefficients for the
	 * ith and jth clusters, which itself depends on the current
	 * temperature.
	 *
	 * @param reaction The reaction
	 * @return The rate
	 */
	double calculateReactionRateConstant(const ProductionReaction& reaction) const;

	/**
	 * Calculate the dissociation constant of the first cluster with respect to
	 * the single-species cluster of the same type based on the current clusters
	 * atomic volume, reaction rate constant, and binding energies.
	 *
	 * Need to be overwritten by daughter classes because of the atomic volume.
	 *
	 * @param reaction The reaction
	 * @return The dissociation constant
	 */
	virtual double calculateDissociationConstant(const DissociationReaction& reaction) const = 0;

	/**
	 * Calculate the binding energy for the dissociation cluster to emit the single
	 * and second cluster.
	 *
	 * @param reaction The reaction
	 * @return The binding energy corresponding to this dissociation
	 */
	virtual double computeBindingEnergy(const DissociationReaction& reaction) const = 0;


public:

    /**
     * Default constructor, deleted because we need params to construct.
     */
    ReactionNetwork() = delete;


	/**
	 * The constructor that takes the performance handler registry.
	 * It initializes the properties and reactants vector.
	 *
	 * @param registry The performance handler registry
	 */
	ReactionNetwork(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);


	/**
	 * The destructor.
	 */
	virtual ~ReactionNetwork() {
	}

	/**
	 * This operation sets the temperature at which the reactants currently
	 * exists. It calls setTemperature() on each reactant.
	 *
	 * This is the simplest way to set the temperature for all reactants.
	 *
	 * @param temp The new temperature
	 */
	virtual void setTemperature(double temp);

	/**
	 * This operation returns the temperature at which the cluster currently exists.
	 *
	 * @return The temperature
	 */
	virtual double getTemperature() const;


	/**
	 * This operation returns all reactants in the network without regard for
	 * their composition or whether they are compound reactants. The list may
	 * or may not be ordered and the decision is left to implementers.
	 *
	 * @return The list of all of the reactants in the network
	 */
    virtual const std::vector<IReactant*> & getAll() const override {

        return allReactants;
    }

	/**
	 * This operation returns all reactants in the network with the given type.
	 * The list may or may not be ordered and the decision is left to
	 * implementers.
	 *
	 * @param type The reactant or compound reactant type
	 * @return The list of all of the reactants in the network or null if the
	 * type is invalid
	 */
    virtual const IReactionNetwork::ReactantVector& getAll(ReactantType type) const override {

        return clusterTypeMap.at(type);
    }

	/**
	 * This operation adds a reactant or a compound reactant to the network.
	 * Adding a reactant to the network does not set the network as the
	 * reaction network for the reactant. This step must be performed
	 * separately to allow for the scenario where the network is generated
	 * entirely before running.
	 *
	 * @param reactant The reactant that should be added to the network
	 */
	virtual void add(std::shared_ptr<IReactant> reactant) {
		return;
	}

	/**
	 * This operation adds a super reactant to the network.
	 * Used with a grouping method.
	 *
	 * @param reactant The reactant that should be added to the network
	 */
	virtual void addSuper(std::shared_ptr<IReactant> reactant) {
		return;
	}

	/**
	 * This operation reinitializes the network.
	 *
	 * It computes the cluster Ids and network size from the allReactants vector.
	 */
	virtual void reinitializeNetwork() {
		return;
	}

#if READY
	/**
	 * This operation returns the names of the reactants in the network.
	 *
	 * @return A vector with one entry for each of the distinct reactant types
	 * in the network
	 */
	const std::vector<Species> & getNames() const override {

        return names;
    }

	/**
	 * This operation returns the names of the compound reactants in the
	 * network.
	 *
	 * @return A vector with one each for each of the distinct compound
	 * reactant types in the network
	 */
	const std::vector<Species> & getCompoundNames() const override {

        return compoundNames;
    }
#endif // READY

	/**
	 * This operation returns the size or number of reactants in the network.
	 *
	 * @return The number of reactants in the network
	 */
	int size() {
		return networkSize;
	}

	/**
	 * This operation returns the size or number of reactants and momentums in the network.
	 *
	 * @return The number of degrees of freedom
	 */
	virtual int getDOF() {
		return networkSize;
	}

	/**
	 * This operation adds a production reaction to the network.
	 *
	 * @param reaction The reaction that should be added to the network
	 * @return The pointer to the reaction that is now in the network
	 */
	virtual std::shared_ptr<ProductionReaction> addProductionReaction(
			std::shared_ptr<ProductionReaction> reaction);

	/**
	 * This operation adds a dissociation reaction to the network.
	 *
	 * @param reaction The reaction that should be added to the network
	 * @return The pointer to the reaction that is now in the network
	 */
	virtual std::shared_ptr<DissociationReaction> addDissociationReaction(
			std::shared_ptr<DissociationReaction> reaction);

	/**
	 * This operation fills an array of doubles with the concentrations of all
	 * of the reactants in the network.
	 *
	 * @param concentrations The array that will be filled with the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. If the array is too small to hold the concentrations, SIGSEGV will
	 * be thrown.
	 */
	void fillConcentrationsArray(double * concentrations);

	/**
	 * This operation updates the concentrations for all reactants in the
	 * network from an array.
	 *
	 * @param concentrations The array of doubles that will be for the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. Properly aligning the array in memory so that this operation
	 * does not overrun is up to the caller.
	 */
	virtual void updateConcentrationsFromArray(double * concentrations);

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 *
	 * Do nothing here, this method need to be implemented in
	 * subclasses.
	 *
	 * @param diagFill The pointer to the vector where the connectivity information is kept
	 */
	virtual void getDiagonalFill(int *diagFill) {
		return;
	}

	/**
	 * Get the total concentration of atoms contained in the network.
	 *
	 * Returns 0.0 here and needs to be implemented by the daughter classes.
	 *
	 * @return The total concentration
	 */
	virtual double getTotalAtomConcentration() {
		return 0.0;
	}

	/**
	 * Get the total concentration of atoms contained in bubbles in the network.
	 *
	 * Returns 0.0 here and needs to be implemented by the daughter classes.
	 *
	 * @return The total concentration
	 */
	virtual double getTotalTrappedAtomConcentration() {
		return 0.0;
	}

	/**
	 * Get the total concentration of vacancies contained in the network.
	 *
	 * Returns 0.0 here and needs to be implemented by the daughter classes.
	 *
	 * @return The total concentration
	 */
	virtual double getTotalVConcentration() {
		return 0.0;
	}

	/**
	 * Get the total concentration of material interstitials in the network.
	 *
	 * Returns 0.0 here and needs to be implemented by the daughter classes.
	 *
	 * @return The total concentration
	 */
	virtual double getTotalIConcentration() {
		return 0.0;
	}

	/**
	 * Calculate all the rate constants for the reactions and dissociations of the network.
	 * Need to be called only when the temperature changes.
	 *
	 * Need to be overwritten by daughter classes.
	 */
	virtual void computeRateConstants() {
		return;
	}

	/**
	 * Compute the fluxes generated by all the reactions
	 * for all the clusters and their momentum.
	 *
	 * Do nothing here, this method need to be implemented in
	 * subclasses.
	 *
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the fluxes are computed used to find the next solution
	 */
	virtual void computeAllFluxes(double *updatedConcOffset) {
		return;
	}

	/**
	 * Compute the partial derivatives generated by all the reactions
	 * for all the clusters and their momentum.
	 *
	 * Do nothing here, this method need to be implemented in
	 * subclasses.
	 *
	 * @param vals The pointer to the array that will contain the values of
	 * partials for the reactions
	 * @param indices The pointer to the array that will contain the indices
	 * of the clusters
	 * @param size The pointer to the array that will contain the number of reactions for
	 * this cluster
	 */
	virtual void computeAllPartials(double *vals, int *indices, int *size) {
		return;
	}

	/**
	 * This operation returns the biggest production rate in the network.
	 *
	 * @return The biggest rate
	 */
	double getBiggestRate() const {
		return biggestRate;
	}

	/**
	 * Are dissociations enabled?
	 * @returns true if reactions are enabled, false otherwise.
	 */
	bool getDissociationsEnabled() const {
		return dissociationsEnabled;
	}

	/**
	 * Enable dissociations.
	 */
	void enableDissociations() {
		dissociationsEnabled = true;
	}

	/**
	 * Disable dissociations.
	 */
	void disableDissociations() {
		dissociationsEnabled = false;
	}

	/**
	 * Number of vacancy clusters in our network.
	 */
	int getNumVClusters() const {
		return numVClusters;
	}

	/**
	 * Number of interstitial clusters in our network.
	 */
	int getNumIClusters() const {
		return numIClusters;
	}

	/**
	 * Number of super clusters in our network.
	 */
	int getNumSuperClusters() const {
		return numSuperClusters;
	}

	/**
	 * Maximum vacancy cluster size in our network.
	 */
    IReactant::SizeType getMaxVClusterSize() const {
		return maxVClusterSize;
	}

	/**
	 * Maximum interstitial cluster size in our network.
	 */
    IReactant::SizeType getMaxIClusterSize() const {
		return maxIClusterSize;
	}

};

}

#endif
