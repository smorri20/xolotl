#ifndef PSI_CLUSTER_REACTION_NETWORK_H
#define PSI_CLUSTER_REACTION_NETWORK_H

// Includes
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "ReactionNetwork.h"
#include "PSISuperCluster.h"
#include "ReactantType.h"



// Using std::unordered_map often gives better performance than std::map,
// but requires us to define our own hash function for more complex types.
using ReactantSizePair = std::pair<xolotlCore::IReactant::SizeType, xolotlCore::IReactant::SizeType>;

namespace std {

template<>
struct hash<ReactantSizePair> {
    size_t operator()(const ReactantSizePair& pr) const {
        return (pr.first << 16) + (pr.second);
    }
};

} // namespace std


namespace xolotlCore {

/**
 *  This class manages the set of reactants and compound reactants (
 *  combinations of normal reactants) for PSI clusters. It also manages a
 *  set of properties that describes the total collection.
 *
 *  This class is a very heavyweight class that should not be abused.
 *
 *  Reactants that are added to this network must be added as with shared_ptrs.
 *  Furthermore, reactants that are added to this network have their ids set to
 *  a network specific id. Reactants should not be shared between separate
 *  instances of a PSIClusterReactionNetwork.
 */
class PSIClusterReactionNetwork: public ReactionNetwork {

private:
    /**
     * Nice name for map supporting quick lookup of supercluster containing 
     * specifc number of He and V.
     */
    using HeVToSuperClusterMap = std::unordered_map<ReactantSizePair, std::reference_wrapper<IReactant> >;


    /**
     * Map supporting quick identification of super cluster containing
     * given number of He and V.
     */
    HeVToSuperClusterMap superClusterLookupMap;

	/**
	 * Calculate the dissociation constant of the first cluster with respect to
	 * the single-species cluster of the same type based on the current clusters
	 * atomic volume, reaction rate constant, and binding energies.
	 *
	 * @param reaction The reaction
	 * @return The dissociation constant
	 */
	double calculateDissociationConstant(const DissociationReaction& reaction) const override;

	/**
	 * Calculate the binding energy for the dissociation cluster to emit the single
	 * and second cluster.
	 *
	 * @param reaction The reaction
	 * @return The binding energy corresponding to this dissociation
	 */
	double computeBindingEnergy(const DissociationReaction& reaction) const override;

	/**
     * Find the super cluster that contains the original cluster with nHe
     * helium atoms and nV vacancies.
	 *
	 * @param nHe the type of the compound reactant
	 * @param nV an array containing the sizes of each piece of the reactant.
	 * @return The super cluster representing the cluster with nHe helium
     * and nV vacancies, or nullptr if no such cluster exists.
	 */
    IReactant * getSuperFromComp(IReactant::SizeType nHe, IReactant::SizeType nV);


public:

    /**
     * Default constructor, deleted to force construction using parameters.
     */
    PSIClusterReactionNetwork() = delete;

	/**
	 * The Constructor
	 *
	 * @param registry The performance handler registry
	 */
	PSIClusterReactionNetwork(
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

    /**
     * Copy constructor, deleted to prevent use.
     */
    PSIClusterReactionNetwork(const PSIClusterReactionNetwork& other) = delete;

	/**
	 * Computes the full reaction connectivity matrix for this network.
	 */
	void createReactionConnectivity();

	/**
	 * Add the dissociation connectivity for the reverse reaction if it is allowed.
	 *
	 * @param emittingReactant The reactant that would emit the pair
	 * @param reaction The reaction we want to reverse
	 * @param a The helium number for the emitting superCluster
	 * @param b The vacancy number for the emitting superCluster
	 * @param c The helium number for the emitted superCluster
	 * @param d The vacancy number for the emitted superCluster
	 *
	 */
	void checkDissociationConnectivity(IReactant * emittingReactant,
			std::shared_ptr<ProductionReaction> reaction, int a = 0, int b = 0,
			int c = 0, int d = 0);

	/**
	 * This operation sets the temperature at which the reactants currently
	 * exists. It calls setTemperature() on each reactant.
	 *
	 * This is the simplest way to set the temperature for all reactants is to
	 * call the ReactionNetwork::setTemperature() operation.
	 *
	 * @param temp The new temperature
	 */
	virtual void setTemperature(double temp);

	/**
	 * This operation reinitializes the network.
	 *
	 * It computes the cluster Ids and network size from the allReactants vector.
	 */
	void reinitializeNetwork();

	/**
	 * This method redefines the connectivities for each cluster in the
	 * allReactans vector.
	 */
	void reinitializeConnectivities();

	/**
	 * This operation updates the concentrations for all reactants in the
	 * network from an array.
	 *
	 * @param concentrations The array of doubles that will be for the
	 * concentrations. This operation does NOT create, destroy or resize the
	 * array. Properly aligning the array in memory so that this operation
	 * does not overrun is up to the caller.
	 */
	void updateConcentrationsFromArray(double * concentrations);

	/**
	 * This operation returns the size or number of reactants and momentums in the network.
	 *
	 * @return The number of degrees of freedom
	 */
	virtual int getDOF() const override {
		return size() + 2 * getAll(ReactantType::PSISuper).size();
	}

	/**
	 * Get the diagonal fill for the Jacobian, corresponding to the reactions.
	 *
	 * @param diagFill The pointer to the vector where the connectivity information is kept
	 */
	void getDiagonalFill(int *diagFill);

	/**
	 * Get the total concentration of atoms contained in the network.
	 *
	 * Here the atoms that are considered are helium atoms.
	 *
	 * @return The total concentration
	 */
	double getTotalAtomConcentration();

	/**
	 * Get the total concentration of atoms contained in bubbles in the network.
	 *
	 * Here the atoms that are considered are helium atoms.
	 *
	 * @return The total concentration
	 */
	double getTotalTrappedAtomConcentration();

	/**
	 * Get the total concentration of vacancies contained in the network.
	 *
	 * @return The total concentration
	 */
	double getTotalVConcentration();

	/**
	 * Get the total concentration of tungsten interstitials in the network.
	 *
	 * @return The total concentration
	 */
	double getTotalIConcentration();

	/**
	 * Calculate all the rate constants for the reactions and dissociations of the network.
	 * Need to be called only when the temperature changes.
	 */
	void computeRateConstants();

	/**
	 * Compute the fluxes generated by all the reactions
	 * for all the clusters and their momentums.
	 *
	 * @param updatedConcOffset The pointer to the array of the concentration at the grid
	 * point where the fluxes are computed used to find the next solution
	 */
	void computeAllFluxes(double *updatedConcOffset);

	/**
	 * Compute the partial derivatives generated by all the reactions
	 * for all the clusters and their momentum.
	 *
	 * @param vals The pointer to the array that will contain the values of
	 * partials for the reactions
	 * @param indices The pointer to the array that will contain the indices
	 * of the clusters
	 * @param size The pointer to the array that will contain the number of reactions for
	 * this cluster
	 */
	virtual void computeAllPartials(double *vals, int *indices, int *size);

    /**
     * Construct the super cluster lookup map, keyed by number of He atoms
     * and vacancies.
     *
     * @param bounds Vector indicating boundaries of intervals to use 
     *               for num Helium and num Vacancies in super clusters.
     *               Assumed to be sorted smallest to largest, and that
     *               last element is one past the last interval's largest 
     *               allowed value.
     */
    void buildSuperClusterMap(const std::vector<IReactant::SizeType>& bounds);
};

} // namespace xolotlCore

#endif
