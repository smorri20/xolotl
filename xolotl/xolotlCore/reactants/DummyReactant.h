#ifndef XCORE_DUMMYREACTANT_H
#define XCORE_DUMMYREACTANT_H

#include "Reactant.h"

namespace xolotlCore {

class DummyReactant : public Reactant {
protected:

	/**
	 * This operation returns the total change in this cluster due to
	 * other clusters dissociating into it.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] flux The flux due to dissociation of other clusters
	 */
    void getDissociationFlux(const double* concs, int i,
                                        Reactant::Flux& flux) const override {
    }

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] flux The flux due to its dissociation
	 */
    void getEmissionFlux(const double* concs, int i,
                                        Reactant::Flux& flux) const override {
    }

	/**
	 * This operation returns the total change in this cluster due to
	 * the production of this cluster by other clusters.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] flux The flux due to this cluster being produced
	 */
    void getProductionFlux(const double* concs, int i,
                                        Reactant::Flux& flux) const override {
    }

	/**
	 * This operation returns the total change in this cluster due to
	 * the combination of this cluster with others.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] flux The flux due to this cluster combining with other clusters
	 */
    void getCombinationFlux(const double* concs, int i,
                                        Reactant::Flux& flux) const override {
    }


public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	DummyReactant() = delete;

	/**
	 * The constructor.
	 *
	 * @param _network The network we will belong to.
	 * @param _name Our human-readable name.
	 * @param _registry The performance handler registry to use
	 */
	DummyReactant(IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> _registry,
			const std::string& _name = "DummyReactant")
      : Reactant(_network, _registry, _name) {
    }

    /**
     * Copy constructor.
     *
     * @param other The object to copy.
     */
    DummyReactant(const Reactant& other)
      : Reactant(other) {
    }

	/**
	 * The destructor
	 */
	virtual ~DummyReactant() {
	}

	/**
	 * This operation returns the list of partial derivatives of this reactant
	 * with respect to all other reactants in the network. The combined lists
	 * of partial derivatives from all of the reactants in the network can be
	 * used to form, for example, a Jacobian.
	 *
     * @param concs Current solution vector for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] partials The partial derivatives for this reactant where
     *  index zero corresponds to the first reactant in the list returned
     *  by the ReactionNetwork::getAll() operation.
     */
	void getPartialDerivatives(const double* concs, int i,
            std::vector<double> & partials) const override {
    }

};

} // namespace xolotlCore

#endif // XCORE_DUMMYREACTANT_H
