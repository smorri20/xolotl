#ifndef FEINTERSTITIALCLUSTER_H
#define FEINTERSTITIALCLUSTER_H

// Includes
#include <sstream>
#include "FeCluster.h"
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial defects.
 */
class FeInterstitialCluster: public FeCluster {

	static std::string buildName(IReactant::SizeType nI) {
		// Set the reactant name appropriately
		std::stringstream nameStream;
		nameStream << "I_" << nI;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	FeInterstitialCluster() = delete;

	/**
	 * The constructor. All FeInterstitialClusters must be initialized with
	 * a size.
	 *
	 * @param nI The number of interstitial defect in this cluster
	 * @param registry The performance handler registry
	 */
	FeInterstitialCluster(int nI, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			FeCluster(_network, registry, buildName(nI)) {

		// Set the size
		size = nI;
		// Update the composition map
		composition[toCompIdx(Species::I)] = size;

		// Set the typename appropriately
		type = ReactantType::I;

		// Compute the reaction radius
		double EightPi = 8.0 * xolotlCore::pi;
		reactionRadius = xolotlCore::ironLatticeConstant
				* pow((3.0 / EightPi) * size, (1.0 / 3.0));

		// Bounds on He and V
		heBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));
		vBounds = IntegerRange<IReactant::SizeType>(
				static_cast<IReactant::SizeType>(0),
				static_cast<IReactant::SizeType>(1));

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	FeInterstitialCluster(const FeInterstitialCluster& other) = delete;

	/**
	 * The Destructor
	 */
	~FeInterstitialCluster() {
	}

	/**
	 * This operation returns the total change in this cluster due its
	 * own dissociation.
	 *
     * @param concs Current solution array for desired grid point.
	 * @param i The location on the grid in the depth direction
	 * @param[out] flux The flux due to its dissociation
	 */
	void getEmissionFlux(const double* concs, int i,
                            Reactant::Flux& flux) const override {

        FeCluster::getEmissionFlux(concs, i, flux);

		// Compute the loss to dislocation sinks
		if (size < 2) {
			// bias * k^2 * D * C
			flux.flux += sinkBias * sinkStrength * diffusionCoefficient[i]
					* getConcentration(concs);
		}
	}

	/**
	 * This operation computes the partial derivatives due to emission
	 * reactions.
	 *
	 * @param partials The vector into which the partial derivatives should be
	 * inserted. This vector should have a length equal to the size of the
	 * network.
	 * @param i The location on the grid in the depth direction
	 */
	void getEmissionPartialDerivatives(const double* concs, int i,
            std::vector<double> & partials) const override {
		// Initial declarations
		FeCluster::getEmissionPartialDerivatives(concs, i, partials);

		// Compute the loss to dislocation sinks
		if (size < 2) {
			// bias * k^2 * D * C
			partials[id - 1] -= sinkBias * sinkStrength * diffusionCoefficient[i];
		}

		return;
	}

};
//end class FeInterstitialCluster

} /* end namespace xolotlCore */
#endif
