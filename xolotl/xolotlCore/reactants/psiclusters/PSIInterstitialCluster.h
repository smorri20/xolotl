#ifndef PSIINTERSTITIALCLUSTER_H
#define PSIINTERSTITIALCLUSTER_H

// Includes
#include <sstream>
#include "PSICluster.h"
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of interstitial defects.
 */
class PSIInterstitialCluster: public PSICluster {

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
	PSIInterstitialCluster() = delete;

	/**
	 * The constructor. All PSIInterstitialClusters must be initialized with
	 * a size.
	 *
	 * @param nI The number of interstitial defect in this cluster
	 * @param registry The performance handler registry
	 */
	PSIInterstitialCluster(int nI, IReactionNetwork& _network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			PSICluster(_network, registry, buildName(nI)) {
		// Set the size
		size = nI;
		// Update the composition map
		composition[toCompIdx(Species::I)] = size;

		// Set the typename appropriately
		type = ReactantType::I;

		// Compute the reaction radius
		constexpr auto EightPi = 8 * xolotlCore::pi;
		constexpr auto aCubed = ipow<3>(xolotlCore::tungstenLatticeConstant);
		double termOne = 1.15 * (sqrt(3.0) / 4)
				* xolotlCore::tungstenLatticeConstant;
		double termTwo = std::cbrt((3 / EightPi) * aCubed * size);
		double termThree = std::cbrt((3 / EightPi) * aCubed);
		reactionRadius = termOne + termTwo - termThree;

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
	PSIInterstitialCluster(const PSIInterstitialCluster& other) = delete;

	/**
	 * The Destructor
	 */
	~PSIInterstitialCluster() {
	}

};
//end class PSIInterstitialCluster

} /* end namespace xolotlCore */
#endif
