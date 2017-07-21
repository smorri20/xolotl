// Includes
#include "HeVCluster.h"
#include "PSIClusterReactionNetwork.h"
#include <iostream>
#include <Constants.h>

using namespace xolotlCore;

HeVCluster::HeVCluster(int numHe, int numV,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(registry), numHe(numHe), numV(numV) {
	// Set the cluster size as the sum of
	// the number of Helium and Vacancies
	size = numHe + numV;

	// Update the composition map
	compositionMap[heType] = numHe;
	compositionMap[vType] = numV;

	// Set the reactant name appropriately
	std::stringstream nameStream;
	nameStream << "He_" << numHe << "V_" << numV;
	name = nameStream.str();
	// Set the typename appropriately
	typeName = heVType;

	// Compute the reaction radius
	// It is the same formula for HeV clusters
	reactionRadius = xolotlCore::tungstenLatticeConstant
			* pow((3.0 * numV) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;

	return;
}

HeVCluster::HeVCluster(HeVCluster &other) :
		PSICluster(other) {
	numHe = other.numHe;
	numV = other.numV;

	return;
}
