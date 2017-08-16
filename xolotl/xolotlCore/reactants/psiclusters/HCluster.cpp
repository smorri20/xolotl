// Includes
#include "HCluster.h"

using namespace xolotlCore;

HCluster::HCluster(int nH,
                    IReactionNetwork& _network,
                    std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		PSICluster(_network, registry) {
	// Set the size appropriately
	size = nH;
	// Set the reactant name appropriately
	name = "Hydrogen";
}

HCluster::~HCluster() { }
