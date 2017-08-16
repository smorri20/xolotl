#include "DCluster.h"

using namespace xolotlCore;

DCluster::DCluster(int nH,
                    IReactionNetwork& _network,
                    std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		HCluster(nH, _network, registry) {
	// Set the reactant name appropriately
	name = "Deuterium";
}
DCluster::~DCluster() {
}
