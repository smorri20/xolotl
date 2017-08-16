#include "TCluster.h"

using namespace xolotlCore;

TCluster::TCluster(int nT,
        IReactionNetwork& _network,
        std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		HCluster(nT, _network, registry) {
	// Set the reactant name appropriately
	name = "Tritium";
}
TCluster::~TCluster() {
}
