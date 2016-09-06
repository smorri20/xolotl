// Includes
#include "Reactant.h"
#include "ReactionNetwork.h"
#include <xolotlPerf.h>
#include <iostream>
#include <math.h>
#include <MathUtils.h>
#include <Constants.h>

// Namespaces
using namespace xolotlCore;

Reactant::Reactant() :
		concentration(0.0), id(0), temperature(0.0), typeName(""), network(
				nullptr), handlerRegistry(nullptr), size(0), formationEnergy(
				0.0), diffusionFactor(0.0), diffusionCoefficient(0.0), migrationEnergy(
				0.0), name("Reactant"), reactionRadius(0.0), thisNetworkIndex(
				-1), biggestRate(0.0) {
	// Setup the composition map.
	compositionMap[heType] = 0;
	compositionMap[vType] = 0;
	compositionMap[iType] = 0;
}

Reactant::Reactant(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		concentration(0.0), id(0), temperature(0.0), typeName(""), network(
				nullptr), handlerRegistry(registry), size(0), formationEnergy(
				0.0), diffusionFactor(0.0), diffusionCoefficient(0.0), migrationEnergy(
				0.0), name("Reactant"), reactionRadius(0.0), thisNetworkIndex(
				-1), biggestRate(0.0) {
	// Setup the composition map.
	compositionMap[heType] = 0;
	compositionMap[vType] = 0;
	compositionMap[iType] = 0;
}

Reactant::Reactant(Reactant &other) :
		concentration(other.concentration), name(other.name), typeName(
				other.typeName), id(other.id), temperature(other.temperature), network(
				other.network), compositionMap(other.compositionMap), handlerRegistry(
				other.handlerRegistry), size(other.size), formationEnergy(
				other.formationEnergy), diffusionFactor(other.diffusionFactor), diffusionCoefficient(
				other.diffusionCoefficient), migrationEnergy(
				other.migrationEnergy), reactionRadius(other.reactionRadius), thisNetworkIndex(
				other.thisNetworkIndex), biggestRate(other.biggestRate), reactionConnectivitySet(
				other.reactionConnectivitySet), dissociationConnectivitySet(
				other.dissociationConnectivitySet) {
	// Setup the composition map.
	compositionMap[heType] = other.compositionMap[heType];
	compositionMap[vType] = other.compositionMap[vType];
	compositionMap[iType] = other.compositionMap[iType];
}

void Reactant::recomputeDiffusionCoefficient(double temp) {
	// Return zero if the diffusion factor is zero.
	if (xolotlCore::equal(diffusionFactor, 0.0)) {
		diffusionCoefficient = 0.0;
	} else {
		// Otherwise use the Arrhenius equation to compute the diffusion
		// coefficient
		double k_b = xolotlCore::kBoltzmann;
		double kernel = -1.0 * migrationEnergy / (k_b * temp);
		diffusionCoefficient = diffusionFactor * exp(kernel);
	}

	return;
}

double Reactant::getConcentration(double distHe, double distV) const {
	return concentration;
}

double Reactant::getHeMomentum() const {
	return 0.0;
}

double Reactant::getVMomentum() const {
	return 0.0;
}

void Reactant::setConcentration(double conc) {
	concentration = conc;

	return;
}

void Reactant::setReactionNetwork(
		std::shared_ptr<IReactionNetwork> reactionNetwork) {
	network = reactionNetwork;

	return;
}

void Reactant::releaseReactionNetwork() {
	network.reset();

	return;
}

void Reactant::setReactionConnectivity(int id) {
	// Add the cluster to the set.
	reactionConnectivitySet.insert(id);

	return;
}

void Reactant::setZerothMomentum(double mom) {
	return;
}

void Reactant::setHeMomentum(double mom) {
	return;
}

void Reactant::setVMomentum(double mom) {
	return;
}

double Reactant::getTotalFlux() {
	return 0.0;
}

double Reactant::getHeMomentFlux() const {
	return 0.0;
}

double Reactant::getVMomentFlux() const {
	return 0.0;
}

void Reactant::setDissociationConnectivity(int id) {
	// Add the cluster to the set.
	dissociationConnectivitySet.insert(id);

	return;
}

std::vector<int> Reactant::getConnectivity() const {
	// The connectivity array by default is filled with
	// zeros.
	int connectivityLength = network->size();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);

	// This reactant should be connected to itself
	connectivity[id - 1] = 1;

	return connectivity;
}

std::vector<double> Reactant::getPartialDerivatives() const {
	// The partial derivatives array by default is filled with
	// zeros.
	int length = network->size();
	std::vector<double> partial_derivatives = std::vector<double>(length, 0.0);

	return partial_derivatives;
}

void Reactant::getPartialDerivatives(std::vector<double> &) const {
	// Do nothing.
	return;
}

const std::string Reactant::getName() const {
	return name;
}

std::string Reactant::getType() const {
	return typeName;
}

const std::map<std::string, int> & Reactant::getComposition() const {
	// Return the composition map
	return compositionMap;
}

void Reactant::setId(int nId) {
	id = nId;

	return;
}

int Reactant::getId() const {
	return id;
}

void Reactant::setTemperature(double temp) {
	temperature = temp;

	// Recompute the diffusion coefficient
	recomputeDiffusionCoefficient(temp);
}

double Reactant::getTemperature() const {
	return temperature;
}

int Reactant::getSize() const {
	// Return this cluster's size
	return size;
}

double Reactant::getFormationEnergy() const {
	return formationEnergy;
}

void Reactant::setFormationEnergy(double energy) {
	// Set the formation energy
	formationEnergy = energy;
	return;
}

double Reactant::getDiffusionFactor() const {
	// Return the diffusion factor
	return diffusionFactor;
}

void Reactant::setDiffusionFactor(const double factor) {
	// Set the diffusion factor
	diffusionFactor = factor;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

double Reactant::getDiffusionCoefficient() const {
	return diffusionCoefficient;
}

double Reactant::getMigrationEnergy() const {
	// Return the migration energy
	return migrationEnergy;
}

void Reactant::setMigrationEnergy(const double energy) {
	// Set the migration energy
	migrationEnergy = energy;
	// Update the diffusion coefficient
	recomputeDiffusionCoefficient(temperature);

	return;
}

double Reactant::getReactionRadius() const {
	return reactionRadius; // Computed by subclasses in constructors.
}

double Reactant::getBiggestRate() const {
	return biggestRate; // Computed by computeRateConstants
}
