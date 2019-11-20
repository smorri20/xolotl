#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <xolotlPerf.h>
#include <DummyHandlerRegistry.h>
#include <AlloyCluster.h>
#include <AlloyClusterReactionNetwork.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <math.h>
#include "SimpleReactionNetwork.h"

using namespace std;
using namespace xolotlCore;
using namespace testUtils;

static std::shared_ptr<xolotlPerf::IHandlerRegistry> registry =
		std::make_shared<xolotlPerf::DummyHandlerRegistry>();

/**
 * This suite is responsible for testing the AlloyCluster.
 */
BOOST_AUTO_TEST_SUITE (AlloyCluster_testSuite)

/** This operation checks the loader. */
BOOST_AUTO_TEST_CASE(checkDiffusionCoefficient) {
	// Get the simple reaction network
	auto network = getSimpleAlloyReactionNetwork(0);
	// Create a cluster
	AlloyCluster cluster(*(network.get()), registry);
	// Add a grid point for the temperature
	cluster.addGridPoints(1);

	// Check E_m = 0.0
	cluster.setMigrationEnergy(0.0);
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), exp(0.0), 0.00001);
	BOOST_REQUIRE_CLOSE(1.0, cluster.getTemperature(0), 0.0001);

	// Make sure the diffusion coefficient is 0.0 if E_m is infinite
	cluster.setMigrationEnergy(numeric_limits<double>::infinity());
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 0.0, 0.000001);

	// Make sure the diffusion coefficient is zero if the diffusion factor is zero
	cluster.setMigrationEnergy(5.0);
	cluster.setDiffusionFactor(0.0);
	cluster.setTemperature(1.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 0.0, 0.000001);

	// Make sure the diffusion coefficient is equal to the diffusion factor
	// if the temperature is infinite
	cluster.setMigrationEnergy(5.0);
	cluster.setDiffusionFactor(1.0);
	cluster.setTemperature(numeric_limits<double>::infinity(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 1.0, 0.000001);

	// Throw something random in there to be certain
	cluster.setMigrationEnergy(0.013);
	cluster.setDiffusionFactor(1.08E10);
	cluster.setTemperature(1500.0, 0);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(0), 9766651101.800613,
			0.0000001);

	return;
}

/**
 * This operation tests the default values returned by select flux routines.
 */
BOOST_AUTO_TEST_CASE(checkDefaultFluxes) {
	// Get the simple reaction network
	auto network = getSimpleAlloyReactionNetwork(0);
	// Create a cluster
	AlloyCluster cluster(*(network.get()), registry);
	// Add a grid point for the temperature
	cluster.addGridPoints(1);

	// Check the default values of the fluxes
	BOOST_REQUIRE_CLOSE(cluster.getProductionFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getCombinationFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getDissociationFlux(0), 0.0, 1e-5);
	BOOST_REQUIRE_CLOSE(cluster.getTotalFlux(0), 0.0, 1e-5);

	return;
}

/**
 * This operation checks the different properties of the vacancy cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyVacancy) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first V reactant (numV=1)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::V, 1);

	// Check the type name
	BOOST_REQUIRE(ReactantType::V == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for V
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Faulted
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (V2)
	auto secondReactant = (AlloyCluster *) network->get(Species::V, 2);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(0.00011995, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -1.16184e-05, 0.000244496, 2.99057e-06,
			-1.73821e-06, -2.40682e-06, 0.000694105, -0.0130582, -0.0135451,
			-0.0139327, -1.25486e-06, -1.35189e-06, -1.38392e-06, -1.40886e-06,
			-1.43006e-06, -1.44894e-06, -1.46619e-06, -1.48222e-06,
			-1.49727e-06, -1.51151e-06, -1.52506e-06, -1.538e-06, -1.55041e-06,
			-1.56235e-06, 5.27689e-10, -1.47469e-06, -1.63663e-06, -1.70616e-06,
			-1.76236e-06, -1.81138e-06, -1.85584e-06, -1.89704e-06,
			-1.93571e-06, -1.97231e-06, -2.00716e-06, -2.0405e-06, -2.07249e-06,
			-2.10328e-06, -2.13299e-06, 9.60458e-10, -3.17038e+07, -2.73337e+07,
			-2.41438e+07, -2.17038e+07, -1.97715e+07, -1.81996e+07,
			-1.68933e+07, -1.57886e+07, -1.48407e+07, -1.40176e+07,
			-1.32952e+07, -1.26555e+07, -1.20845e+07, -1.15714e+07,
			-1.11075e+07, -1.06857e+07, -1.56472e-06, -1.62485e-06,
			-1.67958e-06, -1.73005e-06, -1.77702e-06, -1.82106e-06, -1.8626e-06,
			-1.90197e-06, -1.93943e-06, -1.97521e-06, -2.00947e-06,
			-2.04238e-06, -2.07405e-06, -2.10459e-06, -2.1341e-06, -2.16267e-06,
			0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of V1
	BOOST_REQUIRE_CLOSE(0.14069, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the interstitial cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyInterstitial) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first I reactant (numI=1)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::I, 1);

	// Check the type name
	BOOST_REQUIRE(ReactantType::I == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for I
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Faulted
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (I2)
	auto secondReactant = (AlloyCluster *) network->get(Species::I, 2);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-0.03778277, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { 0.000694105, -0.0130582, -0.0135451, -0.0139327,
			-0.01426, -0.125017, -0.0261139, -0.0270876, -0.0278628, -0.0145448,
			-0.0148006, -0.0150331, -0.0152469, -0.0154455, -0.0156312,
			-0.0158059, -0.0159712, -0.0161283, -0.016278, -0.0164212,
			-0.0165585, -0.0166906, -0.0168179, -0.0169407, -0.0209805,
			-0.0216873, -0.022339, -0.0229455, -0.0235141, -0.0240505,
			-0.0245588, -0.0250425, -0.0255045, -0.0259469, -0.0263718,
			-0.0267807, -0.0271751, -0.0275562, -0.027925, -3.80446e+07,
			-3.28005e+07, -2.89725e+07, -2.60446e+07, -2.37258e+07,
			-2.18396e+07, -2.02719e+07, -1.89463e+07, -1.78089e+07,
			-1.68211e+07, -1.59542e+07, -1.51866e+07, -1.45014e+07,
			-1.38857e+07, -1.3329e+07, 0, -0.0202042, -0.0209805, -0.0216873,
			-0.022339, -0.0229455, -0.0235141, -0.0240505, -0.0245588,
			-0.0250425, -0.0255045, -0.0259469, -0.0263718, -0.0267807,
			-0.0271751, -0.0275562, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of I1
	BOOST_REQUIRE_CLOSE(0.14069, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the void cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyVoid) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first Void reactant (numV=6)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::Void, 6);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Void == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Void
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Faulted
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
			// Frank
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (V1)
	auto secondReactant = (AlloyCluster *) network->get(Species::V, 1);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-7.24284e-07, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -1.35171e-06, -1.41621e-06, -1.46146e-06,
			-1.49748e-06, 1.12235e-06, -0.0145448, -0.0152388, -0.0157257,
			-0.0161133, -1.44857e-06, 2.36006e-08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.62252e+07,
			-3.11156e+07, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of Void6
	BOOST_REQUIRE_CLOSE(0.2556446, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the faulted cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyFaulted) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first Faulted reactant (numV=6)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::Faulted, 6);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Faulted == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Faulted
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Faulted
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (V1)
	auto secondReactant = (AlloyCluster *) network->get(Species::V, 1);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-8.87502e-07, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -1.62485e-06, -1.69184e-06, -1.73863e-06,
			-1.77578e-06, -1.80708e-06, -0.0209805, -0.0218456, -0.0224498,
			-0.0229294, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1.775e-06,
			4.29499e-08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4.27792e+07,
			-3.66568e+07, -3.22129e+07, -2.88311e+07, -2.61653e+07,
			-2.39717e+07, -2.21565e+07, -2.0628e+07, -1.93216e+07, -1.81913e+07,
			-1.72028e+07, -1.63303e+07, -1.5554e+07, -1.48585e+07, -1.42313e+07,
			-1.36627e+07, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of Faulted6
	BOOST_REQUIRE_CLOSE(0.431073253, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the perfect cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyPerfect) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first perfect reactant (numI=5)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::Perfect, 5);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Perfect == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Perfect
	int connectivityExpected[] = {
	// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Faulted
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			// Frank
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (I1)
	auto secondReactant = (AlloyCluster *) network->get(Species::I, 1);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-174178754, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -3.17038e+07, -3.31459e+07, -3.41551e+07,
			-3.49574e+07, -3.56342e+07, -3.80446e+07, -3.97751e+07,
			-4.09861e+07, -4.19489e+07, -3.62252e+07, -3.67534e+07,
			-3.72332e+07, -3.76742e+07, -3.80835e+07, -3.84662e+07,
			-3.88239e+07, -3.91428e+07, -3.94453e+07, -3.97333e+07,
			-4.00083e+07, -4.02719e+07, -4.05249e+07, -4.07685e+07,
			-4.10034e+07, -4.27792e+07, -4.40166e+07, -4.5161e+07, -4.62294e+07,
			-4.72341e+07, -4.81843e+07, -4.90873e+07, -4.99488e+07,
			-5.07735e+07, -5.15652e+07, -5.23272e+07, -5.30622e+07,
			-5.37726e+07, -5.44602e+07, -5.5127e+07, -6.5867e+08, -1.46335e+08,
			-1.40347e+08, -1.36122e+08, -1.33067e+08, -1.30822e+08, -1.2916e+08,
			-1.27931e+08, -1.27029e+08, -1.26381e+08, -1.25931e+08, 0, 0, 0, 0,
			0, -4.14247e+07, -4.27792e+07, -4.40166e+07, -4.5161e+07,
			-4.62294e+07, -4.72341e+07, -4.81843e+07, -4.90873e+07,
			-4.99488e+07, -5.07735e+07, -5.15652e+07, 0, 0, 0, 0, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of Perfect5
	BOOST_REQUIRE_CLOSE(0.32114234, reactant->getReactionRadius(), 0.01);

	return;
}

/**
 * This operation checks the different properties of the frank cluster
 */
BOOST_AUTO_TEST_CASE(checkAlloyFrank) {
	shared_ptr<ReactionNetwork> network = getSimpleAlloyReactionNetwork();

	// Check the reaction connectivity of the first frank reactant (numI=5)
	// Get the connectivity array from the reactant
	auto reactant = (AlloyCluster *) network->get(Species::Frank, 5);

	// Check the type name
	BOOST_REQUIRE(ReactantType::Frank == reactant->getType());
	auto reactionConnectivity = reactant->getConnectivity();

	// Check the connectivity for Frank
	int connectivityExpected[] = {
			// V
			1, 1, 1, 1, 1,
			// I
			1, 1, 1, 1,
			// Void
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Faulted
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Perfect
			1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			// Frank
			1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// Temperature
			0 };

	for (unsigned int i = 0; i < reactionConnectivity.size(); i++) {
		BOOST_REQUIRE_EQUAL(reactionConnectivity[i], connectivityExpected[i]);
	}

	// Add a grid point for the rates
	network->addGridPoints(1);
	// Get one that it combines with (I1)
	auto secondReactant = (AlloyCluster *) network->get(Species::I, 1);
	reactant->setConcentration(0.5);
	secondReactant->setConcentration(0.5);

	// Compute the rate constants that are needed for the flux
	network->setTemperature(1000.0, 0);
	// The flux can pretty much be anything except "not a number" (nan).
	double flux = reactant->getTotalFlux(0);
	BOOST_REQUIRE_CLOSE(-0.01010211, flux, 0.001);

	// The vector of partial derivatives to compare with
	double knownPartials[] = { -1.56472e-06, -1.63127e-06, -1.67779e-06,
			-1.71473e-06, -1.74588e-06, -0.0202042, -0.0210635, -0.0216641,
			0.00572162, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4.14247e+07, -3.55245e+07,
			-3.12398e+07, -2.79554e+07, -2.53459e+07, -2.32325e+07,
			-2.14832e+07, -2.00095e+07, -1.87497e+07, -1.76591e+07,
			-1.67051e+07, 0, 0, 0, 0, 0, -0.0202042, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0 };
	// Get the vector of partial derivatives
	auto partials = reactant->getPartialDerivatives(0);

	// Check the size of the partials
	BOOST_REQUIRE_EQUAL(partials.size(), 72U);

	// Check all the values
	for (unsigned int i = 0; i < partials.size(); i++) {
		BOOST_REQUIRE_CLOSE(partials[i], knownPartials[i], 0.1);
	}

	// Check the radius of Frank5
	BOOST_REQUIRE_CLOSE(0.39351424, reactant->getReactionRadius(), 0.01);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
