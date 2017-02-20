#include "DisplacementHandler.h"
#include <xolotlPerf.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <mpi.h>

namespace xolotlCore {

DisplacementHandler::DisplacementHandler() :
		stepSize(0.0),
		krFluenceAmplitude(0.0),
		thresholdDisplacementEnergy(0),
		displacementIndex(-1),
		interstitialIndex(-1),
		normFactor(0.0){
	return;
}

void DisplacementHandler::initializeDisplacementHandler(PSIClusterReactionNetwork *network,
		int nx, double hx) {
	// Set the step and elementary surface sizes
	stepSize = hx;

	// Compute the norm factor because the fit function has an
	// arbitrary amplitude
	normFactor = 0.0;
	// Loop on the x grid points skipping the first and last because
	// of the boundary conditions
	for (int i = 1; i < nx - 1; i++) {
		// Get the x position
		double x = (double) i * stepSize;

		// Add the the value of the function times the step size
		normFactor += VacancyFitFunction(x) * stepSize;
	}

	// Factor the initial krypton fluence will be multiplied by to get
	// the wanted intensity
	double V_Integ = 1./3. * 3.2630862758 * 10.0;
	double krFluenceNormalized = V_Integ * krFluenceAmplitude / normFactor;

	// The first value should always be 0.0 because of boundary conditions
	initialDisplacementVec.push_back(0.0);

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = 1; i < nx - 1; i++) {
		// Get the x position
		double x = i * stepSize;

		// Compute the fluence value
		double initialDisplacement = krFluenceNormalized * VacancyFitFunction(x);
		// Add it to the vector
		initialDisplacementVec.push_back(initialDisplacement);
	}

	// The last value should always be 0.0 because of boundary conditions
	initialDisplacementVec.push_back(0.0);

	// Set the displacement index corresponding the the single vacancy cluster here
	auto displacementCluster = (PSICluster *) network->get(vType, 1);
	// Check that the helium cluster is present in the network
	if (!displacementCluster) {
		throw std::string(
				"\nThe single vacancy cluster is not present in the network, "
				"cannot use the desorption branch!");
	}
	displacementIndex = displacementCluster->getId() - 1;

	/**************************************************************************************
	 *
	 *
	 *
	 *
	 *************************************************************************************/
	// Set the step and elementary surface sizes
//	stepSize = hx;

	// Compute the norm factor because the fit function has an
	// arbitrary amplitude
	normFactor = 0.0;
	// Loop on the x grid points skipping the first and last because
	// of the boundary conditions
	for (int i = 1; i < nx - 1; i++) {
		// Get the x position
		double x = (double) i * stepSize;

		// Add the the value of the function times the step size
		normFactor += InterstitialFitFunction(x) * stepSize;
	}

	// Factor the initial krypton fluence will be multiplied by to get
	// the wanted intensity
	double I_Integ = 1./3. * 331985259.408 * 1.0e-7;

	double I_krFluenceNormalized = I_Integ * krFluenceAmplitude / normFactor;

	// The first value should always be 0.0 because of boundary conditions
	initialInterstitialVec.push_back(0.0);

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = 1; i < nx - 1; i++) {
		// Get the x position
		double x = i * stepSize;

		// Compute the fluence value
		double initialInterstitial = I_krFluenceNormalized * InterstitialFitFunction(x);

		// Add it to the vector
		initialInterstitialVec.push_back(initialInterstitial);
	}

	// The last value should always be 0.0 because of boundary conditions
	initialInterstitialVec.push_back(0.0);

	// Set the displacement index corresponding the the single vacancy cluster here
	auto interstitialCluster = (PSICluster *) network->get(iType, 1);

	// Check that the helium cluster is present in the network
	if (!interstitialCluster) {
		throw std::string(
				"\nThe single interstitial cluster is not present in the network, "
				"cannot use the desorption branch!");
	}
	interstitialIndex = interstitialCluster->getId() - 1;



	return;
}

std::vector<double> DisplacementHandler::getInitialDisplacementVec() {
	// Compute the displacement vector
	return initialDisplacementVec;
}

std::vector<double> DisplacementHandler::getInitialInterstitialVec() {
	// Compute the displacement vector
	return initialInterstitialVec;
}

int DisplacementHandler::getInitialDisplacementClusterIndex() {
	return displacementIndex;
}

int DisplacementHandler::getInitialInterstitialClusterIndex() {
	return interstitialIndex;
}

void DisplacementHandler::setKrFluenceAmplitude(double krFluence) {
	krFluenceAmplitude = krFluence;
}

double DisplacementHandler::getKrFluenceAmplitude() const {
	return krFluenceAmplitude;
}

void DisplacementHandler::setDispEnergy(int thresholdEnergy) {
	thresholdDisplacementEnergy = thresholdEnergy;
}

int DisplacementHandler::getDispEnergy() const {
	return thresholdDisplacementEnergy;
}

} // end namespace xolotlCore
