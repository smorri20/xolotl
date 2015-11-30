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
		displacementIndex(-1),
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

	// Factor the incident flux will be multiplied by to get
	// the wanted intensity
	double krFluenceNormalized = krFluenceAmplitude / normFactor;

	// The first value should always be 0.0 because of boundary conditions
	initialDisplacementVec.push_back(0.0);

	// Starts a i = 1 because the first value was already put in the vector
	for (int i = 1; i < nx - 1; i++) {
		// Get the x position
		double x = i * stepSize;

		// Compute the flux value
		double initialDisplacement = krFluenceNormalized * VacancyFitFunction(x);
		// Add it to the vector
		initialDisplacementVec.push_back(initialDisplacement);
	}

	// The last value should always be 0.0 because of boundary conditions
	initialDisplacementVec.push_back(0.0);

	// Set the flux index corresponding the the single helium cluster here
	auto displacementCluster = (PSICluster *) network->get(vType, 1);
	// Check that the helium cluster is present in the network
	if (!displacementCluster) {
		throw std::string(
				"\nThe single helium cluster is not present in the network, "
				"cannot use the flux option!");
	}
	displacementIndex = displacementCluster->getId() - 1;

	return;
}

std::vector<double> DisplacementHandler::getInitialDisplacementVec() {
	// Recompute the flux vector if a time profile is used
	return initialDisplacementVec;
}

int DisplacementHandler::getInitialDisplacementClusterIndex() {
	return displacementIndex;
}

void DisplacementHandler::setKrFluenceAmplitude(double krFluence) {
	krFluenceAmplitude = krFluence;
}

double DisplacementHandler::getKrFluenceAmplitude() const {
	return krFluenceAmplitude;
}

} // end namespace xolotlCore
