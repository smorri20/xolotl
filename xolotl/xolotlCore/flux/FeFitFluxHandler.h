#ifndef FEFITFLUXHANDLER_H
#define FEFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for iron material.
 */
class FeFitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 128.0;

		if (x > x1)
			return 0.0;

		// Set fit parameters
		double a0 = -0.00073090;
		double a1 = -0.0029330;
		double b1 = 0.0039810;
		double a2 = -0.0041960;
		double b2 = 0.0084180;
		double a3 = -0.0015640;
		double b3 = 0.0099430;
		double a4 = 0.0025910;
		double b4 = 0.0063010;
		double a5 = 0.0038630;
		double b5 = 0.000880;
		double a6 = 0.0022260;
		double b6 = -0.0017580;
		double a7 = 0.00053690;
		double b7 = -0.0013570;
		double a8 = 0.00001092;
		double b8 = -0.00035910;
		double w = 0.013880;

		// Compute the fit
		double value = a0 + a1 * cos(x * w) + b1 * sin(x * w)
				+ a2 * cos(2.0 * x * w) + b2 * sin(2.0 * x * w)
				+ a3 * cos(3.0 * x * w) + b3 * sin(3.0 * x * w)
				+ a4 * cos(4.0 * x * w) + b4 * sin(4.0 * x * w)
				+ a5 * cos(5.0 * x * w) + b5 * sin(5.0 * x * w)
				+ a6 * cos(6.0 * x * w) + b6 * sin(6.0 * x * w)
				+ a7 * cos(7.0 * x * w) + b7 * sin(7.0 * x * w)
				+ a8 * cos(8.0 * x * w) + b8 * sin(8.0 * x * w);

		return std::max(value, 0.0);
	}

	/**
	 * Function that calculate the damage at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double DamageFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 128.0;

		if (x > x1)
			return 0.0;

		// Set fit parameters
		double a1 = -1.308;
		double b1 = 30.23;
		double c1 = 8.918;
		double a2 = 1.617;
		double b2 = 29.64;
		double c2 = 9.679;
		double a3 = 0.0;
		double b3 = 5.432;
		double c3 = 0.01185;
		double a4 = 0.09679;
		double b4 = 42.08;
		double c4 = 15.34;
		double a5 = 0.4024;
		double b5 = 10.19;
		double c5 = 9.01;
		double a6 = 0.383;
		double b6 = 44.74;
		double c6 = 27.35;
		double a7 = 0.08783;
		double b7 = 6.536;
		double c7 = 2.188;
		double a8 = 0.08043;
		double b8 = 18.23;
		double c8 = 4.971;

		// Compute the fit
		double value = a1 * exp(-pow((x - b1)/c1, 2.0))
		+ a2 * exp(-pow((x - b2)/c2, 2.0)) + a3 * exp(-pow((x - b3)/c3, 2.0))
		+ a4 * exp(-pow((x - b4)/c4, 2.0)) + a5 * exp(-pow((x - b5)/c5, 2.0))
		+ a6 * exp(-pow((x - b6)/c6, 2.0)) + a7 * exp(-pow((x - b7)/c7, 2.0))
		+ a8 * exp(-pow((x - b8)/c8, 2.0));

		return std::max(value, 0.0);
	}

	/**
	 * Vector to hold the damage flux values at each grid
	 * point (x position).
	 */
	std::vector<double> damageVector;

public:

	/**
	 * The constructor
	 */
	FeFitFluxHandler() {
	}

	/**
	 * The Destructor
	 */
	~FeFitFluxHandler() {
	}

	void initializeFluxHandler(PSIClusterReactionNetwork *network,
			int surfacePos, std::vector<double> grid) {
		// Set the grid
		xGrid = grid;

		// Compute the norm factor because the fit function has an
		// arbitrary amplitude
		normFactor = 0.0;
		// Loop on the x grid points skipping the first after the surface position
		// and last because of the boundary conditions
		for (int i = surfacePos + 1; i < xGrid.size() - 1; i++) {
			// Get the x position
			double x = xGrid[i] - xGrid[surfacePos];

			// Add the the value of the function times the step size
			normFactor += FitFunction(x) * (xGrid[i] - xGrid[i-1]);
		}

		// Factor the incident flux will be multiplied by to get
		// the wanted intensity
		double fluxNormalized = fluxAmplitude / normFactor;

		// Clear the flux vector
		incidentFluxVec.clear();
		// The first value corresponding to the surface position should always be 0.0
		incidentFluxVec.push_back(0.0);

		// Starts a i = surfacePos + 1 because the first value was already put in the vector
		for (int i = surfacePos + 1; i < xGrid.size() - 1; i++) {
			// Get the x position
			auto x = xGrid[i] - xGrid[surfacePos];

			// Compute the flux value
			double incidentFlux = fluxNormalized * FitFunction(x);
			// Add it to the vector
			incidentFluxVec.push_back(incidentFlux);
		}

		// The last value should always be 0.0 because of boundary conditions
		incidentFluxVec.push_back(0.0);

		// Set the flux index corresponding the the single helium cluster here
		auto fluxCluster = (PSICluster *) network->get(heType, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single helium cluster is not present in the network, "
							"cannot use the flux option!");
		}
		fluxIndex = fluxCluster->getId() - 1;

		// Compute the norm factor because the damage function has an
		// arbitrary amplitude
		double damFactor = 0.0;
		// Loop on the x grid points skipping the first after the surface position
		// and last because of the boundary conditions
		for (int i = surfacePos + 1; i < xGrid.size() - 1; i++) {
			// Get the x position
			double x = xGrid[i] - xGrid[surfacePos];

			// Add the the value of the function times the step size
			damFactor += DamageFunction(x) * (xGrid[i] - xGrid[i-1]);
		}

		// Factor the incident flux will be multiplied by to get
		// the wanted intensity
		fluxNormalized = fluxAmplitude / damFactor;

		// Clear the damage vector
		damageVector.clear();
		// The first value corresponding to the surface position should always be 0.0
		damageVector.push_back(0.0);

		// Starts a i = surfacePos + 1 because the first value was already put in the vector
		for (int i = surfacePos + 1; i < xGrid.size() - 1; i++) {
			// Get the x position
			auto x = xGrid[i] - xGrid[surfacePos];

			// Compute the flux value
			double incidentFlux = fluxNormalized * DamageFunction(x);
			// Add it to the vector
			damageVector.push_back(incidentFlux);
		}

		// The last value should always be 0.0 because of boundary conditions
		damageVector.push_back(0.0);

		return;
	}

	/**
	 * This operation returns the damage flux vector.
	 * \see IFluxHandler.h
	 */
	std::vector<double> getDamageFluxVec(double currentTime, int surfacePos) {
		return damageVector;
	}

};
//end class FeFitFluxHandler

}

#endif
