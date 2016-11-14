#ifndef FITFLUXHANDLER_H
#define FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>
#include <fstream>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for tungsten material with any surface orientation.
 */
class FitFluxHandler: public FluxHandler {
private:

	/**
		 * Function that calculate the flux at a given position x (in nm).
		 * This function is not normalized. Valid for all surface orientations.
		 *
		 * @param x The position where to evaluate he fit
		 * @return The evaluated value
	*/
	double FitFunction(double x) {

		double p[2];	// Weibull parameters


//		// Uncomment to enter the Weibull parameters for the selected surface orientation
//		p[0] = 1.415022317091653;
//		p[1] = 1.918091190461391;
//		// Compute the Weibull fit
//		double value = p[0]/p[1] * pow(x/p[1],(p[0]-1.0)) * exp(-pow(x/p[1],p[0]));


		// Uncomment to read the Weibull parameters from a file
		ifstream parameters;
		parameters.open("/home/ocekmer/Workspaces/Sourceforge/uqtkXolotl/Step2_SurrogateConstruction/FitParameters.dat");  // Location of the Weibull-parameter file.
		for (int i=0;i<2;i++)
			parameters >> p[i];
		// Compute the Weibull fit
		double value = p[0]/p[1] * pow(x/p[1],(p[0]-1.0)) * exp(-pow(x/p[1],p[0]));
		parameters.close();

		return value;

	}

public:

	/**
	 * The constructor
	 */
	FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~FitFluxHandler() {}

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(IReactionNetwork *network,
			int surfacePos, std::vector<double> grid) {
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Set the flux index corresponding the the single helium cluster here
		auto fluxCluster = network->get(heType, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single helium cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndex = fluxCluster->getId() - 1;

		return;
	}

};
//end class FitFluxHandler

}

#endif
