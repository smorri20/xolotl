#ifndef FITFLUXHANDLER_H
#define FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

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


		// Uncomment to enter the Weibull parameters for the selected surface orientation
		p[0] = 1.394689526243859;
		p[1] = 1.877732240485492;

//		// Uncomment to read the Weibull parameters from a file
//		ifstream parameters;
//		parameters.open("/home/ocekmer/Workspaces/UQTk_Xolotl/UQTk/Step2_SurrogateConstruction/FitParameters.dat");
//		double p[2];
//		for (int i=0;i<2;i++)
//			parameters >> p[i];

		// Compute the Weibull fit
		double value = p[0]/p[1] * pow(x/p[1],(p[0]-1.0)) * exp(-pow(x/p[1],p[0]));

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

};
//end class FitFluxHandler

}

#endif
