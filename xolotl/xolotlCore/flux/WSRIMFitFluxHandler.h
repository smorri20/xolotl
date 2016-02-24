#ifndef WSRIMFITFLUXHANDLER_H
#define WSRIMFITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for an amorphous tungsten material.
 */
class WSRIMFitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be amorphous.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1) return 0.0;

		double x2 = 2.0 * x/1.000100000000000122e+01-1.0;

		// Compute the fit
		double value = 0.236888
				- 0.455819 * x2
				+ 0.206163 * 1./2. * (3.0 * pow(x2, 2)-1.0)
				+ 0.222464 * 1./2. * (5.0 * pow(x2, 3)-3.0 * x2)
				- 0.426399 * 1./8. * (35.0 * pow(x2, 4)-30.0 * pow(x2, 2) + 3.0)
				+ 0.313419 * 1./8. * (63.0 * pow(x2, 5)-70.0 * pow(x2, 3) + 15.0 * x2)
				- 0.091608 * 1./16. * (231.0 * pow(x2, 6)-315.0 * pow(x2, 4) + 105.0 * (x2, 2) - 5.0);

		return value;
	}

public:

	/**
	 * The constructor
	 */
	WSRIMFitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~WSRIMFitFluxHandler() {}

};
//end class WSRIMFitFluxHandler

}

#endif
