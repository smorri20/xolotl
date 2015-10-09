#ifndef W100FITFLUXHANDLER_H
#define W100FITFLUXHANDLER_H

#include "FluxHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (100) oriented tungsten material.
 */
class W100FitFluxHandler: public FluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (100).
	 *
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
//		// Value at which the flux goes to 0
//		double x1 = 10.0;
//
//		if (x > x1) return 0.0;
//
//		// Compute the fit
//		double value = 7.00876507 + 0.6052078 * x - 3.01711048 * pow(x, 2)
//		+ 1.36595786 * pow(x, 3) - 0.295595 * pow(x, 4) + 0.03597462 * pow(x, 5)
//		- 0.0025142 * pow(x, 6) + 0.0000942235 * pow(x, 7) - 0.0000014679 * pow(x, 8);
//
//		return value;
		if (x > 7.9 && x < 8.1) return 1.0;
		if (x > 11.9 && x < 12.1) return 1.0;
		return 0.0;
	}

public:

	/**
	 * The constructor
	 */
	W100FitFluxHandler() {}

	/**
	 * The Destructor
	 */
	~W100FitFluxHandler() {}

};
//end class W100FitFluxHandler

}

#endif
