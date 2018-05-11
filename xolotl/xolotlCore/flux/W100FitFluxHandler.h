#ifndef W100FITFLUXHANDLER_H
#define W100FITFLUXHANDLER_H

#include "PSIFluxHandler.h"
#include <array>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (100) oriented tungsten material.
 */
class W100FitFluxHandler: public PSIFluxHandler {
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
		double x1 = 10.0;

		if (x > x1)
			return 0.0;

		// Compute the fit
		std::array<double, 9> coeffs {
            7.00876507,
            0.6052078,
            -3.01711048,
			1.36595786,
            -0.295595,
			0.03597462,
            -0.0025142,
			0.0000942235,
            -0.0000014679
        };
		return computePolynomial<double, 9>(coeffs, x);
	}

public:

	/**
	 * The constructor
	 */
	W100FitFluxHandler() {
	}

	/**
	 * The Destructor
	 */
	~W100FitFluxHandler() {
	}

};
//end class W100FitFluxHandler

}

#endif
