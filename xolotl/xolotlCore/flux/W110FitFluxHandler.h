#ifndef W110FITFLUXHANDLER_H
#define W110FITFLUXHANDLER_H

#include "PSIFluxHandler.h"
#include <array>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (110) oriented tungsten material.
 */
class W110FitFluxHandler: public PSIFluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (110).
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Value at which the flux goes to 0
		double x1 = 10.0;

		if (x > x1)
			return 0.0;

		// Compute the fit
        std::array<double, 9> coeffs {
            7.93260868,
            1.49429886,
            -4.48320209,
			1.97014869,
            -0.407986353,
			0.0454535058,
            -0.0026618556,
			0.0000678768532,
            -0.000000271171991
        };
		return computePolynomial<double, 9>(coeffs, x);
	}

public:

	/**
	 * The constructor
	 */
	W110FitFluxHandler() {
	}

	/**
	 * The Destructor
	 */
	~W110FitFluxHandler() {
	}

};
//end class W110FitFluxHandler

}

#endif
