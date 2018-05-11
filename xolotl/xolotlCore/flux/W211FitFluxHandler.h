#ifndef W211FITFLUXHANDLER_H
#define W211FITFLUXHANDLER_H

#include "PSIFluxHandler.h"
#include <array>

namespace xolotlCore {

/**
 * This class realizes the IFluxHandler interface to calculate the incident helium flux
 * for a (211) oriented tungsten material.
 */
class W211FitFluxHandler: public PSIFluxHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized. The surface is supposed to be (211).
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
            4.07203818, 
            5.34773722,
            -4.98297871,
			1.55833787,
            -0.234772157,
			0.0165912511,
            -2.38031874e-04,
			-3.18871642e-05,
            1.27931311e-06
        };
		return computePolynomial<double, 9>(coeffs, x);
	}

public:

	/**
	 * The constructor
	 */
	W211FitFluxHandler() {
	}

	/**
	 * The Destructor
	 */
	~W211FitFluxHandler() {
	}

};
//end class W211FitFluxHandler

}

#endif
