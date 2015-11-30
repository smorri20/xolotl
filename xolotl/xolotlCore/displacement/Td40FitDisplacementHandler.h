#ifndef TD40FITDISPLACEMENTHANDLER_H
#define TD40FITDISPLACEMENTHANDLER_H

#include "DisplacementHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IDisplacementHandler interface to calculate the initial vacancy distribution
 * in tungsten material with a threshold energy of 40 eV.
 */
class Td40FitDisplacementHandler: public DisplacementHandler {
private:

	/**
	 * Function that calculate the flux at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double VacancyFitFunction(double x) {
		// Value at which the vacancy goes to 0
		double x1 = 8.6;

		if (x > x1) return 0.0;

		double x2 = 2.0 * x/1.000100000000000122e+01-1.0;

		// Compute the fit
		double value = 0.095807 - 0.166471 * x2 + 0.054180 * pow(x2, 2)
		+ 0.075083 * pow(x2, 3) - 0.100688 * pow(x2, 4) + 0.052661 * pow(x2, 5)
		- 0.007299 * pow(x2, 6);

		return std::max(value,0.0);
	}

public:

	/**
	 * The constructor
	 */
	Td40FitDisplacementHandler() {}

	/**
	 * The Destructor
	 */
	~Td40FitDisplacementHandler() {}

};
//end class Td40FitDisplacementHandler

}

#endif
