#ifndef TD120FITDISPLACEMENTHANDLER_H
#define TD120FITDISPLACEMENTHANDLER_H

#include "DisplacementHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IDisplacementHandler interface to calculate the initial vacancy distribution
 * in tungsten material with a threshold energy of 120 eV.
 */
class Td120FitDisplacementHandler: public DisplacementHandler {
private:

	/**
	 * Function that calculates the vacancy at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double VacancyFitFunction(double x) {
		// Value at which the vacancy goes to 0
		double x1 = 7.6;

		if (x > x1) return 0.0;

		double x2 = 2.0 * x/1.000100000000000122e+01-1.0;

		// Compute the fit
		double value = 0.016989 - 0.031020 * x2 + 0.012823 * pow(x2, 2)
		+ 0.011519 * pow(x2, 3) - 0.019257 * pow(x2, 4) + 0.013605 * pow(x2, 5)
		- 0.00511 * pow(x2, 6);

		return std::max(value,0.0);
	}

public:

	/**
	 * The constructor
	 */
	Td120FitDisplacementHandler() {}

	/**
	 * The Destructor
	 */
	~Td120FitDisplacementHandler() {}

};
//end class Td120FitDisplacementHandler

}

#endif
