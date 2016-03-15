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
	 * Function that calculates the vacancy at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double VacancyFitFunction(double x) {
		// Value at which the vacancy goes to 0
		double x1 = 8.6;

		if (x > x1) return 0.0;

		double x2 = 2.0 * x/10.001-1.0;


		// Compute the fit
		double value = 0.303483
				- 0.528152 * x2
				+ 0.172618 * 1./2. * (3.0 * pow(x2, 2)-1.0)
				+ 0.238057 * 1./2. * (5.0 * pow(x2, 3)-3.0 * x2)
				- 0.318498 * 1./8. * (35.0 * pow(x2, 4)-30.0 * pow(x2, 2) + 3.0)
				+ 0.168097 * 1./8. * (63.0 * pow(x2, 5)-70.0 * pow(x2, 3) + 15.0 * x2)
				- 0.021070 * 1./16. * (231.0 * pow(x2, 6)-315.0 * pow(x2, 4) + 105.0 * (x2, 2) - 5.0);

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
