#ifndef TD80FITDISPLACEMENTHANDLER_H
#define TD80FITDISPLACEMENTHANDLER_H

#include "DisplacementHandler.h"
#include <cmath>

namespace xolotlCore {

/**
 * This class realizes the IDisplacementHandler interface to calculate the initial vacancy distribution
 * in tungsten material with a threshold energy of 80 eV.
 */
class Td80FitDisplacementHandler: public DisplacementHandler {
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
		double p0, p1;
		double value;

//		// Uncomment to read the parameters from a file
//		ifstream parameters;
//		parameters.open("/home/ocekmer/Workspaces/UQTk-Xolotl/Step2_SurrogateConstruction/FitParameters.dat");
//		double p[2];
//		for (int i=0;i<2;i++)
//			parameters >> p[i];
//		z = (x-p[0])/p[1];
//		value = 1/p[1] * exp(-(z+exp(-z)));
//		parameters.close();

		p0 = 1.600418928759454;
		p1 = 2.248882336566957;

		value=p0/p1 * pow(x/p1,(p0-1.0)) * exp(-pow(x/p1,p0));

		return value;
	}

	/**
	 * Function that calculates the vacancy at a given position x (in nm).
	 * This function is not normalized.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	double InterstitialFitFunction(double x) {
		// Value at which the vacancy goes to 0
		double p0, p1;
		double value;

//		// Uncomment to read the parameters from a file
//		ifstream parameters;
//		parameters.open("/home/ocekmer/Workspaces/UQTk-Xolotl/Step2_SurrogateConstruction/FitParameters.dat");
//		double p[2];
//		for (int i=0;i<2;i++)
//			parameters >> p[i];
//		z = (x-p[0])/p[1];
//		value = 1/p[1] * exp(-(z+exp(-z)));
//		parameters.close();

		p0 = 1.664901915671198;
		p1 = 2.463576251524528;

		value=p0/p1 * pow(x/p1,(p0-1.0)) * exp(-pow(x/p1,p0));

		return value;
	}

public:

	/**
	 * The constructor
	 */
	Td80FitDisplacementHandler() {}

	/**
	 * The Destructor
	 */
	~Td80FitDisplacementHandler() {}

};
//end class Td80FitDisplacementHandler

}

#endif
