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
	 * @param x The position where to evaluate he fit
	 * @return The evaluated value
	 */
	double FitFunction(double x) {
		// Initial declarations
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

		p0 = 1.904709327255052;
		p1 = 3.311156742832611;

		value=p0/p1 * pow(x/p1,(p0-1.0)) * exp(-pow(x/p1,p0));

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

	/**
	 * Compute and store the incident flux values at each grid point.
	 * \see IFluxHandler.h
	 */
	void initializeFluxHandler(IReactionNetwork *network,
			int surfacePos, std::vector<double> grid) {
		// Call the general method
		FluxHandler::initializeFluxHandler(network, surfacePos, grid);

		// Set the flux index corresponding the the single helium cluster here
		auto fluxCluster = network->get(heType, 1);
		// Check that the helium cluster is present in the network
		if (!fluxCluster) {
			throw std::string(
					"\nThe single helium cluster is not present in the network, "
					"cannot use the flux option!");
		}
		fluxIndex = fluxCluster->getId() - 1;

		return;
	}

};
//end class WSRIMFitFluxHandler

}

#endif
