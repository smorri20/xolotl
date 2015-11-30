#ifndef DISPLACEMENTHANDLER_H
#define DISPLACEMENTHANDLER_H

#include "IDisplacementHandler.h"
#include <vector>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * flux calculations.
 */
class DisplacementHandler: public IDisplacementHandler {

protected:

	/**
	 * Vector to hold the incident flux values at each grid
	 * point (x position).
	 */
	std::vector<double> initialDisplacementVec;

	/**
	 * Step size between each grid point in the x direction.
	 */
	double stepSize;

	/**
	 * The amplitude of the flux.
	 */
	double krFluenceAmplitude;

	/**
	 * The index of the cluster.
	 */
	int displacementIndex;

	/**
	 * Value of the fit function integrated on the grid.
	 */
	double normFactor;

	/**
	 * Function that calculates the flux at a given position x (in nm).
	 * It needs to be implemented by the daughter classes.
	 *
	 * @param x The position where to evaluate the fit
	 * @return The evaluated value
	 */
	virtual double VacancyFitFunction(double x) {return 0.0;}

public:

	DisplacementHandler();

	~DisplacementHandler() {
	}

	/**
	 * Compute and store the incident flux values at each grid point.
     * \see IFluxHandler.h
	 */
	virtual void initializeDisplacementHandler(PSIClusterReactionNetwork *network,
			int nx, double hx);

	/**
	 * This operation returns the incident flux vector.
     * \see IFluxHandler.h
	 */
	virtual std::vector<double> getInitialDisplacementVec();

	/**
	 * This operation returns the index of the cluster that is irradiating
	 * the material.
     * \see IFluxHandler.h
	 */
	virtual int getInitialDisplacementClusterIndex();


	/**
	 * This operation sets the factor to change the intensity of the flux.
     * \see IFluxHandler.h
	 */
	virtual void setKrFluenceAmplitude(double krFluence);

	/**
	 * This operation gets the factor that changes the flux intensity/amplitude.
     * \see IFluxHandler.h
	 */
	virtual double getKrFluenceAmplitude() const;

};
//end class FluxHandler

}

#endif
