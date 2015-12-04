#ifndef DISPLACEMENTHANDLER_H
#define DISPLACEMENTHANDLER_H

#include "IDisplacementHandler.h"
#include <vector>
#include <memory>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for handling the initial
 * displacement calculations.
 */
class DisplacementHandler: public IDisplacementHandler {

protected:

	/**
	 * Vector to hold the initial displacement values at each grid
	 * point (x position).
	 */
	std::vector<double> initialDisplacementVec;

	/**
	 * Step size between each grid point in the x direction.
	 */
	double stepSize;

	/**
	 * The amplitude of the krypton fluence.
	 */
	double krFluenceAmplitude;

	/**
	 * The threshold displacement energy for tungsten.
	 */
	int thresholdDisplacementEnergy;

	/**
	 * The index of the cluster.
	 */
	int displacementIndex;

	/**
	 * Value of the fit function integrated on the grid.
	 */
	double normFactor;

	/**
	 * Function that calculates the displacement at a given position x (in nm).
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
	 * Compute and store the initial displacement values at each grid point.
     * \see IDisplacementHandler.h
	 */
	virtual void initializeDisplacementHandler(PSIClusterReactionNetwork *network,
			int nx, double hx);

	/**
	 * This operation returns the initial displacement vector.
     * \see IDisplacementHandler.h
	 */
	virtual std::vector<double> getInitialDisplacementVec();

	/**
	 * This operation returns the index of the vacancy cluster.
     * \see IDisplacementHandler.h
	 */
	virtual int getInitialDisplacementClusterIndex();


	/**
	 * This operation sets the factor to change the intensity of the krypton fluence.
     * \see IDisplacementHandler.h
	 */
	virtual void setKrFluenceAmplitude(double krFluence);

	/**
	 * This operation gets the factor that changes the krypton fluence intensity/amplitude.
     * \see IDisplacementHandler.h
	 */
	virtual double getKrFluenceAmplitude() const;

	/**
	 * This operation sets the factor to change the threshold energy.
     * \see IDisplacementHandler.h
	 */
	virtual void setDispEnergy(int thresholdEnergy);

	/**
	 * This operation gets the factor that changes the threshold energy.
     * \see IDisplacementHandler.h
	 */
	virtual int getDispEnergy() const;

};
//end class DisplacementHandler

}

#endif
