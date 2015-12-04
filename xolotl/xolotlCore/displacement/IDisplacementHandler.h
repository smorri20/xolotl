#ifndef IDISPLACEMENTHANDLER_H
#define IDISPLACEMENTHANDLER_H

#include <vector>
#include <string>
#include <PSIClusterReactionNetwork.h>

namespace xolotlCore {

/**
 * Realizations of this interface are responsible for handling the initial
 * displacement calculations.
 */
class IDisplacementHandler {

public:

	virtual ~IDisplacementHandler() { }

	/**
	 * Compute and store the initial vacancy values at each grid point.
	 *
	 * @param network The reaction network
	 * @param nx The total number of grid points that will be used on the x axis
	 * @param hx The step size between grid points on the x axis
	 */
	virtual void initializeDisplacementHandler(PSIClusterReactionNetwork *network,
			int nx, double hx) = 0;

	/**
	 * This operation returns the initial displacement vector.
	 *
	 * @param currentTime The time
	 * @return The incident flux vector
	 */
	virtual std::vector<double> getInitialDisplacementVec() = 0;

	/**
	 * This operation returns the index of the vacancy cluster.
	 *
	 * @return The index of the incident flux cluster
	 */
	virtual int getInitialDisplacementClusterIndex() = 0;

	/**
	 * This operation sets the factor to change the intensity of the krypton fluence.
	 *
	 * @param flux The flux intensity
	 */
	virtual void setKrFluenceAmplitude(double krFluence) = 0;

	/**
	 * This operation gets the factor to change the intensity of the krypton fluence.
	 *
	 * @return The flux amplitude
	 */
	virtual double getKrFluenceAmplitude() const = 0;

	/**
	 * This operation sets the factor to change the threshold energy.
	 *
	 * @param flux The flux intensity
	 */
	virtual void setDispEnergy(int thresholdEnergy) = 0;

	/**
	 * This operation gets the factor that changes the threshold energy.
	 * intensity/amplitude.
	 *
	 * @return The flux amplitude
	 */
	virtual int getDispEnergy() const = 0;

}; //end class IFluxHandler

}

#endif
