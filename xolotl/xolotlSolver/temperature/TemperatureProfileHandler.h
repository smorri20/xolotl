#ifndef TEMPERATUREPROFILEHANDLER_H
#define TEMPERATUREPROFILEHANDLER_H

#include "ITemperatureHandler.h"
#include <string>

namespace xolotlSolver{

/**
 * Realizations of this interface are responsible for handling the incident (incoming)
 * and outgoing Temperature calculations.
 */
class TemperatureProfileHandler: public ITemperatureHandler {

private:

	std::string tempFile;

	/**
	 * The default constructor is private because the TemperatureProfileHandler
	 * must be initialized with an input file
	 */
	TemperatureProfileHandler() :
		tempFile("")
	{ }

protected:

	std::vector<double> time;
	std::vector<double> temp;

public:

	/**
	 * The constructor
	 */
	TemperatureProfileHandler(std::string profileFileName)
		: tempFile(profileFileName)
	{ }

	/**
	 * The Destructor
	 */
	virtual ~TemperatureProfileHandler() { }

	void initializeTempData(const char* filename);

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * @param position        The position
	 * @param currentTime     The time
	 * @return temperature   The temperature
	 */
	virtual double getTemperature(std::vector<double> position, double currentTime) const;

}; //end class TemperatureProfileHandler

}

#endif
