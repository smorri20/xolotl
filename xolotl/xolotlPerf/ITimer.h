#ifndef ITIMER_H
#define ITIMER_H

// Include
#include <string>

using namespace std;

namespace xolotlPerf {

/**
 * Realizations of this interface are responsible for the collection
 * of performance timing statistics.
 */
class ITimer {

public:

	/**
	 * ITimer constructor that takes the argument timerName
	 * to distinguish specific timer.
	 *
	 * @param timerName The ITimer's name
	 */
	ITimer(std::string timerName) { }

	/**
	 * The destructor
	 */
	virtual ~ITimer() { }


    /**
     * This operations starts the ITimer.
     */
    virtual void start() = 0;

    /**
     * This operation stops the ITimer.
     */
    virtual void stop() = 0;

	/**
	 * This operation returns the name of the ITimer.
	 *
	 * @return The name of this ITimer
	 */
	virtual const std::string getName() const = 0;

    /**
     * This operation returns the value of the ITimer.
     */
    virtual double getValue() = 0;

    /**
     * This operation returns the units of the ITimer.
     */
    virtual long getUnits() const = 0;

};
//end class ITimer

}//end namespace xolotlPerf

#endif
