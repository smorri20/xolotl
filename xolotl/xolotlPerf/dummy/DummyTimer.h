#ifndef DUMMYTIMER_H
#define DUMMYTIMER_H

// Includes
#include <string>
#include "ITimer.h"

using namespace std;

namespace xolotlPerf{

/**
 * The DummyTimer class is instantiated by the DummerHandlerRegistry class
 * and realizes the DummyTimer interface.
 */
class DummyTimer : public ITimer
{
private:

	/**
	 * The name of this ITimer.
	 */
	std::string name;

	/**
	 * The default constructor is declared as private since Timers
	 *  must be initialized with a name.
	 */
	DummyTimer():ITimer("") { }

public:

	/**
	 * DummyTimer constructor that takes the argument timerName
	 * to distinguish specific DummyTimer.
	 *
	 * @param timerName The DummyTimer's name
	 */
	DummyTimer(std::string timerName);

	/**
	 * The destructor.
	 */
	~DummyTimer();

    /**
     * This operations starts the ITimer.
     */
	void start();

    /**
     * This operation stops the ITimer.
     */
	void stop();

	/**
	 * This operation returns the name of the ITimer.
	 *
	 * @return The name of this ITimer
	 */
	const std::string getName() const;

    /**
     * This operation returns the value of the DummyTimer.
     */
    double getValue();

	/**
	 * This operation returns the units of the GPTLTimer.
	 */
    long getUnits() const;

};  //end class DummyTimer

}  //end namespace xolotlPerf

#endif
