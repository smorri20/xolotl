#ifndef GPTLTIMER_H
#define GPTLTIMER_H
#include "gptl.h"
#include "ITimer.h"
#include <string>

namespace xolotlPerf {

/**
 * The GPTLTimer class is instantiated by the StandardHandlerRegistry class
 * and realizes the ITimer interface to access the timing statistics found
 * General Purpose Timing Library (GPTL).
 */
class GPTLTimer: public ITimer {

private:

	/**
	 * The name of this ITimer.
	 */
	std::string name;

	/**
	 * The value of this ITimer.
	 */
	double value;

	/**
	 * The default constructor is declared private since all Timers
	 *  must be initialized with a name.
	 */
	GPTLTimer():ITimer(""), name("private"), value(0.0) { }

public:

	/**
	 * GPTLEventCounter constructor that takes the argument
	 * timerName
	 *
	 * @param timerName The GPTLEventCounter's name
	 */
	GPTLTimer(std::string timerName);

	/**
	 * The destructor
	 */
	~GPTLTimer();

    /**
     * This operations starts the ITimer.
     */
	void start();

    /**
     * This operation stops the ITimer.
     */
	void stop();

	/**
	 * This operation returns the name of the GPTLTimer.
	 *
	 * @return The name of this ITimer
	 */
	const std::string getName() const;

    /**
     * This operation returns the value of the GPTLTimer.
     */
	double getValue();

	/**
	 * This operation returns the units of the GPTLTimer.
	 *
	 * NOTE:  wall -- wallclock time (seconds)
	 * 		  usr -- user CPU time (seconds)
	 * 		  sys -- system CPU time (seconds)
	 *
	 */
	long getUnits() const;


};
//end class GPTLTimer

}//end namespace xolotlPerf

#endif
