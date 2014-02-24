#ifndef GPTLTIMER_H
#define GPTLTIMER_H
#include "gptl.h"
#include "papi.h"

#include "Timer.h"
#include <string>

namespace xolotlPerf {

/**
 * The GPTLTimer class is instantiated by the StandardHandlerRegistry class
 * and realizes the Timer interface to access the timing statistics found
 * General Purpose Timing Library (GPTL).
 */
class GPTLTimer: public Timer {

private:

	/**
	 * The name of this Timer.
	 */
	std::string name;

	/**
	 * The value of this Timer.
	 */
	double value;

	/**
	 * The default constructor is declared private since all EventCounters
	 *  must be initialized with a name.
	 */
//	GPTLTimer():Timer("") { }
//	GPTLTimer();

public:

	/**
	 * GPTLEventCounter constructor that takes the argument name
	 *
	 * @param aname The GPTLEventCounter's name
	 */
	GPTLTimer(std::string aname);

	/**
	 * The destructor
	 */
	~GPTLTimer();

	// This operations starts the Timer.
	void start();

	// This operation stops the Timer.
	void stop();

	/**
	 * This operation returns the name of the GPTLTimer.
	 *
	 * @return The name of this Timer
	 */
	const std::string getName() const;

	// This operation returns the value of the Timer.
//	long long getValue() const;
	double getValue();

	// This operation returns the units of the Timer.
	long getUnits() const;


};
//end class GPTLTimer

}//end namespace xolotlPerf

#endif
