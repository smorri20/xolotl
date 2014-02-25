#ifndef GPTLHARDWARECOUNTER_H
#define GPTLHARDWARECOUNTER_H
#include "gptl.h"
#include "papi.h"
#include "HardwareCounter.h"
#include "HardwareQuantities.h"

namespace xolotlPerf{

// The GPTLHardwareCounter class is instantiated by the StandardHandlerRegistry class
// and realizes the HardwareCounter interface to gather hardware performance counter data
// found by utilizing the PAPI (Performance Application Programming Interface) library via the
// General Purpose Timing Library (GPTL).
class GPTLHardwareCounter : public HardwareCounter
{

private:

		/**
		 * The name of this HardwareCounter.
		 */
		std::string name;

		/**
		 * The vector of hardware counter values corresponding to the PAPI
		 * quantities being monitored by the GPTLHardwareCounter
		 */
		std::vector<int> values;

		/**
		 * The vector of quantities the HardwareCounter will monitor
		 */
//		std::vector<HardwareQuantities> quantities;

public:

	/**
	 * HardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param aname The HardwareCounter's name
	 * @param hquantities The vector of quantities the GPTLHardwareCounter will monitor
	 */
	GPTLHardwareCounter(std::string aname, const std::vector<HardwareQuantities> &hquantities);

	/**
	 * The copy constructor.
	 * @param other The GPTLHardwareCounter to copy
	 */
//	GPTLHardwareCounter(const GPTLHardwareCounter &other);

	/**
	 * The destructor
	 */
    ~GPTLHardwareCounter();

    /**
     * This operation returns a list of values of the, initially specified,
     * PAPI preset quantities monitored by the GPTLHardwareCounter.
     */
    std::vector<int> getValues() const;

    /**
     * This operation returns the name of the GPTLHardwareCounter.
     * @return the name
     */
    const std::string getName() const;

    // This operation increments the GPTLHardwareCounter.
    void increment();


};  //end class GPTLHardwareCounter

}  //end namespace xolotlPerf

#endif
