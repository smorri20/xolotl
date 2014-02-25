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
		 * The vector of quantities the HardwareCounter will monitor
		 */
//		std::vector<HardwareQuantities> quantities;

public:

	/**
	 * HardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param aname The HardwareCounter's name
	 * @param hquantities The HardwareCounter's list of quantities
	 */
//	GPTLHardwareCounter(std::string aname, std::vector<HardwareQuantities> hquantities);
	GPTLHardwareCounter(std::string aname);

	/**
	 * The copy constructor.
	 * @param other The GPTLHardwareCounter to copy
	 */
	GPTLHardwareCounter(const GPTLHardwareCounter &other);


    ~GPTLHardwareCounter();

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
