#ifndef GPTLHARDWARECOUNTER_H
#define GPTLHARDWARECOUNTER_H

#include "gptl.h"
#include "IHardwareCounter.h"
#include "HardwareQuantities.h"
#include <unordered_map>

namespace xolotlPerf{

/**
 * The GPTLHardwareCounter class is instantiated by the StandardHandlerRegistry class
 * and realizes the IHardwareCounter interface to gather hardware performance counter data
 * found by utilizing the PAPI (Performance Application Programming Interface) library via the
 * General Purpose Timing Library (GPTL).
 */

class GPTLHardwareCounter : public IHardwareCounter
{

private:

		/**
		 * The name of this IHardwareCounter.
		 */
		std::string name;

		/**
		 * The hardware quantities this GPTLHardwareCounter monitors.
		 */
		std::vector<HardwareQuantities> quantities;

		/**
		 * The default constructor is private because HardwareCounters must
		 * always be given a name and a vector of quantities to
		 * be monitored.
		 */
		GPTLHardwareCounter():IHardwareCounter("", std::vector<HardwareQuantities>(0)) { }

protected:

		/**
		 * A map that is used to map the hardware quantities to their
		 * corresponding PAPI counterparts
		 */
//		static std::unordered_map<HardwareQuantities, int> hardwareQuantitiesMap;

public:

	/**
	 * IHardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param counterName The IHardwareCounter's name
	 * @param counterQuantities The vector of quantities the GPTLHardwareCounter will monitor
	 */
	GPTLHardwareCounter(std::string counterName, const std::vector<HardwareQuantities> &counterQuantities);

	/**
	 * The destructor
	 */
    ~GPTLHardwareCounter();

    /**
     * This operation returns a list of values of the, initially specified,
     * hardware quantities monitored by the GPTLHardwareCounter.
     */
//    std::vector<int> getValues() const;

    /**
     * This operation returns the name of the GPTLHardwareCounter.
     *
     * @return the name
     */
    const std::string getName() const;

	/**
	 * This operation returns the list of hardware
	 * quantities monitored by the GPTLHardwareCounter.
	 */
	std::vector<HardwareQuantities> getHardwareQuantities() const;


};  //end class GPTLHardwareCounter

}  //end namespace xolotlPerf

#endif
