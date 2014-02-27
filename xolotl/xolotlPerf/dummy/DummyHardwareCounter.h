#ifndef DUMMYHARDWARECOUNTER_H
#define DUMMYHARDWARECOUNTER_H

#include <string>
#include <vector>
#include <memory>
#include "IHardwareCounter.h"
#include "HardwareQuantities.h"


namespace xolotlPerf{

/**
 * The DummyHardwareCounter class is instantiated by the DummyHandlerRegistry class
 * and realizes the DummyHardwareCounter interface.
 */
class DummyHardwareCounter : public IHardwareCounter
{

private:

	/**
	 * The name of this IHardwareCounter.
	 */
	std::string name;

	/**
	 * The hardware quantities this DummyHardwareCounter monitors.
	 */
	std::vector<HardwareQuantities> quantities;

	/**
	 * The default constructor is private because HardwareCounters
	 * must always be given a name and a vector of quantities to
	 * be monitored.
	 */
	DummyHardwareCounter():IHardwareCounter("", std::vector<HardwareQuantities>(0)) { }

public:

	/**
	 * DummyHardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param counterName The DummyHardwareCounter's name
	 * @param counterQuantities The vector of quantities the DummyHardwareCounter will monitor
	 */
	DummyHardwareCounter(std::string counterName, const std::vector<HardwareQuantities> &counterQuantities);

	/**
	 * The destructor
	 */
	~DummyHardwareCounter();

    /**
     * This operation returns a list of values of the, initially specified,
     * hardware quantities monitored by the DummyHardwareCounter.
     */
//    std::vector<int> getValues() const;

    /**
     * This operation returns the name of the DummyHardwareCounter.
     *
     * @return the name
     */
    const std::string getName() const;

	/**
	 * This operation returns the list of hardware
	 * quantities monitored by the GPTLHardwareCounter.
	 */
	std::vector<HardwareQuantities> getHardwareQuantities() const;

};  //end class DummyHardwareCounter

}  //end namespace xolotlPerf

#endif
