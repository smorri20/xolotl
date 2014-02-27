#ifndef IHARDWARECOUNTER_H
#define IHARDWARECOUNTER_H

#include <string>
#include <vector>
#include <memory>
#include "HardwareQuantities.h"

using namespace std;

namespace xolotlPerf{

/**
 * Realizations of this interface are responsible for the
 * collection of hardware performance counter data.
 */
class IHardwareCounter {

private:

	/**
	 * The default constructor is private because HardwareCounters must
	 * always be given a name and a vector of quantities to
	 * be monitored.
	 */
	IHardwareCounter() { }

public:

	/**
	 * IHardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param counterName The IHardwareCounter's name
	 * @param counterQuantities The vector of quantities the IHardwareCounter will monitor
	 */
	IHardwareCounter(std::string counterName, const std::vector<HardwareQuantities> &counterQuantities) { }

	/**
	 * The destructor
	 */
	virtual ~IHardwareCounter() { }

    /**
     * This operation returns a list of values of the, initially specified,
     * hardware quantities monitored by the IHardwareCounter.
     */
//    virtual std::vector<int> getValues() const = 0;

    /**
     * This operation returns the name.
     *
     * @return the name
     */
    virtual const std::string getName() const = 0;

	/**
	 * This operation returns the list of hardware
	 * quantities monitored by the IHardwareCounter.
	 */
	virtual std::vector<HardwareQuantities> getHardwareQuantities() const = 0;

};  //end class IHardwareCounter

}  //end namespace xolotlPerf

#endif
