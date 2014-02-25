#ifndef HARDWARECOUNTER_H
#define HARDWARECOUNTER_H

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
class HardwareCounter {

//private:

	/**
	 * The default constructor is private because HardwareCounters must
	 * always be initialized with a name and a vector of quantities to
	 * be monitored.
	 */
//	HardwareCounter();

public:

	/**
	 * HardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param aname The HardwareCounter's name
	 * @param hquantities The vector of quantities the HardwareCounter will monitor
	 */
	HardwareCounter(std::string aname, const std::vector<HardwareQuantities> &hquantities) {}

	/**
	 * The destructor
	 */
	virtual ~HardwareCounter() {}

    /**
     * This operation returns a list of values of the, initially specified,
     * PAPI preset quantities monitored by the HardwareCounter.
     */
    virtual std::vector<int> getValues() const = 0;

    /**
     * This operation returns the name.
     * @return the name
     */
    virtual const std::string getName() const = 0;

    // This operation increments the HardwareCounter.
    virtual void increment() = 0;

};  //end class HardwareCounter

}  //end namespace xolotlPerf

#endif
