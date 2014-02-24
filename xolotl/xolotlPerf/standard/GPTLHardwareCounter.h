#ifndef GPTLHARDWARECOUNTER_H
#define GPTLHARDWARECOUNTER_H
#include "gptl.h"
#include "papi.h"
#include "HardwareCounter.h"

namespace xolotlPerf{

// The GPTLHardwareCounter class is instantiated by the StandardHandlerRegistry class
// and realizes the HardwareCounter interface to gather hardware performance counter data
// found by utilizing the PAPI (Performance Application Programming Interface) library via the
// General Purpose Timing Library (GPTL).
class GPTLHardwareCounter : public HardwareCounter
{

    public:

        ~GPTLHardwareCounter(); 

};  //end class GPTLHardwareCounter

}  //end namespace xolotlPerf

#endif
