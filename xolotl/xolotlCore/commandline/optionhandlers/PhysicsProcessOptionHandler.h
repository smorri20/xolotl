#ifndef PHYSICSPROCESSOPTIONHANDLER_H
#define PHYSICSPROCESSOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"

namespace xolotlCore {

/**
 * PhysicsProcessOptionHandler handles the physics process option.
 * This option is used when a user wants to only use advection, diffusion,
 * bubble bursting, or nothing, additionally to the cluster reactions in Xolotl.
 * By default all the processes are used.
 */
class PhysicsProcessOptionHandler : public OptionHandler {
public:

	/**
	 * The default constructor
	 */
    PhysicsProcessOptionHandler() :
    	OptionHandler("process",
    			"process <process>                   "
    			"This option allows the user to only use one of the "
    			"process additionally to the reactions.  \n"
    			"                                      The options are as follows: "
    			"{diff, advec, burst, none, all}") {}

	/**
	 * The destructor
	 */
    ~PhysicsProcessOptionHandler() {}

    /**
     * This method will set the IOptions processName to the value given as the argument.
     *
     * @param opt The pointer to the option that will be modified.
     * @param arg The name of the process.
     */
    bool handler(IOptions *opt, std::string arg) {
    	// Set the process name
    	opt->setPhysicsProcess(arg);
    	return true;
    }

};//end class PhysicsProcessOptionHandler

} /* namespace xolotlCore */

#endif
