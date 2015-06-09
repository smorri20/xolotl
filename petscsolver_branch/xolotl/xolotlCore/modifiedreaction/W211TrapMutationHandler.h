#ifndef W211TRAPMUTATIONHANDLER_H
#define W211TRAPMUTATIONHANDLER_H

#include "TrapMutationHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface for a (211) oriented
 * tungsten material.
 */
class W211TrapMutationHandler: public TrapMutationHandler {
private:

	/**
	 * Method initializing the depth vector, the size vector,
	 * and desorption information.
	 */
	void initializeDepthSize() {
		depthVec = {0.6, 0.8, 1.1, 1.1, 1.1, 1.1, 1.1};
		sizeVec = {1, 1, 1, 1, 1, 1, 2};

		// He1 desorpts with 35%
		desorp = Desorption(1, 0.35);

		return;
	}

public:

	/**
	 * The constructor
	 */
	W211TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~W211TrapMutationHandler() {}

};
//end class W211TrapMutationHandler

}

#endif
