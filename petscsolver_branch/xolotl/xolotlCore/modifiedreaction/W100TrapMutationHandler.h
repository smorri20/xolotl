#ifndef W100TRAPMUTATIONHANDLER_H
#define W100TRAPMUTATIONHANDLER_H

#include "TrapMutationHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface for a (100) oriented
 * tungsten material.
 */
class W100TrapMutationHandler: public TrapMutationHandler {
private:

	/**
	 * Method initializing the depth and size vectors.
	 */
	void initializeDepthSize() {
		depthVec = {-0.1, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
		sizeVec = {2, 3, 4, 5, 6, 7, 8};

		return;
	}

public:

	/**
	 * The constructor
	 */
	W100TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~W100TrapMutationHandler() {}

};
//end class W100TrapMutationHandler

}

#endif
