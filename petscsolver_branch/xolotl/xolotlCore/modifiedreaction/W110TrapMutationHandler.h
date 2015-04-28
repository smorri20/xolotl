#ifndef W110TRAPMUTATIONHANDLER_H
#define W110TRAPMUTATIONHANDLER_H

#include "TrapMutationHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface for a (110) oriented
 * tungsten material.
 */
class W110TrapMutationHandler: public TrapMutationHandler {
private:

	/**
	 * Method initializing the depth and size vectors.
	 */
	void initializeDepthSize() {
		depthVec = {-0.1, 0.5, 0.7, 0.8, 0.9, 1.0, 1.1};
		sizeVec = {2, 3, 4, 5, 6, 7, 8};

		return;
	}

public:

	/**
	 * The constructor
	 */
	W110TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~W110TrapMutationHandler() {}

};
//end class W110TrapMutationHandler

}

#endif
