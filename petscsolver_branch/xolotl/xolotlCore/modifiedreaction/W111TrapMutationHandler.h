#ifndef W111TRAPMUTATIONHANDLER_H
#define W111TRAPMUTATIONHANDLER_H

#include "TrapMutationHandler.h"

namespace xolotlCore {

/**
 * This class realizes the ITrapMutationHandler interface responsible for the modified
 * trap-mutation of small helium clusters close to the surface for a (111) oriented
 * tungsten material.
 */
class W111TrapMutationHandler: public TrapMutationHandler {
private:

	/**
	 * Method initializing the depth and size vectors.
	 */
	void initializeDepthSize() {
		depthVec = {-0.1, 0.9, 1.1, 1.3, 1.5};
		sizeVec = {1, 2, 3, 4, 8};

		return;
	}

public:

	/**
	 * The constructor
	 */
	W111TrapMutationHandler() {}

	/**
	 * The Destructor
	 */
	~W111TrapMutationHandler() {}

};
//end class W111TrapMutationHandler

}

#endif
