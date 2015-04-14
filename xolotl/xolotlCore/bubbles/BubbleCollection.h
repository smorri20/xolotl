#ifndef BUBBLECOLLECTION_H
#define BUBBLECOLLECTION_H

// Includes
#include "Bubble.h"

namespace xolotlCore {

/**
 * This class gathers all the bubble that exist on the grid.
 */
class BubbleCollection {
private:

	/**
	 * The list of all the bubbles.
	 */
	std::forward_list<Bubble *> bubbleList;

public:

	/**
	 * The constructor.
	 */
	BubbleCollection() {}

	/**
	 * The destructor.
	 */
	~BubbleCollection() {}

	/**
	 * Method adding a bubble to the collection.
	 *
	 * @param bubble The bubble to add
	 */
	void addBubble(Bubble * bubble);

	/**
	 * Method getting a pointer to the bubble list.
	 *
	 * @return The pointer to the bubble list
	 */
	std::forward_list<Bubble *> * getBubbleList();

	/**
	 * Method checking if the given grid point is in a bubble.
	 *
	 * @param point The grid point to check
	 * @return True if the grid point is in a bubble
	 */
	bool isPointABubble(int point);

	/**
	 * Method getting the total helium quantity in the bubbles.
	 *
	 * @return The helium quantity in the bubbles
	 */
	double getHeliumQuantity();

	/**
	 * Method checking if two bubbles are next to each other.
	 *
	 * @param firstBubble The first bubble
	 * @param secondBubble The second bubble
	 * @return True if they are next to each other
	 */
	bool getNeighbor(Bubble * firstBubble, Bubble * secondBubble);
};

} // end namespace xolotlCore

#endif
