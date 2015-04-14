#ifndef BUBBLE_H
#define BUBBLE_H

// Includes
#include <string>
#include <forward_list>
#include <iostream>

namespace xolotlCore {

/**
 * A bubble is a portion of the grid that is not treated with cluster
 * dynamics anymore.
 *
 * It can expend on the grid by cumulating helium clusters and it is
 * created when the vacancy concentration at a grid point reaches the
 * tungsten density.
 */
class Bubble {
public:

	/**
	 * This is a class defining the notion of interface between the
	 * inside of the bubble and outside.
	 * It will also keep track of the quantity of helium that crossed it.
	 */
	class Interface {
		public:

			/**
			 * The inside grid point
			 */
			int insidePoint;

			/**
			 * The outside grid point
			 */
			int outsidePoint;

			/**
			 * The quantity of helium
			 */
			double heliumQuantity;

			/**
			 * The flux of helium
			 */
			double heliumFlux;

			//! The constructor
			Interface(int in, int out, double he, double heFlux = 0.0)
			: insidePoint(in), outsidePoint(out), heliumQuantity(he),
			  heliumFlux(heFlux) {}
		};

private:

	/**
	 * A list of all the grid points located inside the bubble.
	 */
	std::forward_list<int> gridPointList;

	/**
	 * A list of all the interfaces constituting the full
	 * interface of the bubble.
	 */
	std::forward_list<Interface *> interfaceList;

	/**
	 * The quantity of helium that can be contained in one cell.
	 */
	double heliumReference;

	/**
	 * The default constructor is private because bubbles
	 * must always be created with a grid point,a helium quantity,
	 * and a cell size.
	 */
	Bubble() {}

public:

	/**
	 * The constructor.
	 *
	 * @param point The initial grid point where the bubble is formed.
	 * @param he The total helium quantity
	 * @param size The size of the cell (grid point)
	 */
	Bubble(int point, double he, double size);

	/**
	 * The destructor
	 */
	~Bubble() {}

	/**
	 * Method getting the total helium quantity in the bubble.
	 *
	 * @return The helium quantity in the bubble
	 */
	double getHeliumQuantity();

	/**
	 * Method getting a pointer to the grid point list.
	 *
	 * @return The pointer to the grid point list
	 */
	std::forward_list<int> * getGridPointList();

	/**
	 * Method getting a pointer to the interface list.
	 *
	 * @return The pointer to the interface list
	 */
	std::forward_list<Interface *> * getInterfaceList();

	/**
	 * Method adding a point to the grid point list.
	 *
	 * @param The point to add
	 */
	void addGridPoint(int point);

};

} // end namespace xolotlCore

#endif
