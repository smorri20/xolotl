/*
 * MathUtils.h
 *
 * Various math utilities.
 *
 *  Created on: May 20, 2014
 *      Author: Jay Jay Billings
 */

#ifndef MATHUTILS_H_
#define MATHUTILS_H_

#include <limits>
#include <cmath>
#include <numeric>

namespace xolotlCore {

/**
 * This function tests two doubles to see if they are equal.
 * @param a The first double
 * @param b The second double
 * @return True if the doubles are equal to within machine precision, false
 * otherwise.
 *
 * TODO careful - this test is effective only if a and b are small (e.g.,
 * less than 1).  It is not an effective test if the values are large.
 * A relative error test would be more effective.
 */
inline bool equal(double a, double b) {
	return std::fabs(b - a) < std::numeric_limits<double>::epsilon();
}

/**
 * This operation computes the Legendre polynomial, P_n (x), of degree n.
 *
 * The Legendre polynomials of degree 0 and 1 are P_0(x)=1.0 and P_1(x)=x,
 * respectively. With these conditions, the Legendre polynomials satisfy the
 * following recurrence relation: (n+1)*P_(n+1)(x) = (2n+1)*x*P_n(x) -
 * n*P_(n-1)(x)
 *
 * @param x
 *            The x value of the function
 * @param order
 *            The order of the polynomial
 */
inline double legendrePolynomial(double x, int degree) {
	// For degree 0 the return value is 1.0
	if (degree == 0)
		return 1.0;
	// For degree1 the return value i x
	if (degree == 1)
		return x;

	// Initialize the polynomials orders for the loop
	double Pn2 = 1.0, Pn1 = x, Pn = 0.0;
	// Loop on the wanted degree
	for (int n = 1; n < degree; n++) {
		// Compute the polynomial at the current order
		Pn = (((2.0 * (double) n + 1.0) * x * Pn1) - ((double) n * Pn2))
				/ ((double) n + 1.0);
		// Update the polynomials orders
		Pn2 = Pn1;
		Pn1 = Pn;
	}

	return Pn;
}

/**
 * This operation computes the Nth order Legendre polynomials
 *
 * f(x) = c0*P_0(x) + c1*P_1(x) + ... + cN*P_N(x)
 *
 * for a coefficient set {c0,c1,...,cN}.
 *
 * @param x
 *            The x value of the function.
 * @param coeffs
 *            The coefficients array.
 */
template<uint32_t N>
inline double computeNthOrderLegendre(double x, const std::array<double, N+1>& coeffs) {
    int currDegree = 0;
    auto valAtX = std::accumulate(coeffs.begin(), coeffs.end(), 0.0,
        [x,&currDegree](double running, double currCoeff) {
            return running + (currCoeff * legendrePolynomial(x, currDegree++));
        });
    return valAtX;
}

}

#endif /* MATHUTILS_H_ */
