#ifndef XCORE_NESUPERFLUX_H
#define XCORE_NESUPERFLUX_H

#include "Flux.h"


namespace xolotlCore {

// Flux for a NE super cluster.
struct NESuperFlux : public Flux {

    double xeMoment;

    NESuperFlux(double _total = 0,
                double _xeMoment = 0)
      : Flux(_total),
        xeMoment(_xeMoment) {

    }

    NESuperFlux(const NESuperFlux& other)
      : Flux(other.total),
        xeMoment(other.xeMoment) {

    }

    NESuperFlux& operator+=(const NESuperFlux& other) {

        Flux::operator+=(other);
        xeMoment += other.xeMoment;

        return *this;
    }

    NESuperFlux& operator-=(const NESuperFlux& other) {
    
        Flux::operator-=(other);
        xeMoment -= other.xeMoment;

        return *this;
    }
};

/**
 * Add two NESuperFluxes.
 */
inline
NESuperFlux operator+(const NESuperFlux& a, const NESuperFlux& b) {

    NESuperFlux ret(a);
    return ret += b;
}


/**
 * Subtract two NESuperFluxes.
 */
inline
NESuperFlux operator-(const NESuperFlux& a, const NESuperFlux& b) {

    NESuperFlux ret(a);
    return ret -= b;
}

} // namespace xolotlCore

#endif // XCORE_NESUPERFLUX_H
