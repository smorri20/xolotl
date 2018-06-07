#ifndef XCORE_PSISUPERFLUX_H
#define XCORE_PSISUPERFLUX_H

#include "Flux.h"


namespace xolotlCore {

// Flux for a PSI super cluster.
struct PSISuperFlux : public Flux {

    double heMoment;
    double vMoment;

    PSISuperFlux(double _total = 0,
                double _heMoment = 0,
                double _vMoment = 0)
      : Flux(_total),
        heMoment(_heMoment),
        vMoment(_vMoment) {

    }

    PSISuperFlux(const PSISuperFlux& other)
      : Flux(other.total),
        heMoment(other.heMoment),
        vMoment(other.vMoment) {

    }

    PSISuperFlux& operator+=(const PSISuperFlux& other) {

        Flux::operator+=(other);
        heMoment += other.heMoment;
        vMoment += other.vMoment;

        return *this;
    }

    PSISuperFlux& operator-=(const PSISuperFlux& other) {
    
        Flux::operator-=(other);
        heMoment -= other.heMoment;
        vMoment -= other.vMoment;

        return *this;
    }
};

/**
 * Add two PSISuperFluxes.
 */
inline
PSISuperFlux operator+(const PSISuperFlux& a, const PSISuperFlux& b) {

    PSISuperFlux ret(a);
    return ret += b;
}


/**
 * Subtract two PSISuperFluxes.
 */
inline
PSISuperFlux operator-(const PSISuperFlux& a, const PSISuperFlux& b) {

    PSISuperFlux ret(a);
    return ret -= b;
}

} // namespace xolotlCore

#endif // XCORE_PSISUPERFLUX_H
