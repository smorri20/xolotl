#ifndef XCORE_FESUPERFLUX_H
#define XCORE_FESUPERFLUX_H

#include "Flux.h"


namespace xolotlCore {

// Flux for a Fe super cluster.
struct FeSuperFlux : public Flux {

    double heMoment;
    double vMoment;

    FeSuperFlux(double _total = 0,
                double _heMoment = 0,
                double _vMoment = 0)
      : Flux(_total),
        heMoment(_heMoment),
        vMoment(_vMoment) {

    }

    FeSuperFlux(const FeSuperFlux& other)
      : Flux(other.total),
        heMoment(other.heMoment),
        vMoment(other.vMoment) {

    }

    FeSuperFlux& operator+=(const FeSuperFlux& other) {

        Flux::operator+=(other);
        heMoment += other.heMoment;
        vMoment += other.vMoment;

        return *this;
    }

    FeSuperFlux& operator-=(const FeSuperFlux& other) {
    
        Flux::operator-=(other);
        heMoment -= other.heMoment;
        vMoment -= other.vMoment;

        return *this;
    }
};

/**
 * Add two FeSuperFluxes.
 */
inline
FeSuperFlux operator+(const FeSuperFlux& a, const FeSuperFlux& b) {

    FeSuperFlux ret(a);
    return ret += b;
}


/**
 * Subtract two FeSuperFluxes.
 */
inline
FeSuperFlux operator-(const FeSuperFlux& a, const FeSuperFlux& b) {

    FeSuperFlux ret(a);
    return ret -= b;
}

} // namespace xolotlCore

#endif // XCORE_FESUPERFLUX_H
