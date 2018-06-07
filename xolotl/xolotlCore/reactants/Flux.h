#ifndef XCORE_FLUX_H
#define XCORE_FLUX_H

namespace xolotlCore {

// Basic flux at a reactant.
struct Flux {

    // Total flux at reactant.
    double total;

    // Construct a Flux value.
    Flux(double val = 0)
      : total(val) {

    }

    Flux(const Flux& other)
      : total(other.total) {

    }

    // Add a Flux to ours.
    Flux& operator+=(const Flux& other) {
        total += other.total;
        return *this;
    }

    // Subtract a Flux from ours.
    Flux& operator-=(const Flux& other) {
        total -= other.total;
        return *this;
    }

    // Scale a Flux.
    Flux& operator*=(double factor) {
        
        total *= factor;
        return *this;
    }
};

/**
 * Add two Fluxes.
 */
inline
Flux operator+(const Flux& a, const Flux& b) {
    Flux ret(a);
    return ret += b;
}


/**
 * Subtract two Fluxes.
 */
inline
Flux operator-(const Flux& a, const Flux& b) {
    Flux ret(a);
    return ret -= b;
}

/**
 * Scale a flux.
 */
inline
Flux operator*(const Flux& a, double factor) {
    Flux ret(a);
    return ret *= factor;
}

} // end namespace xolotlCore

#endif // XCORE_FLUX_H
