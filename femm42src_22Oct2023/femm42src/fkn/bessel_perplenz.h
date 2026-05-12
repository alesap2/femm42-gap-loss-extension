#pragma once
// bessel_perplenz.h
// Inline complex Bessel helpers for the perpendicular Lenz μ⊥(ω) model.
// Included by matprop.cpp, prob2big.cpp, prob4big.cpp.
//
// J0(z) = Σ_{n=0}^∞ (-1)^n (z/2)^{2n} / (n!)^2          [power series, 60 terms]
// J1(z) = Σ_{n=0}^∞ (-1)^n (z/2)^{2n+1} / (n! (n+1)!)
// Precise to ~1e-12 rel. for |z| ≤ 20; caller must guard |za| > 20.

#ifndef BESSEL_PERPLENZ_GUARD
#define BESSEL_PERPLENZ_GUARD

static inline CComplex cBesselJ0(CComplex z)
{
    CComplex sum  = 1.0;
    CComplex term = 1.0;
    CComplex z2o4 = -0.25 * z * z;
    for (int n = 1; n <= 60; n++) {
        term *= z2o4 / double(n * n);
        sum  += term;
    }
    return sum;
}

static inline CComplex cBesselJ1(CComplex z)
{
    CComplex term = 0.5 * z;        // n=0: z/2
    CComplex sum  = term;
    CComplex z2o4 = -0.25 * z * z;
    for (int n = 0; n < 60; n++) {
        term *= z2o4 / double((n + 1) * (n + 2));
        sum  += term;
    }
    return sum;
}

// shape(za) = 2 J₁(za) / (za J₀(za))
//   → 1 as za → 0  (no screening)
//   → 0 for |za| ≥ 20  (full Lenz screening)
static inline CComplex PerpLenzShape(CComplex za)
{
    if (abs(za) < 1e-10) return 1.0;
    if (abs(za) >= 20.0)  return 0.0;
    CComplex j0 = cBesselJ0(za);
    if (abs(j0) < 1e-30)  return 0.0;
    return 2.0 * cBesselJ1(za) / (za * j0);
}

#endif // BESSEL_PERPLENZ_GUARD
