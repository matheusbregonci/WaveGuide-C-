#pragma once

#include <cmath>

// Bessel function of the first kind, J_n(x), for integer n >= 0.
//
// Why not std::cyl_bessel_j: it is one of the C++17 "mathematical special
// functions", and libc++ does not implement them. libstdc++ (GCC) does, which
// is why the desktop build works — but the WebAssembly build goes through
// Emscripten, which uses libc++, so std::cyl_bessel_j does not compile there at
// all. Rather than #ifdef between two different numerical routines and ship a
// binary whose physics depends on the toolchain, everything uses this one.
//
// Domain that actually matters here: the circular-guide models evaluate
// J_n(k_c * rho) with n <= 3 and k_c * rho <= p_max ~ 11.62 (the largest
// tabulated Bessel zero). The ascending series is accurate and cheap over that
// range; it is NOT suitable for large arguments, where the alternating terms
// cancel catastrophically and an asymptotic form would be needed instead.

namespace waveguide {

inline double besselJ(int n, double x)
{
    if (n < 0) {
        // J_{-n} = (-1)^n J_n
        const double v = besselJ(-n, x);
        return ((-n) & 1) ? -v : v;
    }
    if (x < 0.0) {
        // J_n(-x) = (-1)^n J_n(x)
        const double v = besselJ(n, -x);
        return (n & 1) ? -v : v;
    }

    //            inf    (-1)^k          / x \ (2k+n)
    //  J_n(x) =  sum  ------------      | - |
    //            k=0   k! (n+k)!        \ 2 /
    //
    // Built by RATIO rather than by evaluating powers and factorials directly:
    // (x/2)^(2k+n) and (n+k)! both overflow long before the series is done,
    // while their quotient stays small.
    const double h = 0.5 * x;

    // term_0 = (x/2)^n / n!
    double term = 1.0;
    for (int i = 1; i <= n; ++i) term *= h / i;

    double sum = term;
    const double h2 = h * h;
    for (int k = 1; k < 80; ++k) {
        // term_k / term_{k-1} = -(x/2)^2 / (k (n+k))
        term *= -h2 / (double(k) * double(n + k));
        sum += term;
        // The terms fall off super-exponentially once k > x/2, so this exits
        // after ~x/2 + a handful of iterations.
        if (std::fabs(term) < 1e-18 * std::fabs(sum)) break;
    }
    return sum;
}

// J'_n(x). The identity J'_n = (J_{n-1} - J_{n+1})/2 needs J_{-1} = -J_1 at
// n = 0, which reduces to J'_0 = -J_1.
inline double besselJPrime(int n, double x)
{
    if (n == 0) return -besselJ(1, x);
    return 0.5 * (besselJ(n - 1, x) - besselJ(n + 1, x));
}

} // namespace waveguide
