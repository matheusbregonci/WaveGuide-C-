#pragma once

// Phase 1 of the frequency-domain microstrip solver: a 2D quasi-static
// electrostatic cross-section solve. It solves the Laplace equation
// div(eps * grad(phi)) = 0 over the y-z cross-section twice -- once with the
// real dielectric and once with air replacing it -- and gets the per-unit-length
// capacitances C and C_air. From those:
//
//     eps_eff = C / C_air ,   Z0 = 1 / ( c0 * sqrt(C * C_air) )
//
// which are the classic quasi-static transmission-line parameters. This validates
// the "assemble a finite-difference operator and solve a sparse system" pipeline
// on a real, verifiable result (comparable to the Hammerstad formula) before the
// full 3D complex FDFD wave solve of the later phases.

namespace waveguide
{
    struct XsecResult
    {
        double eeff = 0.0;   // effective permittivity (C / C_air)
        double Z0   = 0.0;   // characteristic impedance [ohm]
        double C    = 0.0;   // capacitance per length, with dielectric [F/m]
        double Cair = 0.0;   // capacitance per length, air only [F/m]
        int    nz   = 0;     // grid cells across the cross-section (width)
        int    ny   = 0;     // grid cells (vertical)
        int    iters = 0;    // CG iterations of the dielectric solve
        bool   ok   = false;
    };

    // W: strip width, H: substrate height (metres). epsr: substrate permittivity.
    // cellsPerH: grid resolution (cells across the substrate height).
    XsecResult solveMicrostripXsec(double W, double H, double epsr, int cellsPerH);

    // Hammerstad closed-form eps_eff / Z0 for comparison (fills eeff, Z0 only).
    XsecResult hammerstadMicrostrip(double W, double H, double epsr);
}
