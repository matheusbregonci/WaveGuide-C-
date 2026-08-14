#pragma once

#include "Types.hpp"

#include <array>
#include <utility>
#include <vector>

// Abstract source of an electromagnetic field sampled on demand. Everything in
// the visualization (particle cloud, cross-section quiver/streamlines, 3D field
// lines, slice planes) only needs to (a) evaluate the field at a point, (b) know
// the domain extent, and (c) know the mode's cutoff/resonant frequency. Hiding
// the analytic rectangular/circular models behind this interface lets a future
// numerical solver (FDFD/FEM over an arbitrary voxel geometry) drop in without
// touching any of the rendering code.
//
// Coordinates passed to fieldVector/transverseField/inside are MODEL coordinates
// (the same frame the concrete models sample in): x in [0,a] or [-R,R], etc.
// The renderer centers them for display.

namespace waveguide
{

    class FieldSource
    {
    public:
        virtual ~FieldSource() = default;

        // Real 3D Cartesian field vector (E or H per fieldKind()) at time `phase`.
        virtual std::array<double, 3> fieldVector(double x, double y, double z,
                                                  double phase = 0.0) const = 0;

        // Real transverse (Vx, Vy) components of the chosen field.
        virtual std::pair<double, double> transverseField(double x, double y, double z,
                                                           double phase = 0.0) const = 0;

        // Regular-grid particle cloud of |field|, one particle per kept cell.
        virtual std::vector<Particle> sampleGrid(int nx, int ny, int nz,
                                                 bool cutawayOn = true,
                                                 float minIntensity = 0.05f,
                                                 double phase = 0.0) const = 0;

        // Axis-aligned bounding extent (meters).
        virtual Bounds bounds() const = 0;

        // Axis-aligned sampling domain in MODEL coordinates.
        //
        // bounds() gives only the SIZE, and the models disagree about where the
        // origin sits: the rectangular guide spans [0,a]x[0,b]x[0,d], while the
        // cylinder spans [-R,R]x[-R,R]x[0,L]. Anything that generates sample
        // points must ask for the real min/max instead of assuming [0,extent] --
        // feeding [0,2R] coordinates to the cylinder's inside() silently keeps
        // only the quarter disc near the corner and offsets it from the mesh.
        //
        // Default is the [0,extent] convention (rectangular, microstrip, FDTD,
        // numerical); centred models override it.
        virtual void domain(double& x0, double& x1,
                            double& y0, double& y1,
                            double& z0, double& z1) const
        {
            const Bounds b = bounds();
            x0 = 0.0; x1 = double(b.width);
            y0 = 0.0; y1 = double(b.height);
            z0 = 0.0; z1 = double(b.depth);
        }

        // Peak |field| over the domain (cached by the model; normalization ref).
        virtual double peakField() const = 0;

        virtual FieldKind fieldKind() const = 0;

        // Whether fieldVector() returns real SI units (V/m for E, A/m for H) or
        // an arbitrary amplitude. Analytic modes are eigenfunctions: their scale
        // is free until it is pinned to something physical (transported power,
        // stored energy, drive level). Sources that have not pinned it must
        // report false so the colour bar labels the axis "u.a." rather than
        // printing a number in A/m that means nothing.
        virtual bool physicalUnits() const { return false; }

        // Transverse cutoff wavenumber k_c [rad/m] of the current mode.
        virtual double cutoffWavenumber() const = 0;

        // Resonant frequency f_mnl [Hz] when closed; transverse cutoff otherwise.
        virtual double resonantFrequency() const = 0;

        virtual double epsilonRel() const = 0;
        virtual double muRel() const = 0;

        // Whether the model coordinate (x,y,z) lies inside the guide/cavity
        // (used to cull the cross-section and 3D field to the real shape).
        virtual bool inside(double x, double y, double z) const = 0;
    };

} // namespace waveguide
