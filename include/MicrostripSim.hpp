#pragma once

#include "FieldSource.hpp"

#include <cstdint>
#include <vector>

// Open-region FDTD for microstrip-type structures (adapted from the reference
// "microstrip code": Yee grid, per-region material). Unlike FdtdSim (a hollow
// PEC-walled metallic guide), this is an OPEN problem: a PEC ground plane + a
// dielectric substrate + a thin PEC trace sitting in air, with absorbing outer
// boundaries. Materials are resolved per Yee face so eps_r varies in space and
// PEC nodes are forced to zero.
//
// Two upgrades over the reference solver (its README's "próximos passos"):
//   * CPML (convolutional PML, Roden-Gedney) on all 6 outer faces instead of a
//     1st-order Mur ABC -- far lower spurious reflection, which is what a
//     precise S-parameter sweep needs.
//   * Optional material losses: each primitive carries an electric conductivity
//     sigma [S/m] folded into the E update (lossy-FDTD update coefficients), so
//     a lossy substrate (tan delta) or finite-conductivity fill dissipates
//     energy. sigma defaults to 0, so lossless scenes are unchanged.
//
// Vertical axis is y (ground at the bottom, substrate above it, air on top);
// the trace runs in the x-z plane. The source is a quasi-TEM modal sheet: a
// vertical Ey drive (strip -> ground voltage) over the substrate height and the
// strip width, launching a clean traveling wave down the line.

namespace waveguide
{
    class MicrostripSim : public FieldSource
    {
    public:
        enum Mat { Air = 0, Dielectric = 1, Pec = 2 };

        // A material primitive. kind 0 = axis-aligned box, 1 = z-axis cylinder
        // (a via), 2 = disk in the x-z plane extruded along y (a round trace patch
        // for the 2D designer). Coordinates in meters, in the domain frame [0,dom].
        struct Prim
        {
            int    kind = 0;
            double xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0; // box
            double cx = 0, cy = 0, cz = 0, radius = 0, rinner = 0, zlo = 0, zhi = 0; // cyl / disk (rinner>0 = ring)
            int    mat = Air;
            double epsr = 1.0;
            double sigma = 0.0; // electric conductivity [S/m]; 0 = lossless
        };

        // domX/Y/Z: domain size (m); dx: cubic cell. VERTICAL is y (ground plane
        // in the x-z plane, substrate stacked along y, air above), matching the
        // app's "up". prims: painter-order scene (last containing prim wins).
        // The source is a quasi-TEM modal sheet at x=srcX: a vertical Ey drive
        // (strip -> ground voltage) over the substrate y in [subY0,subY1] and the
        // strip width z in [srcZ0,srcZ1]; fcHz drives a continuous wave.
        MicrostripSim(double domX, double domY, double domZ, double dx,
                      const std::vector<Prim>& prims,
                      double srcX, double srcZ0, double srcZ1,
                      double subY0, double subY1, double fcHz);

        void step(int n);
        void syncField();
        void reset();
        double simTime() const { return t_; }
        long   stepCount() const { return nstep_; }

        // Pick which field the visualization samples: E (default) or H. The next
        // syncField() re-projects the chosen field to the cell-centered buffers.
        void setDisplayField(FieldKind k) { displayKind_ = k; }

        // --- Quasi-TEM transmission-line probes -----------------------------
        // Measured on a transverse plane a few cells past the launch port. The
        // strip current is I = closed-loop integral of H around the trace; the
        // port voltage is V = -integral of E through the substrate under the
        // strip. Per-unit-length L'/C' come from the local field energies, which
        // hold even under standing waves: Wm' = 1/2 L' I^2, We' = 1/2 C' V^2.
        double stripCurrent() const { return iPk_; }         // I [A]
        double portVoltage() const { return vPk_; }          // V [V]
        double inductancePerLength() const;                  // L' [H/m]
        double capacitancePerLength() const;                 // C' [F/m]
        double lineImpedance() const;                        // Z0 [ohm]

        // Surface-current magnitude |Js| ~ |H_tangential| (Hx, Hz) at a point,
        // sampled from the raw magnetic field regardless of which display field
        // (E or H) is selected. Used to colour the conductor "skin".
        double surfaceCurrent(double x, double y, double z) const;

        // Move the launch (source) and sense (measurement) ports along the line
        // (x, metres). Both are live-adjustable and clamped to the physical
        // region (outside the PML). Relocating the sense port clears its running
        // averages so L'/I re-settle at the new plane.
        void setSourcePlaneX(double x);
        void setProbePlaneX(double x);
        double sourcePlaneX() const { return (srcI_ + 0.5) * dx_; } // [m]
        double probePlaneX()  const { return iP_ * dx_; }           // [m]
        // Port markers as {x, y, z, r} in metres: centre of the launch (source) /
        // sense sheet plus a radius = half the local copper strip width. Drawn as
        // spheres at the microstrip extremities.
        std::array<double, 4> sourceMarker() const;
        std::array<double, 4> senseMarker() const;

        // Retune the launch frequency [Hz]; clears the field so the new CW
        // wave starts clean (used by the frequency sweep).
        void setFrequency(double fcHz);

        // Conformal (Dey-Mittra) PEC boundary: weight the Faraday loop by the
        // partially-filled cell so a curved trace (ring) cuts cells instead of
        // staircasing. Toggle at runtime to A/B the effect.
        void setConformal(bool on) { conformal_ = on; }
        bool conformal() const { return conformal_; }

        // ---- S-parameters (reflection / transmission) ----------------------
        // A lock-in at the drive frequency accumulates the complex voltage and
        // current phasors at two reference planes: port 1 (just past the launch)
        // and port 2 (the sense plane). Each is split into forward/backward TEM
        // waves, V+- = (V +- Z0*I)/2. With the trace running into the CPML at
        // both ends (approximately matched), S11 = V1-/V1+ (reflection) and
        // S21 = V2+/V1+ (transmission). Start a window once the field settles,
        // integrate a few periods, then read the complex S-parameters.
        void sparamStartWindow();               // reset the lock-in and integrate
        void s11(double& re, double& im) const; // reflection at port 1
        void s21(double& re, double& im) const; // transmission port 1 -> port 2

        // ---- Broadband pulse sweep (whole band in ONE run) -----------------
        // A single Gaussian pulse excites every frequency at once; a running DFT
        // of the port V/I at each requested frequency fills S11/S21 across the
        // band after one ring-down -- ~10-40x fewer time steps than driving a
        // separate CW steady state per point. Call startPulseSweep(freqs), step()
        // until pulsePast() && pulseDecay() is small, then read pulseS(m,...).
        void startPulseSweep(const std::vector<double>& freqsHz);
        int  pulseFreqCount() const { return int(pf_.size()); }
        // z0ext > 0 overrides the reference impedance for the wave split. Pass the
        // THRU line's pulseZ0(m) so a reflective filter (tiny port-2 current) is
        // decomposed against the real feed-line Z0 instead of a garbage local
        // |V2|/|I2|, which otherwise blows |S| far past 1 (non-passive).
        void pulseS(int m, double& s11r, double& s11i,
                           double& s21r, double& s21i, double z0ext = 0.0) const;
        // Characteristic impedance of THIS run's feed at frequency m (|V2|/|I2| at
        // port 2). Meaningful only when port 2 carries a clean forward wave (the
        // THRU); on a reflective filter it is unreliable -- use the THRU's value.
        double pulseZ0(int m) const;
        bool   pulsePast()  const { return pulseMode_ && t_ > 2.0 * pulseT0_; }
        double pulseDecay() const { return eNow_ / eEver_; }

        // The material primitives (ground / substrate / trace / vias) so the
        // renderer can draw the physical structure around the field.
        const std::vector<Prim>& prims() const { return prims_; }

        // ---- FieldSource ----
        std::array<double, 3> fieldVector(double x, double y, double z, double phase = 0.0) const override;
        std::pair<double, double> transverseField(double x, double y, double z, double phase = 0.0) const override;
        std::vector<Particle> sampleGrid(int nx, int ny, int nz, bool cutawayOn = true,
                                         float minIntensity = 0.05f, double phase = 0.0) const override;
        Bounds bounds() const override { return {float(domX_), float(domY_), float(domZ_)}; }
        double peakField() const override { return peak_; }
        FieldKind fieldKind() const override { return displayKind_; }
        double cutoffWavenumber() const override { return 0.0; }
        double resonantFrequency() const override { return fc_; }
        double epsilonRel() const override { return epsrSub_; }
        double muRel() const override { return 1.0; }
        bool inside(double x, double y, double z) const override;

    private:
        std::size_t idx(int i, int j, int k) const { return (std::size_t(i) * ny_ + j) * nz_ + k; }
        void stepH();
        void stepE();
        void setupCpml();
        // Build the graded (non-uniform) grid: fill the per-axis node coordinates
        // and the primary/dual edge-length metric. hf is the fine spacing; the air
        // above subTopY and the x/z margins outside the trace bbox coarsen to ~hc.
        void buildGradedGrid(double hf, double subTopY);
        // Map a physical coordinate on one axis to the floor cell index using the
        // graded node coordinates (linear search is fine at these sizes). Used by
        // the constructor's Y-index setup and sampleCell.
        static int coordToIndex(const std::vector<double>& node, double p);
        void injectSource();
        void probe();   // sample I, V and the transverse field energies
        // Port voltage & current at plane i, using a per-port transverse probe box.
        // The y voltage path (ground -> strip) and the loop's y span are shared;
        // only the z window (kc, kLoop0, kLoop1) varies per port, so port 1 and
        // port 2 can hug feed-line traces of different width / z-offset and never
        // straddle a gap.
        void planeVI(int i, int kc, int kLoop0, int kLoop1, double& V, double& I) const;
        void setupLockin();                              // phasor rotation from fc_, dt_
        double sparamZ0() const;                         // Z0 for the wave decomposition
        // Place the two reference planes on the input/output feed extremities and
        // rebuild each port's transverse probe box from the local trace geometry.
        void layoutPorts();
        // z-extent [kz0,kz1] of the trace (PEC above ground) at x=xm; false if a
        // gap (no copper) sits at that plane.
        bool traceZExtentAtX(double xm, int& kz0, int& kz1) const;
        // Nearest x-plane to i that actually carries copper (walks off a gap).
        int  snapToTrace(int i) const;
        double sampleCell(const std::vector<double>& A, double x, double y, double z) const;

        FieldKind displayKind_ = FieldKind::Electric; // field exposed to the viz

        // Transmission-line probe: plane index, voltage path (ground->strip) and
        // the Ampere loop around the trace, plus peak-held quantities.
        int iP_ = 0, jGnd_ = 0, jStrip_ = 0, kc_ = 0;
        int jLoop0_ = 0, jLoop1_ = 0, kLoop0_ = 0, kLoop1_ = 0;
        // Port-1 transverse box differs from port 2 only in z (the two feeds may
        // have different strip width / z-offset).
        int kc1_ = 0, kLoop0_1_ = 0, kLoop1_1_ = 0;
        // Feed-line x-extents (cells): the input feed ends at feedInEndI_, the
        // output feed starts at feedOutStartI_. The filter body lies between them;
        // the reference planes are kept on the feeds, off that body.
        int feedInEndI_ = 0, feedOutStartI_ = 0;
        double iPk_ = 0, vPk_ = 0, wmPk_ = 0, wePk_ = 0;

        // S-parameter lock-in: port-1 plane, phasor rotation at fc_, and the
        // accumulated V/I phasors at port 1 and port 2 (= the sense plane iP_).
        int iPort1_ = 0;
        // Fixed reference impedance for the forward/backward wave split: the feed
        // line's quasi-static (Hammerstad) Z0. Using a fixed Z0 instead of the
        // local V/I at the sense plane keeps S-params sane when a reference plane
        // is near a gap, where I -> 0 would blow up a measured Z0. 0 = fall back.
        double zRef_ = 0.0;
        bool sparamActive_ = false;
        double phc_ = 1.0, phs_ = 0.0, rotc_ = 1.0, rots_ = 0.0;
        double v1r_ = 0, v1i_ = 0, i1r_ = 0, i1i_ = 0;
        double v2r_ = 0, v2i_ = 0, i2r_ = 0, i2i_ = 0;

        // Broadband pulse sweep: Gaussian-pulse timing, per-frequency recursive
        // DFT of the port V/I (port 1 and the sense plane iP_), and decay
        // tracking (eNow_/eEver_ falls as the transient rings down).
        bool pulseMode_ = false;
        double pulseT0_ = 0, pulseTau_ = 0, pulseFc_ = 0;
        double eNow_ = 0.0, eEver_ = 1e-30;
        std::vector<double> pf_, protC_, protS_, ppc_, pps_;   // freqs + rotation
        std::vector<double> pV1r_, pV1i_, pI1r_, pI1i_;        // [freq] port-1 DFT
        std::vector<double> pV2r_, pV2i_, pI2r_, pI2i_;        // [freq] port-2 DFT

        int nx_ = 0, ny_ = 0, nz_ = 0;
        double dx_ = 0, dt_ = 0;   // dx_ = FINE reference spacing (CFL floor, defaults)
        double domX_ = 0, domY_ = 0, domZ_ = 0;

        // --- Graded (non-uniform) grid metric --------------------------------
        // Integer-node coordinates per axis (m). Cell (i,j,k) spans primary edges
        // dxp_[i] x dyp_[j] x dzp_[k]. stepH (Faraday) differentiates E across
        // primary edges -> uses idxp_/idyp_/idzp_ = 1/dxp etc. stepE (Ampere)
        // differentiates H across DUAL edges (half-node spacing) -> idxd_/idyd_/
        // idzd_ = 1/(0.5*(dxp[i]+dxp[i-1])). A uniform grid makes every entry 1/dx_.
        std::vector<double> xe_, ye_, ze_;
        std::vector<double> dxp_, dyp_, dzp_;              // primary edge lengths (probes/energy)
        std::vector<double> idxp_, idyp_, idzp_;           // 1/primary (stepH)
        std::vector<double> idxd_, idyd_, idzd_;           // 1/dual    (stepE)
        double fc_ = 0, t0_ = 0, spread_ = 0, t_ = 0;
        double epsrSub_ = 1.0;
        long nstep_ = 0;
        double peak_ = 1e-30;

        std::vector<double> Ex_, Ey_, Ez_, Hx_, Hy_, Hz_;
        std::vector<double> epsx_, epsy_, epsz_;          // kept for the energy probe
        // Precomputed lossy-FDTD update coefficients per E face: E = ca*E + cb*curl.
        // Baking these once removes ~12 divisions per cell from the stepE hot loop.
        std::vector<double> caEx_, cbEx_, caEy_, cbEy_, caEz_, cbEz_;
        std::vector<std::uint8_t> pecx_, pecy_, pecz_;

        // --- Conformal (Dey-Mittra) PEC boundary, trace layer only ------------
        // fEx_/fEz_: unblocked LENGTH fraction of the Ex / Ez edge. aHy_: unblocked
        // AREA fraction of the Hy Faraday loop face (the x-z loop that sees the
        // trace outline). 1.0 away from the copper edge, so the update reduces to
        // ordinary FDTD there. Computed only across the trace y-range (elsewhere
        // the sub-sampling would cost minutes for no benefit).
        bool conformal_ = true;
        std::vector<float> fEx_, fEz_, aHy_;

        // --- CPML (convolutional PML) on all 6 faces ------------------------
        // Per-axis 1D stretch profiles, evaluated at E-node (integer) and
        // H-node (half-integer) positions: recursion coefficients b, a and the
        // inverse stretch 1/kappa. Outside the PML slabs: a=0, 1/kappa=1.
        int npx_ = 0, npy_ = 0, npz_ = 0;                 // PML thickness per axis (cells)
        std::vector<double> beX_, aeX_, keX_, bhX_, ahX_, khX_; // x, E-node / H-node
        std::vector<double> beY_, aeY_, keY_, bhY_, ahY_, khY_; // y
        std::vector<double> beZ_, aeZ_, keZ_, bhZ_, ahZ_, khZ_; // z
        // Auxiliary convolution memory (psi) -- one per (field component,
        // derivative direction). Full-grid for clarity; nonzero only inside the
        // PML slabs, so at typical microstrip grid sizes (~1e5 cells) the extra
        // memory is a few MB. A slab-only layout would trade code for memory.
        std::vector<double> psiExy_, psiExz_, psiEyz_, psiEyx_, psiEzx_, psiEzy_;
        std::vector<double> psiHxy_, psiHxz_, psiHyz_, psiHyx_, psiHzx_, psiHzy_;

        int srcI_ = 0, srcK0_ = 0, srcK1_ = 0, srcJ0_ = 0, srcJ1_ = 0; // Ey sheet
        std::vector<Prim> prims_;          // the material scene (for drawing)

        std::vector<double> cx_, cy_, cz_; // cell-centered E for the visualization
    };

} // namespace waveguide
