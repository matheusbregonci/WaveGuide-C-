#pragma once

#include "FieldSource.hpp"
#include "Geometry.hpp" // VoxelMask

#include <cstdint>
#include <vector>

// Time-domain (FDTD) electromagnetic simulation on an arbitrary voxel geometry
// (Phase 5). A wave is launched from an input port and propagates through the
// circuit; PEC walls come from the mask and graded absorbing layers at the
// ports let the wave leave instead of reflecting.
//
// Ports are the real OPENINGS of the solid (where an arm end meets the bounding
// box), each a rectangular patch on one bounding-box face, not the whole face.

namespace waveguide
{
    // A port opening: a transverse rectangle on one bounding-box face.
    //   axis: 0=x, 1=y, 2=z (face normal);  side: 0=min, 1=max
    //   role: 0=wall (inactive), 1=input (source+absorb), 2=output (absorb)
    //   [uMin,uMax]x[vMin,vMax]: transverse extent in MODEL coords (meters).
    //   The transverse axes are (y,z) for x-normal, (x,z) for y-normal,
    //   (x,y) for z-normal.
    struct FdtdPort
    {
        int axis = 2, side = 0, role = 0;
        double uMin = 0, uMax = 0, vMin = 0, vMax = 0;
    };

    class FdtdSim : public FieldSource
    {
    public:
        // A broadband Gaussian pulse covering [fminHz, fmaxHz] is launched; the
        // running DFT at the ports yields S-parameters over that band.
        FdtdSim(const VoxelMask& mask, double fminHz, double fmaxHz,
                const std::vector<FdtdPort>& ports,
                double epsR = 1.0, double muR = 1.0);

        void step(int n);
        void syncField();
        void reset();
        double simTime() const { return t_; }
        long stepCount() const { return nstep_; }

        // ---- Port monitors ----
        static const int kHist = 2000;
        int portCount() const { return int(ports_.size()); }
        const FdtdPort& port(int i) const { return ports_[i]; }
        const std::vector<float>& portHistory(int i) const { return hist_[i]; }
        int portHistoryPos(int i) const { return histPos_[i]; }
        double portRms(int i) const;

        // ---- S-parameter spectra (running DFT) ----
        int numFreqs() const { return int(freqs_.size()); }
        double freqHz(int m) const { return (m >= 0 && m < int(freqs_.size())) ? freqs_[m] : 0.0; }
        void portSpectrum(int p, int m, double& re, double& im) const;
        void sourceSpectrum(int m, double& re, double& im) const;

        // Sweep convergence: the DFT is only valid once the excitation pulse has
        // passed and the field has rung down. decayRatio() is the present field
        // amplitude over the largest seen (so < ~1e-3 means the transient is
        // essentially over and the S-curve has settled).
        bool pulsePast() const { return t_ > 2.0 * t0_; }
        double decayRatio() const { return emaxNow_ / emaxEver_; }

        // ---- CW steady-state sweep (VNA-style, one frequency at a time) ----
        // Switch to a continuous sinusoid at fHz, reset the fields, and start a
        // fresh steady-state measurement. Step() then drives CW and runs a
        // per-port lock-in in windows of a few periods; cwSteady() flips true
        // when two consecutive windows agree (the field reached steady state).
        void setCW(double fHz);
        bool cwSteady() const { return cwSteady_; }
        double cwResidual() const { return cwResidual_; }
        int cwBlocks() const { return cwBlocks_; }
        void cwPortAmp(int p, double& re, double& im) const;
        void cwSourceAmp(double& re, double& im) const;

        // ---- FieldSource ----
        std::array<double, 3> fieldVector(double x, double y, double z, double phase = 0.0) const override;
        std::pair<double, double> transverseField(double x, double y, double z, double phase = 0.0) const override;
        std::vector<Particle> sampleGrid(int nx, int ny, int nz, bool cutawayOn = true,
                                         float minIntensity = 0.05f, double phase = 0.0) const override;
        Bounds bounds() const override { return {float(W_), float(H_), float(D_)}; }
        double peakField() const override { return peak_; }
        FieldKind fieldKind() const override { return FieldKind::Electric; }
        double cutoffWavenumber() const override { return 0.0; }
        double resonantFrequency() const override { return f0_; }
        double epsilonRel() const override { return epsR_; }
        double muRel() const override { return muR_; }
        bool inside(double x, double y, double z) const override;

    private:
        int iEx(int i, int j, int k) const { return (k * (ny_ + 1) + j) * nx_ + i; }
        int iEy(int i, int j, int k) const { return (k * ny_ + j) * (nx_ + 1) + i; }
        int iEz(int i, int j, int k) const { return (k * (ny_ + 1) + j) * (nx_ + 1) + i; }
        int iHx(int i, int j, int k) const { return (k * ny_ + j) * (nx_ + 1) + i; }
        int iHy(int i, int j, int k) const { return (k * (ny_ + 1) + j) * nx_ + i; }
        int iHz(int i, int j, int k) const { return (k * ny_ + j) * nx_ + i; }
        std::size_t gidx(int i, int j, int k) const { return (std::size_t(k) * ny_ + j) * nx_ + i; }
        bool solid(int i, int j, int k) const;
        double sampleCell(const std::vector<double>& A, double ux, double uy, double uz) const;

        int nx_ = 0, ny_ = 0, nz_ = 0;
        double dx_ = 0, dy_ = 0, dz_ = 0, dt_ = 0;
        double W_ = 0, H_ = 0, D_ = 0;
        double f0_ = 0, t0_ = 0, tau_ = 0, t_ = 0;
        long nstep_ = 0;
        double peak_ = 1e-30;
        double emaxNow_ = 0.0, emaxEver_ = 1e-30; // field-decay tracking
        double epsR_ = 1.0, muR_ = 1.0;           // uniform filling medium
        std::vector<std::uint8_t> occ_;

        std::vector<double> Ex_, Ey_, Ez_, Hx_, Hy_, Hz_;
        std::vector<std::uint8_t> exFree_, eyFree_, ezFree_;
        std::vector<double> exLoss_, eyLoss_, ezLoss_;

        std::vector<int> srcComp_, srcIdx_;
        std::vector<double> srcAmp_;

        std::vector<FdtdPort> ports_;
        struct Monitor { int comp = 1; std::vector<int> idx; std::vector<float> w; };
        std::vector<Monitor> mon_;
        std::vector<std::vector<float>> hist_;
        std::vector<int> histPos_;
        std::vector<double> portSig_;   // this step's per-port modal signal

        // Running (recursive) DFT of each port signal and of the source drive,
        // evaluated at freqs_ over the sweep band. pc_/ps_ hold cos/sin(2*pi*f*t),
        // advanced each step by the fixed rotation rotC_/rotS_ (no per-step trig).
        std::vector<double> freqs_, rotC_, rotS_, pc_, ps_;
        std::vector<std::vector<double>> Xr_, Xi_;   // [port][freq]
        std::vector<double> Xsr_, Xsi_;              // [freq] source spectrum

        // CW steady-state lock-in. srcMode_: 0 = broadband pulse, 1 = CW sweep.
        // Accumulators are indexed [0..nports-1] = ports, [nports] = source.
        int srcMode_ = 0;
        double fCW_ = 0.0;
        double cwRampEnd_ = 0.0, cwBlockLen_ = 0.0, cwBlockStart_ = 0.0;
        double cwPc_ = 1.0, cwPs_ = 0.0, cwRotC_ = 1.0, cwRotS_ = 0.0;
        std::vector<double> cwI_, cwQ_;              // current window per port(+src)
        std::vector<double> cwPrevRe_, cwPrevIm_;    // previous window phasor
        std::vector<double> cwCurRe_, cwCurIm_;      // last finalized window phasor
        bool cwSteady_ = false;
        double cwResidual_ = 1e9;
        int cwBlocks_ = 0;

        std::vector<double> cx_, cy_, cz_;
    };

} // namespace waveguide
