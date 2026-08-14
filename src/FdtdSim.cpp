#include "FdtdSim.hpp"
#include "Colormap.hpp"

#include <algorithm>
#include <cmath>

namespace waveguide
{
    namespace
    {
        constexpr double kPi = 3.14159265358979323846;
        constexpr double kMu0 = 4.0 * kPi * 1e-7;
        constexpr double kEps0 = 8.8541878128e-12;
        const double kC0 = 1.0 / std::sqrt(kMu0 * kEps0);

        // Shared definition (Colormap.hpp).
        void fireColormap(float t, float& r, float& g, float& b)
        {
            fireColor(t, r, g, b);
        }
    } // namespace

    bool FdtdSim::solid(int i, int j, int k) const
    {
        if (i < 0 || i >= nx_ || j < 0 || j >= ny_ || k < 0 || k >= nz_) return false;
        return occ_[gidx(i, j, k)] != 0;
    }

    FdtdSim::FdtdSim(const VoxelMask& mask, double fminHz, double fmaxHz,
                     const std::vector<FdtdPort>& ports, double epsR, double muR)
    {
        epsR_ = (epsR > 1e-6) ? epsR : 1.0;
        muR_  = (muR  > 1e-6) ? muR  : 1.0;
        nx_ = mask.nx; ny_ = mask.ny; nz_ = mask.nz;
        W_ = mask.box.sizeX(); H_ = mask.box.sizeY(); D_ = mask.box.sizeZ();
        dx_ = W_ / std::max(1, nx_); dy_ = H_ / std::max(1, ny_); dz_ = D_ / std::max(1, nz_);
        dt_ = 0.98 / (kC0 * std::sqrt(1.0/(dx_*dx_) + 1.0/(dy_*dy_) + 1.0/(dz_*dz_)));
        occ_ = mask.occ;

        // Broadband Gaussian pulse covering [fmin, fmax]: centered at fc with a
        // time width tau chosen so the spectrum spans the band, so ONE run
        // excites every frequency and the port DFTs give S-params over the band.
        if (fmaxHz <= fminHz) { const double c = std::max(fminHz, kC0/(2.0*std::max(W_,H_))); fminHz = 0.5*c; fmaxHz = 1.5*c; }
        const double fc = 0.5 * (fminHz + fmaxHz);
        const double bw = std::max(1e6, fmaxHz - fminHz);
        f0_ = fc;
        tau_ = 1.0 / (kPi * bw * 0.5);   // Gaussian spectral half-width ~ bw/2
        t0_ = 4.0 * tau_;                // pulse center delay

        // DFT frequency grid over the sweep band.
        const int nf = 200;
        freqs_.resize(nf); rotC_.resize(nf); rotS_.resize(nf);
        pc_.assign(nf, 1.0); ps_.assign(nf, 0.0);
        Xsr_.assign(nf, 0.0); Xsi_.assign(nf, 0.0);
        for (int m = 0; m < nf; ++m) {
            freqs_[m] = fminHz + (fmaxHz - fminHz) * m / double(nf - 1);
            const double a = 2.0 * kPi * freqs_[m] * dt_;
            rotC_[m] = std::cos(a); rotS_[m] = std::sin(a);
        }

        Ex_.assign(nx_ * (ny_ + 1) * (nz_ + 1), 0.0);
        Ey_.assign((nx_ + 1) * ny_ * (nz_ + 1), 0.0);
        Ez_.assign((nx_ + 1) * (ny_ + 1) * nz_, 0.0);
        Hx_.assign((nx_ + 1) * ny_ * nz_, 0.0);
        Hy_.assign(nx_ * (ny_ + 1) * nz_, 0.0);
        Hz_.assign(nx_ * ny_ * (nz_ + 1), 0.0);

        // PEC/DOF flags: an E edge is a DOF only if all 4 cells touching it are
        // solid (else it lies on a PEC wall -> tangential E = 0).
        exFree_.assign(Ex_.size(), 0); eyFree_.assign(Ey_.size(), 0); ezFree_.assign(Ez_.size(), 0);
        for (int k = 0; k <= nz_; ++k) for (int j = 0; j <= ny_; ++j) for (int i = 0; i < nx_; ++i)
            exFree_[iEx(i,j,k)] = (solid(i,j-1,k-1)&&solid(i,j,k-1)&&solid(i,j-1,k)&&solid(i,j,k)) ? 1 : 0;
        for (int k = 0; k <= nz_; ++k) for (int j = 0; j < ny_; ++j) for (int i = 0; i <= nx_; ++i)
            eyFree_[iEy(i,j,k)] = (solid(i-1,j,k-1)&&solid(i,j,k-1)&&solid(i-1,j,k)&&solid(i,j,k)) ? 1 : 0;
        for (int k = 0; k < nz_; ++k) for (int j = 0; j <= ny_; ++j) for (int i = 0; i <= nx_; ++i)
            ezFree_[iEz(i,j,k)] = (solid(i-1,j-1,k)&&solid(i,j-1,k)&&solid(i-1,j,k)&&solid(i,j,k)) ? 1 : 0;

        // Keep only active port openings.
        for (const FdtdPort& p : ports) if (p.role != 0) ports_.push_back(p);
        mon_.resize(ports_.size());
        hist_.assign(ports_.size(), std::vector<float>(kHist, 0.0f));
        histPos_.assign(ports_.size(), 0);
        Xr_.assign(ports_.size(), std::vector<double>(freqs_.size(), 0.0));
        Xi_.assign(ports_.size(), std::vector<double>(freqs_.size(), 0.0));
        portSig_.assign(ports_.size(), 0.0);

        const int npml = std::max(4, std::min({nx_, ny_, nz_}) / 5);
        const double lossMax = 0.045;
        const int off = npml + 2;

        // Is cell (i,j,k) within a port's transverse rectangle?
        auto inRect = [&](const FdtdPort& p, int i, int j, int k) -> bool {
            double up, vp;
            if (p.axis == 2) { up=(i+0.5)*dx_; vp=(j+0.5)*dy_; }
            else if (p.axis == 0) { up=(j+0.5)*dy_; vp=(k+0.5)*dz_; }
            else { up=(i+0.5)*dx_; vp=(k+0.5)*dz_; }
            return up>=p.uMin && up<=p.uMax && vp>=p.vMin && vp<=p.vMax;
        };
        // Absorbing loss only near each opening (its transverse rect).
        auto lossAt = [&](int i, int j, int k) -> double {
            const int ci=std::clamp(i,0,nx_-1), cj=std::clamp(j,0,ny_-1), ck=std::clamp(k,0,nz_-1);
            double L = 0.0;
            for (const FdtdPort& p : ports_) {
                const int idxA=(p.axis==0)?ci:(p.axis==1)?cj:ck;
                const int dimA=(p.axis==0)?nx_:(p.axis==1)?ny_:nz_;
                const int dist=(p.side==0)? idxA : (dimA-1-idxA);
                if (dist>=npml || !inRect(p, ci, cj, ck)) continue;
                const double r=double(npml-dist)/npml; L=std::max(L, lossMax*r*r);
            }
            return L;
        };
        exLoss_.assign(Ex_.size(), 0.0); eyLoss_.assign(Ey_.size(), 0.0); ezLoss_.assign(Ez_.size(), 0.0);
        for (int k = 0; k <= nz_; ++k) for (int j = 0; j <= ny_; ++j) for (int i = 0; i < nx_; ++i)
            exLoss_[iEx(i,j,k)] = lossAt(i,j,k);
        for (int k = 0; k <= nz_; ++k) for (int j = 0; j < ny_; ++j) for (int i = 0; i <= nx_; ++i)
            eyLoss_[iEy(i,j,k)] = lossAt(i,j,k);
        for (int k = 0; k < nz_; ++k) for (int j = 0; j <= ny_; ++j) for (int i = 0; i <= nx_; ++i)
            ezLoss_[iEz(i,j,k)] = lossAt(i,j,k);

        // Monitor per opening (transverse E over its cells, modal half-sine
        // weight); input openings also drive the source.
        for (std::size_t pi = 0; pi < ports_.size(); ++pi)
        {
            const FdtdPort& p = ports_[pi];
            Monitor& mon = mon_[pi];
            const int axis=p.axis, side=p.side;
            const int dim=(axis==0)?nx_:(axis==1)?ny_:nz_;
            const int s=(side==0)? off : (dim-1-off);
            if (s<0 || s>=dim) continue;
            const double uspan=std::max(1e-12, p.uMax-p.uMin), vspan=std::max(1e-12, p.vMax-p.vMin);
            if (axis==2) { mon.comp=1; // Ey, half-sine across x
                for (int j=0;j<ny_;++j) for (int i=1;i<nx_;++i) {
                    if (!solid(i,j,s) || !inRect(p,i,j,s)) continue;
                    const int e=iEy(i,j,s); if(!eyFree_[e]) continue;
                    mon.idx.push_back(e); mon.w.push_back(float(std::sin(kPi*(((i+0.5)*dx_-p.uMin)/uspan))));
                }
            } else if (axis==0) { mon.comp=1; // Ey, half-sine across z
                for (int k=1;k<nz_;++k) for (int j=0;j<ny_;++j) {
                    if (!solid(s,j,k) || !inRect(p,s,j,k)) continue;
                    const int e=iEy(s,j,k); if(!eyFree_[e]) continue;
                    mon.idx.push_back(e); mon.w.push_back(float(std::sin(kPi*(((k+0.5)*dz_-p.vMin)/vspan))));
                }
            } else { mon.comp=0; // Ex, half-sine across z
                for (int k=1;k<nz_;++k) for (int i=0;i<nx_;++i) {
                    if (!solid(i,s,k) || !inRect(p,i,s,k)) continue;
                    const int e=iEx(i,s,k); if(!exFree_[e]) continue;
                    mon.idx.push_back(e); mon.w.push_back(float(std::sin(kPi*(((k+0.5)*dz_-p.vMin)/vspan))));
                }
            }
            if (p.role==1)
                for (std::size_t m=0;m<mon.idx.size();++m)
                { srcComp_.push_back(mon.comp); srcIdx_.push_back(mon.idx[m]); srcAmp_.push_back(mon.w[m]); }
        }

        cx_.assign(std::size_t(nx_) * ny_ * nz_, 0.0); cy_ = cx_; cz_ = cx_;
    }

    void FdtdSim::reset()
    {
        std::fill(Ex_.begin(), Ex_.end(), 0.0); std::fill(Ey_.begin(), Ey_.end(), 0.0); std::fill(Ez_.begin(), Ez_.end(), 0.0);
        std::fill(Hx_.begin(), Hx_.end(), 0.0); std::fill(Hy_.begin(), Hy_.end(), 0.0); std::fill(Hz_.begin(), Hz_.end(), 0.0);
        t_ = 0.0; nstep_ = 0; peak_ = 1e-30; emaxNow_ = 0.0; emaxEver_ = 1e-30;
        for (std::size_t p = 0; p < hist_.size(); ++p) { std::fill(hist_[p].begin(), hist_[p].end(), 0.0f); histPos_[p] = 0; }
        // Restart the running DFT (phasors back to t=0, accumulators cleared).
        std::fill(pc_.begin(), pc_.end(), 1.0); std::fill(ps_.begin(), ps_.end(), 0.0);
        std::fill(Xsr_.begin(), Xsr_.end(), 0.0); std::fill(Xsi_.begin(), Xsi_.end(), 0.0);
        for (auto& v : Xr_) std::fill(v.begin(), v.end(), 0.0);
        for (auto& v : Xi_) std::fill(v.begin(), v.end(), 0.0);
        cwSteady_ = false; cwResidual_ = 1e9; cwBlocks_ = 0;
        cwPc_ = 1.0; cwPs_ = 0.0;
        std::fill(cwI_.begin(), cwI_.end(), 0.0);
        std::fill(cwQ_.begin(), cwQ_.end(), 0.0);
    }

    void FdtdSim::portSpectrum(int p, int m, double& re, double& im) const
    {
        if (p < 0 || p >= int(Xr_.size()) || m < 0 || m >= int(freqs_.size())) { re = im = 0.0; return; }
        re = Xr_[p][m]; im = Xi_[p][m];
    }

    void FdtdSim::sourceSpectrum(int m, double& re, double& im) const
    {
        if (m < 0 || m >= int(freqs_.size())) { re = im = 0.0; return; }
        re = Xsr_[m]; im = Xsi_[m];
    }

    void FdtdSim::setCW(double fHz)
    {
        reset();                  // zero fields/time so each point is a fresh run
        srcMode_ = 1;
        fCW_ = (fHz > 0.0) ? fHz : 1e9;
        f0_ = fCW_;               // used by the ramp / reporting
        const double Tp = 1.0 / fCW_;
        cwRampEnd_   = 2.0 * Tp;   // smooth source turn-on before measuring
        cwBlockLen_  = 3.0 * Tp;   // lock-in window = 3 periods
        cwBlockStart_ = cwRampEnd_;
        const double a = 2.0 * kPi * fCW_ * dt_;
        cwRotC_ = std::cos(a); cwRotS_ = std::sin(a);
        cwPc_ = 1.0; cwPs_ = 0.0;
        const int np = int(ports_.size()) + 1;   // +1 for the source
        cwI_.assign(np, 0.0); cwQ_.assign(np, 0.0);
        cwPrevRe_.assign(np, 0.0); cwPrevIm_.assign(np, 0.0);
        cwCurRe_.assign(np, 0.0);  cwCurIm_.assign(np, 0.0);
        cwSteady_ = false; cwResidual_ = 1e9; cwBlocks_ = 0;
    }

    void FdtdSim::cwPortAmp(int p, double& re, double& im) const
    {
        if (p < 0 || p >= int(ports_.size())) { re = im = 0.0; return; }
        re = cwCurRe_[p]; im = cwCurIm_[p];
    }

    void FdtdSim::cwSourceAmp(double& re, double& im) const
    {
        const int s = int(cwCurRe_.size()) - 1;
        if (s < 0) { re = im = 0.0; return; }
        re = cwCurRe_[s]; im = cwCurIm_[s];
    }

    double FdtdSim::portRms(int i) const
    {
        if (i < 0 || i >= int(hist_.size()) || hist_[i].empty()) return 0.0;
        double s = 0.0; for (float v : hist_[i]) s += double(v) * v;
        return std::sqrt(s / hist_[i].size());
    }

    void FdtdSim::step(int n)
    {
        // Update coefficients scaled by the uniform filling medium.
        const double cE = dt_ / (kEps0 * epsR_);
        const double cH = dt_ / (kMu0 * muR_);
        auto gEx = [&](int i,int j,int k){ if(i<0||i>=nx_||j<0||j>ny_||k<0||k>nz_) return 0.0; return Ex_[iEx(i,j,k)]; };
        auto gEy = [&](int i,int j,int k){ if(i<0||i>nx_||j<0||j>=ny_||k<0||k>nz_) return 0.0; return Ey_[iEy(i,j,k)]; };
        auto gEz = [&](int i,int j,int k){ if(i<0||i>nx_||j<0||j>ny_||k<0||k>=nz_) return 0.0; return Ez_[iEz(i,j,k)]; };
        auto gHx = [&](int i,int j,int k){ if(i<0||i>nx_||j<0||j>=ny_||k<0||k>=nz_) return 0.0; return Hx_[iHx(i,j,k)]; };
        auto gHy = [&](int i,int j,int k){ if(i<0||i>=nx_||j<0||j>ny_||k<0||k>=nz_) return 0.0; return Hy_[iHy(i,j,k)]; };
        auto gHz = [&](int i,int j,int k){ if(i<0||i>=nx_||j<0||j>=ny_||k<0||k>nz_) return 0.0; return Hz_[iHz(i,j,k)]; };

        for (int s = 0; s < n; ++s)
        {
            // H update
            for (int k = 0; k < nz_; ++k) for (int j = 0; j < ny_; ++j) for (int i = 0; i <= nx_; ++i)
                Hx_[iHx(i,j,k)] -= cH * ((gEz(i,j+1,k)-gEz(i,j,k))/dy_ - (gEy(i,j,k+1)-gEy(i,j,k))/dz_);
            for (int k = 0; k < nz_; ++k) for (int j = 0; j <= ny_; ++j) for (int i = 0; i < nx_; ++i)
                Hy_[iHy(i,j,k)] -= cH * ((gEx(i,j,k+1)-gEx(i,j,k))/dz_ - (gEz(i+1,j,k)-gEz(i,j,k))/dx_);
            for (int k = 0; k <= nz_; ++k) for (int j = 0; j < ny_; ++j) for (int i = 0; i < nx_; ++i)
                Hz_[iHz(i,j,k)] -= cH * ((gEy(i+1,j,k)-gEy(i,j,k))/dx_ - (gEx(i,j+1,k)-gEx(i,j,k))/dy_);

            // E update (curl H)
            for (int k = 0; k <= nz_; ++k) for (int j = 0; j <= ny_; ++j) for (int i = 0; i < nx_; ++i)
                Ex_[iEx(i,j,k)] += cE * ((gHz(i,j,k)-gHz(i,j-1,k))/dy_ - (gHy(i,j,k)-gHy(i,j,k-1))/dz_);
            for (int k = 0; k <= nz_; ++k) for (int j = 0; j < ny_; ++j) for (int i = 0; i <= nx_; ++i)
                Ey_[iEy(i,j,k)] += cE * ((gHx(i,j,k)-gHx(i,j,k-1))/dz_ - (gHz(i,j,k)-gHz(i-1,j,k))/dx_);
            for (int k = 0; k < nz_; ++k) for (int j = 0; j <= ny_; ++j) for (int i = 0; i <= nx_; ++i)
                Ez_[iEz(i,j,k)] += cE * ((gHy(i,j,k)-gHy(i-1,j,k))/dx_ - (gHx(i,j,k)-gHx(i,j-1,k))/dy_);

            // soft source: either a broadband Gaussian pulse (srcMode_ 0, one
            // wavepacket exciting the whole band) or a continuous wave with a
            // smooth ramp (srcMode_ 1, VNA-style steady-state at one frequency).
            double src;
            if (srcMode_ == 1) {
                const double env = (t_ < cwRampEnd_) ? 0.5 * (1.0 - std::cos(kPi * t_ / cwRampEnd_)) : 1.0;
                src = env * std::sin(2.0 * kPi * fCW_ * t_);
            } else {
                const double tg = (t_ - t0_) / tau_;
                src = std::exp(-tg * tg) * std::sin(2.0 * kPi * f0_ * t_);
            }
            for (std::size_t s2 = 0; s2 < srcIdx_.size(); ++s2)
            {
                const double a = srcAmp_[s2] * src;
                if (srcComp_[s2] == 0) Ex_[srcIdx_[s2]] += a;
                else if (srcComp_[s2] == 1) Ey_[srcIdx_[s2]] += a;
                else Ez_[srcIdx_[s2]] += a;
            }

            // PEC: zero tangential E on walls; absorbing damping near ports.
            for (std::size_t e = 0; e < Ex_.size(); ++e) Ex_[e] = exFree_[e] ? Ex_[e] * (1.0 - exLoss_[e]) : 0.0;
            for (std::size_t e = 0; e < Ey_.size(); ++e) Ey_[e] = eyFree_[e] ? Ey_[e] * (1.0 - eyLoss_[e]) : 0.0;
            for (std::size_t e = 0; e < Ez_.size(); ++e) Ez_[e] = ezFree_[e] ? Ez_[e] * (1.0 - ezLoss_[e]) : 0.0;

            // Record the modal signal at each port (mode-weighted transverse E).
            const int nPorts = int(mon_.size());
            for (int pi = 0; pi < nPorts; ++pi) {
                double sig = 0.0;
                if (!mon_[pi].idx.empty()) {
                    const std::vector<double>& E = (mon_[pi].comp == 0) ? Ex_ : (mon_[pi].comp == 1) ? Ey_ : Ez_;
                    for (std::size_t m = 0; m < mon_[pi].idx.size(); ++m) sig += mon_[pi].w[m] * E[mon_[pi].idx[m]];
                }
                portSig_[pi] = sig;
                hist_[pi][histPos_[pi]] = float(sig);
                histPos_[pi] = (histPos_[pi] + 1) % kHist;
            }

            if (srcMode_ == 0) {
                // Broadband pulse: fold every port + the source into the running
                // multi-frequency DFT, then advance the DFT phasors one step.
                const int nf = int(freqs_.size());
                for (int pi = 0; pi < nPorts; ++pi) {
                    double* Xr = Xr_[pi].data(); double* Xi = Xi_[pi].data();
                    const double sdt = portSig_[pi] * dt_;
                    for (int m = 0; m < nf; ++m) { Xr[m] += sdt * pc_[m]; Xi[m] -= sdt * ps_[m]; }
                }
                const double srcdt = src * dt_;
                for (int m = 0; m < nf; ++m) {
                    Xsr_[m] += srcdt * pc_[m]; Xsi_[m] -= srcdt * ps_[m];
                    const double c = pc_[m], s = ps_[m];
                    pc_[m] = c * rotC_[m] - s * rotS_[m];
                    ps_[m] = s * rotC_[m] + c * rotS_[m];
                }
            } else {
                // CW steady state: single-frequency lock-in over windows of a few
                // periods; two consecutive windows agreeing => steady state.
                if (t_ >= cwRampEnd_) {
                    for (int pi = 0; pi < nPorts; ++pi) {
                        cwI_[pi] += portSig_[pi] * cwPc_ * dt_;
                        cwQ_[pi] += portSig_[pi] * cwPs_ * dt_;
                    }
                    cwI_[nPorts] += src * cwPc_ * dt_;
                    cwQ_[nPorts] += src * cwPs_ * dt_;

                    if (t_ - cwBlockStart_ >= cwBlockLen_) {
                        const int np = nPorts + 1;
                        double res = 0.0;
                        const double norm = 2.0 / cwBlockLen_;
                        for (int p = 0; p < np; ++p) {
                            const double re = norm * cwI_[p], im = -norm * cwQ_[p];
                            if (p < nPorts) { // convergence measured on ports only
                                const double d = std::hypot(re - cwPrevRe_[p], im - cwPrevIm_[p]);
                                const double mag = std::hypot(re, im) + 1e-12;
                                res = std::max(res, d / mag);
                            }
                            cwCurRe_[p] = re; cwCurIm_[p] = im;
                        }
                        if (cwBlocks_ >= 1 && res < 1e-2) cwSteady_ = true;
                        cwResidual_ = res;
                        cwPrevRe_ = cwCurRe_; cwPrevIm_ = cwCurIm_;
                        std::fill(cwI_.begin(), cwI_.end(), 0.0);
                        std::fill(cwQ_.begin(), cwQ_.end(), 0.0);
                        cwBlockStart_ += cwBlockLen_;
                        ++cwBlocks_;
                    }
                }
                const double c = cwPc_, s = cwPs_;
                cwPc_ = c * cwRotC_ - s * cwRotS_;
                cwPs_ = s * cwRotC_ + c * cwRotS_;
            }

            t_ += dt_; ++nstep_;
        }
    }

    void FdtdSim::syncField()
    {
        auto gEx = [&](int i,int j,int k){ if(i<0||i>=nx_||j<0||j>ny_||k<0||k>nz_) return 0.0; return Ex_[iEx(i,j,k)]; };
        auto gEy = [&](int i,int j,int k){ if(i<0||i>nx_||j<0||j>=ny_||k<0||k>nz_) return 0.0; return Ey_[iEy(i,j,k)]; };
        auto gEz = [&](int i,int j,int k){ if(i<0||i>nx_||j<0||j>ny_||k<0||k>=nz_) return 0.0; return Ez_[iEz(i,j,k)]; };
        double mx = 1e-30;
        for (int k = 0; k < nz_; ++k) for (int j = 0; j < ny_; ++j) for (int i = 0; i < nx_; ++i)
        {
            const std::size_t g = gidx(i,j,k);
            cx_[g] = 0.25*(gEx(i,j,k)+gEx(i,j+1,k)+gEx(i,j,k+1)+gEx(i,j+1,k+1));
            cy_[g] = 0.25*(gEy(i,j,k)+gEy(i+1,j,k)+gEy(i,j,k+1)+gEy(i+1,j,k+1));
            cz_[g] = 0.25*(gEz(i,j,k)+gEz(i+1,j,k)+gEz(i,j+1,k)+gEz(i+1,j+1,k));
            const double m = std::sqrt(cx_[g]*cx_[g]+cy_[g]*cy_[g]+cz_[g]*cz_[g]);
            if (m > mx) mx = m;
        }
        // running peak (slow decay) so the cloud normalization stays stable.
        peak_ = std::max(mx, peak_ * 0.995);
        if (peak_ < 1e-30) peak_ = 1e-30;
        // True (non-decaying) peak vs present amplitude, for sweep convergence.
        emaxNow_ = mx;
        emaxEver_ = std::max(emaxEver_, mx);
    }

    bool FdtdSim::inside(double x, double y, double z) const
    {
        const int i = int(x / std::max(1e-12, dx_)), j = int(y / std::max(1e-12, dy_)), k = int(z / std::max(1e-12, dz_));
        return solid(i, j, k);
    }

    double FdtdSim::sampleCell(const std::vector<double>& A, double ux, double uy, double uz) const
    {
        double fx = ux/dx_ - 0.5, fy = uy/dy_ - 0.5, fz = uz/dz_ - 0.5;
        fx = std::clamp(fx, 0.0, double(nx_-1)); fy = std::clamp(fy, 0.0, double(ny_-1)); fz = std::clamp(fz, 0.0, double(nz_-1));
        const int i0 = std::min(int(fx), std::max(0,nx_-2)), j0 = std::min(int(fy), std::max(0,ny_-2)), k0 = std::min(int(fz), std::max(0,nz_-2));
        const int i1 = std::min(i0+1,nx_-1), j1 = std::min(j0+1,ny_-1), k1 = std::min(k0+1,nz_-1);
        const double tx=fx-i0, ty=fy-j0, tz=fz-k0;
        auto V=[&](int i,int j,int k){ return A[gidx(i,j,k)]; };
        const double c00=V(i0,j0,k0)*(1-tx)+V(i1,j0,k0)*tx, c10=V(i0,j1,k0)*(1-tx)+V(i1,j1,k0)*tx;
        const double c01=V(i0,j0,k1)*(1-tx)+V(i1,j0,k1)*tx, c11=V(i0,j1,k1)*(1-tx)+V(i1,j1,k1)*tx;
        return (c00*(1-ty)+c10*ty)*(1-tz) + (c01*(1-ty)+c11*ty)*tz;
    }

    std::array<double,3> FdtdSim::fieldVector(double x, double y, double z, double) const
    {
        if (cx_.empty()) return {0,0,0};
        return {sampleCell(cx_,x,y,z), sampleCell(cy_,x,y,z), sampleCell(cz_,x,y,z)};
    }

    std::pair<double,double> FdtdSim::transverseField(double x, double y, double z, double) const
    {
        const auto v = fieldVector(x,y,z,0.0); return {v[0], v[1]};
    }

    std::vector<Particle> FdtdSim::sampleGrid(int nx, int ny, int nz, bool cutawayOn, float minIntensity, double) const
    {
        std::vector<Particle> out;
        if (nx<=1||ny<=1||nz<=1||cx_.empty()) return out;
        const double peak = peak_;
        const float halfW=float(W_)*0.5f, halfH=float(H_)*0.5f, halfD=float(D_)*0.5f;
        out.reserve(std::size_t(nx)*ny*nz);
        for (int i=0;i<nx;++i){ const double ux=W_*(i+0.5)/nx;
        for (int j=0;j<ny;++j){ const double uy=H_*(j+0.5)/ny;
        for (int k=0;k<nz;++k){ const double uz=D_*(k+0.5)/nz;
            if(!inside(ux,uy,uz)) continue;
            const auto g=fieldVector(ux,uy,uz,0.0);
            const double mag=std::sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
            const float t=float(std::min(1.0,mag/peak));
            if(t<minIntensity) continue;
            const float cx=float(ux)-halfW, cy=float(uy)-halfH, cz=float(uz)-halfD;
            if(cutawayOn && cx>0.0f && cy>0.0f && cz>0.0f) continue;
            float r,gg,bb; fireColormap(t,r,gg,bb);
            Particle p; p.x=cx;p.y=cy;p.z=cz;p.r=r;p.g=gg;p.b=bb;p.intensity=t;
            out.push_back(p);
        }}}
        return out;
    }

} // namespace waveguide
