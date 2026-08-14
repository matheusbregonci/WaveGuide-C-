#include "MicrostripSim.hpp"
#include "MicrostripXsec.hpp"
#include "Colormap.hpp"

#include <algorithm>
#include <cmath>

namespace waveguide
{
    namespace
    {
        constexpr double kPi  = 3.14159265358979323846;
        constexpr double kMu0 = 4.0 * kPi * 1e-7;
        constexpr double kEps0 = 8.8541878128e-12;
        const double kC0 = 1.0 / std::sqrt(kMu0 * kEps0);

        // Shared definition (Colormap.hpp).
        void fireColormap(float t, float& r, float& g, float& b)
        {
            fireColor(t, r, g, b);
        }

        // Painter-order material lookup: the last primitive containing the point
        // wins (so a via added after the substrate cuts through it).
        void sampleMat(const std::vector<MicrostripSim::Prim>& prims,
                       double x, double y, double z, int& mat, double& epsr, double& sigma)
        {
            mat = MicrostripSim::Air; epsr = 1.0; sigma = 0.0;
            for (const MicrostripSim::Prim& p : prims) {
                bool in = false;
                if (p.kind == 0)
                    in = x >= p.xmin && x <= p.xmax && y >= p.ymin && y <= p.ymax &&
                         z >= p.zmin && z <= p.zmax;
                else if (p.kind == 2) {   // disk/ring in x-z plane, extruded along y (trace patch)
                    const double dxr = x - p.cx, dzr = z - p.cz;
                    const double d2 = dxr*dxr + dzr*dzr;
                    in = y >= p.ymin && y <= p.ymax &&
                         d2 <= p.radius*p.radius && d2 >= p.rinner*p.rinner;
                }
                else {                    // kind 1: z-axis cylinder (a via)
                    const double dxr = x - p.cx, dyr = y - p.cy;
                    in = z >= p.zlo && z <= p.zhi && (dxr*dxr + dyr*dyr) <= p.radius*p.radius;
                }
                if (in) { mat = p.mat; epsr = p.epsr; sigma = p.sigma; }
            }
        }

        // Fill the 1D CPML stretch/recursion profiles for one axis of `n` nodes
        // with a PML of `npml` cells at each end. `offset` is 0 for E-node
        // (integer) positions, 0.5 for H-node (half-integer) positions. Grading
        // is polynomial (order m); sigma_max is the Gedney optimum. Outside the
        // PML: a=0 and invk=1, so the field update reduces to ordinary FDTD.
        void fillCpmlProfile(int n, int npmlLo, int npmlHi, double offset, double dx, double dt,
                             double eps0, double eta0,
                             std::vector<double>& b, std::vector<double>& a,
                             std::vector<double>& invk)
        {
            b.assign(n, 0.0); a.assign(n, 0.0); invk.assign(n, 1.0);
            if (npmlLo <= 0 && npmlHi <= 0) return;
            const double m = 3.0, kappaMax = 5.0, alphaMax = 0.05;
            const double sigmaMax = 0.8 * (m + 1.0) / (eta0 * dx);
            const double dLo = double(npmlLo), dHi = double(npmlHi);
            for (int i = 0; i < n; ++i) {
                const double pos = i + offset;      // node position in cells
                double rho = 0.0;                   // 0 interior .. 1 outer boundary
                if (npmlLo > 0 && pos <= dLo)                  rho = (dLo - pos) / dLo;
                else if (npmlHi > 0 && pos >= (n - 1) - dHi)   rho = (pos - ((n - 1) - dHi)) / dHi;
                rho = std::clamp(rho, 0.0, 1.0);
                if (rho <= 0.0) continue;           // interior: keep a=0, invk=1
                const double rm = std::pow(rho, m);
                const double sigma = sigmaMax * rm;
                const double kappa = 1.0 + (kappaMax - 1.0) * rm;
                const double alpha = alphaMax * (1.0 - rho);
                const double bb = std::exp(-(sigma / kappa + alpha) * dt / eps0);
                const double denom = kappa * (sigma + kappa * alpha);
                b[i]    = bb;
                a[i]    = (denom != 0.0) ? sigma * (bb - 1.0) / denom : 0.0;
                invk[i] = 1.0 / kappa;
            }
        }

        // Build monotone node coordinates on [0,L]: fine spacing hf inside the
        // padded band [f0,f1], coarsening to hc outside with a bounded neighbour
        // ratio r (>1). Guarantees node.front()=0 and node.back()=L. A uniform
        // grid is the special case hf==hc.
        std::vector<double> gradedAxis(double L, double f0, double f1,
                                       double hf, double hc, double r)
        {
            hf = std::max(hf, 1e-9); hc = std::max(hc, hf);
            f0 = std::clamp(f0, 0.0, L); f1 = std::clamp(f1, f0, L);
            std::vector<double> n; n.reserve(std::size_t(L / hf) + 8);
            n.push_back(0.0);
            double prev = hc;
            int guard = 0;
            while (n.back() < L - 1e-12 && guard++ < 2000000) {
                const double p = n.back();
                const double target = (p >= f0 && p <= f1) ? hf : hc;
                double h = std::clamp(target, prev / r, prev * r);   // limit neighbour ratio
                h = std::clamp(h, hf, hc);
                if (p + h > L - 0.5 * hf) h = L - p;                 // land exactly on L
                n.push_back(p + h);
                prev = h;
            }
            if (n.size() < 2) n = {0.0, L};
            n.back() = L;
            return n;
        }
    } // namespace

    MicrostripSim::MicrostripSim(double domX, double domY, double domZ, double dx,
                                 const std::vector<Prim>& prims,
                                 double srcX, double srcZ0, double srcZ1,
                                 double subY0, double subY1, double fcHz)
    {
        prims_ = prims;
        domX_ = domX; domY_ = domY; domZ_ = domZ; dx_ = dx;
        // Graded grid: X and Z stay uniform at the fine spacing dx; Y is graded so
        // the air above the substrate coarsens while the substrate/trace/fringing
        // stay fine. fineTopY = substrate top + a few substrate heights of fringe.
        const double fineTopY = std::min(domY, subY1 + std::max(3.0 * (subY1 - subY0), 8.0 * dx));
        buildGradedGrid(dx, fineTopY);
        nx_ = int(xe_.size()); ny_ = int(ye_.size()); nz_ = int(ze_.size());
        // CFL from the smallest cell on each axis (idxp_ max = 1/min-edge).
        const double ix = *std::max_element(idxp_.begin(), idxp_.end());
        const double iy = *std::max_element(idyp_.begin(), idyp_.end());
        const double iz = *std::max_element(idzp_.begin(), idzp_.end());
        dt_ = 0.99 / (kC0 * std::sqrt(ix*ix + iy*iy + iz*iz));
        fc_ = fcHz > 0 ? fcHz : kC0 / (10.0 * domX);
        // Continuous wave with a smooth turn-on ramp (a few periods), so the line
        // fills with a steady traveling wave and the field around the strip stays
        // visible (rather than a single pulse that passes and vanishes).
        t0_ = 3.0 / std::max(fc_, 1e6);
        spread_ = 1.0 / (2.0 * kPi * std::max(fc_, 1e6));

        const std::size_t N = std::size_t(nx_) * ny_ * nz_;
        Ex_.assign(N, 0.0); Ey_.assign(N, 0.0); Ez_.assign(N, 0.0);
        Hx_.assign(N, 0.0); Hy_.assign(N, 0.0); Hz_.assign(N, 0.0);
        epsx_.assign(N, 1.0); epsy_.assign(N, 1.0); epsz_.assign(N, 1.0);
        caEx_.assign(N, 1.0); cbEx_.assign(N, 0.0);
        caEy_.assign(N, 1.0); cbEy_.assign(N, 0.0);
        caEz_.assign(N, 1.0); cbEz_.assign(N, 0.0);
        pecx_.assign(N, 0); pecy_.assign(N, 0); pecz_.assign(N, 0);

        // Voxelize: sample the material at each E-component's half-cell offset and
        // bake the lossy update coefficients ca=(1-h)/(1+h), cb=(dt/eps)/(1+h)
        // with h = sigma*dt/(2*eps). PEC faces are lossless (node forced to 0).
        auto coef = [&](double er, double sg, double& ca, double& cb) {
            const double h = 0.5 * dt_ * sg / (kEps0 * er);
            ca = (1.0 - h) / (1.0 + h);
            cb = (dt_ / (kEps0 * er)) / (1.0 + h);
        };
        epsrSub_ = 1.0;
        for (int i = 0; i < nx_; ++i)
            for (int j = 0; j < ny_; ++j)
                for (int k = 0; k < nz_; ++k) {
                    const std::size_t id = idx(i, j, k);
                    int mat; double er, sg;
                    // Half-cell offsets follow the graded metric (X,Z uniform).
                    const double xh = xe_[i] + 0.5*dxp_[i];
                    const double yh = ye_[j] + 0.5*dyp_[j];
                    const double zh = ze_[k] + 0.5*dzp_[k];
                    sampleMat(prims, xh, ye_[j], ze_[k], mat, er, sg);
                    epsx_[id] = (mat == Pec) ? 1.0 : er; coef(epsx_[id], (mat==Pec)?0.0:sg, caEx_[id], cbEx_[id]); pecx_[id] = (mat == Pec);
                    sampleMat(prims, xe_[i], yh, ze_[k], mat, er, sg);
                    epsy_[id] = (mat == Pec) ? 1.0 : er; coef(epsy_[id], (mat==Pec)?0.0:sg, caEy_[id], cbEy_[id]); pecy_[id] = (mat == Pec);
                    sampleMat(prims, xe_[i], ye_[j], zh, mat, er, sg);
                    epsz_[id] = (mat == Pec) ? 1.0 : er; coef(epsz_[id], (mat==Pec)?0.0:sg, caEz_[id], cbEz_[id]); pecz_[id] = (mat == Pec);
                    if (mat == Dielectric) epsrSub_ = er;
                }

        // --- Conformal (Dey-Mittra) fill fractions on the trace layer ---------
        // Sub-sample the geometry to get, per cell, the unblocked length fraction of
        // the Ex/Ez edges and the unblocked area fraction of the Hy loop face. Only
        // the y-band spanned by the trace-layer copper is scanned; everywhere else
        // stays 1.0 (ordinary FDTD). This is what removes the ring's staircasing.
        fEx_.assign(N, 1.0f); fEz_.assign(N, 1.0f); aHy_.assign(N, 1.0f);
        {
            double tyLo = 1e30, tyHi = -1e30;
            for (const Prim& p : prims_)
                if (p.mat == Pec && p.ymin > 1e-9) { tyLo = std::min(tyLo, p.ymin); tyHi = std::max(tyHi, p.ymax); }
            if (tyHi > tyLo) {
                const int j0 = std::max(0,       coordToIndex(ye_, tyLo) - 1);
                const int j1 = std::min(ny_ - 1, coordToIndex(ye_, tyHi) + 1);
                const int NS = 8;                       // sub-samples per edge / per axis on the face
                auto isPec = [&](double x, double y, double z) {
                    int mat; double er, sg; sampleMat(prims_, x, y, z, mat, er, sg);
                    return mat == Pec;
                };
                for (int i = 0; i < nx_ - 1; ++i)
                    for (int j = j0; j <= j1; ++j)
                        for (int k = 0; k < nz_ - 1; ++k) {
                            const std::size_t id = idx(i, j, k);
                            const double x0 = xe_[i], x1 = xe_[i] + dxp_[i];
                            const double z0 = ze_[k], z1 = ze_[k] + dzp_[k];
                            const double yv = ye_[j];
                            int freeL = 0;
                            for (int s = 0; s < NS; ++s)               // Ex edge (along x)
                                if (!isPec(x0 + (x1-x0)*(s+0.5)/NS, yv, z0)) ++freeL;
                            fEx_[id] = float(freeL) / float(NS);
                            freeL = 0;
                            for (int s = 0; s < NS; ++s)               // Ez edge (along z)
                                if (!isPec(x0, yv, z0 + (z1-z0)*(s+0.5)/NS)) ++freeL;
                            fEz_[id] = float(freeL) / float(NS);
                            int freeA = 0;                             // Hy loop face (x-z)
                            for (int a = 0; a < NS; ++a)
                                for (int b = 0; b < NS; ++b)
                                    if (!isPec(x0 + (x1-x0)*(a+0.5)/NS, yv, z0 + (z1-z0)*(b+0.5)/NS)) ++freeA;
                            aHy_[id] = float(freeA) / float(NS*NS);
                        }
            }
        }

        // CPML absorbing boundaries on all 6 faces (replaces the Mur ABC).
        psiExy_.assign(N, 0.0); psiExz_.assign(N, 0.0);
        psiEyz_.assign(N, 0.0); psiEyx_.assign(N, 0.0);
        psiEzx_.assign(N, 0.0); psiEzy_.assign(N, 0.0);
        psiHxy_.assign(N, 0.0); psiHxz_.assign(N, 0.0);
        psiHyz_.assign(N, 0.0); psiHyx_.assign(N, 0.0);
        psiHzx_.assign(N, 0.0); psiHzy_.assign(N, 0.0);
        setupCpml();

        srcI_  = std::clamp(int(srcX  / dx_), 1, nx_ - 2);   // X uniform
        srcK0_ = std::clamp(int(srcZ0 / dx_), 1, nz_ - 2);   // Z uniform
        srcK1_ = std::clamp(int(srcZ1 / dx_), srcK0_, nz_ - 2);
        srcJ0_ = std::clamp(coordToIndex(ye_, subY0), 1, ny_ - 2);   // Y graded
        srcJ1_ = std::clamp(std::max(srcJ0_, coordToIndex(ye_, subY1)), srcJ0_, ny_ - 2);

        // Transmission-line probe geometry. The y voltage path (ground -> strip)
        // and the loop's y span are shared by both ports (same substrate); only
        // the z window (kc, kLoop0, kLoop1) varies per port, so each reference
        // plane hugs its own feed-line trace. Seeded from the source footprint;
        // layoutPorts() refines each port from the local trace below.
        jGnd_   = std::clamp(srcJ0_, 1, ny_ - 2);
        jStrip_ = std::clamp(srcJ1_, jGnd_ + 1, ny_ - 2);
        kc_     = std::clamp((srcK0_ + srcK1_) / 2, 1, nz_ - 2);
        jLoop0_ = std::clamp(srcJ1_ - 2, 1, ny_ - 2);
        jLoop1_ = std::clamp(srcJ1_ + 3, jLoop0_ + 1, ny_ - 2);
        kLoop0_ = std::clamp(srcK0_ - 2, 1, nz_ - 2);
        kLoop1_ = std::clamp(srcK1_ + 2, kLoop0_ + 1, nz_ - 2);
        kc1_ = kc_; kLoop0_1_ = kLoop0_; kLoop1_1_ = kLoop1_;

        // Feed-line extents: the lead-in trace touches x=0, the lead-out touches
        // x=domX (both are trace-layer PEC, i.e. ymin above the ground plane). The
        // filter body lives between feedInEndI_ and feedOutStartI_.
        feedInEndI_ = 0; feedOutStartI_ = nx_;
        {
            double inEndX = 0.0, outStartX = domX_;
            bool haveIn = false, haveOut = false;
            for (const Prim& p : prims_) {
                if (p.mat != Pec || p.kind != 0 || p.ymin <= 1e-9) continue;
                if (p.xmin <= 1.5 * dx_) { inEndX = std::max(inEndX, p.xmax); haveIn = true; }
                if (p.xmax >= domX_ - 1.5 * dx_) { outStartX = std::min(outStartX, p.xmin); haveOut = true; }
            }
            if (haveIn)  feedInEndI_    = std::clamp(int(inEndX / dx_),     1, nx_ - 2);
            if (haveOut) feedOutStartI_ = std::clamp(int(outStartX / dx_), 1, nx_ - 2);
        }

        // Port 1 (input reference) past the source near-field; port 2 (sense)
        // defaults to the output feed extremity. layoutPorts() snaps both onto
        // copper and rebuilds their probe boxes from the local trace.
        iP_     = std::clamp(nx_ - std::max(2, npx_ + 2) - 1, 1, nx_ - 2);
        iPort1_ = std::clamp(srcI_ + 2 * (jStrip_ - jGnd_) + 3, 1, nx_ - 2);
        layoutPorts();

        // Fixed reference impedance for the wave split = the feed line's quasi-
        // static (Hammerstad) Z0. Robust where a measured V/I would diverge (a
        // sense plane near a gap has I -> 0). Frequency-independent to first order,
        // so it persists across a sweep / setFrequency().
        {
            const double Wfeed = std::max(1e-9, srcZ1 - srcZ0);
            const double Hsub  = std::max(1e-9, subY1 - subY0);
            const XsecResult xr = hammerstadMicrostrip(Wfeed, Hsub, epsrSub_);
            zRef_ = (xr.ok && xr.Z0 > 0.0) ? xr.Z0 : 0.0;
        }
        setupLockin();

        cx_.assign(N, 0.0); cy_ = cx_; cz_ = cx_;
    }

    void MicrostripSim::buildGradedGrid(double hf, double fineTopY)
    {
        const double hc = 3.0 * hf;   // coarse air spacing (conservative)
        const double r  = 1.3;        // max neighbour ratio (accuracy/stability)
        // X, Z stay EXACTLY uniform at hf (a graded-generator sliver on those axes
        // would shrink the min cell and drop dt, cancelling the speedup). Only Y is
        // graded; its coarse top sliver is harmless since dt is set by the fine hf.
        auto uniformAxis = [&](double L) {
            const int n = std::max(4, int(L / hf));
            std::vector<double> v(std::size_t(n), 0.0);
            for (int i = 0; i < n; ++i) v[std::size_t(i)] = i * hf;
            return v;
        };
        xe_ = uniformAxis(domX_);
        ze_ = uniformAxis(domZ_);
        ye_ = gradedAxis(domY_, 0.0, fineTopY, hf, hc, r);

        // Primary edge lengths dp[i]=node[i+1]-node[i] (for stepH), dual edges
        // 0.5*(dp[i]+dp[i-1]) (for stepE); store inverses for the hot loops.
        auto metric = [&](const std::vector<double>& node, std::vector<double>& dp,
                          std::vector<double>& idp, std::vector<double>& idd) {
            const int n = int(node.size());
            dp.assign(std::size_t(n), hf); idp.assign(std::size_t(n), 1.0/hf); idd.assign(std::size_t(n), 1.0/hf);
            for (int i = 0; i < n - 1; ++i) dp[i] = node[i+1] - node[i];
            if (n >= 2) dp[n-1] = dp[n-2];
            for (int i = 0; i < n; ++i) idp[i] = 1.0 / std::max(dp[i], 1e-12);
            for (int i = 0; i < n; ++i) {
                const double dd = (i >= 1) ? 0.5*(dp[i] + dp[i-1]) : dp[0];
                idd[i] = 1.0 / std::max(dd, 1e-12);
            }
        };
        metric(xe_, dxp_, idxp_, idxd_);
        metric(ye_, dyp_, idyp_, idyd_);
        metric(ze_, dzp_, idzp_, idzd_);
    }

    int MicrostripSim::coordToIndex(const std::vector<double>& node, double p)
    {
        const int n = int(node.size());
        if (n == 0) return 0;
        if (p <= node.front()) return 0;
        if (p >= node.back())  return n - 1;
        int i = 0;                                  // node sizes are small; linear scan
        while (i + 1 < n && node[i+1] <= p) ++i;
        return i;
    }

    void MicrostripSim::setupLockin()
    {
        const double w = 2.0 * kPi * std::max(fc_, 1e6);
        rotc_ = std::cos(w * dt_);
        rots_ = std::sin(w * dt_);
        phc_ = 1.0; phs_ = 0.0;
    }

    void MicrostripSim::reset()
    {
        std::fill(Ex_.begin(), Ex_.end(), 0.0); std::fill(Ey_.begin(), Ey_.end(), 0.0); std::fill(Ez_.begin(), Ez_.end(), 0.0);
        std::fill(Hx_.begin(), Hx_.end(), 0.0); std::fill(Hy_.begin(), Hy_.end(), 0.0); std::fill(Hz_.begin(), Hz_.end(), 0.0);
        for (std::vector<double>* p : {&psiExy_, &psiExz_, &psiEyz_, &psiEyx_, &psiEzx_, &psiEzy_,
                                       &psiHxy_, &psiHxz_, &psiHyz_, &psiHyx_, &psiHzx_, &psiHzy_})
            std::fill(p->begin(), p->end(), 0.0);
        iPk_ = vPk_ = wmPk_ = wePk_ = 0.0;
        sparamActive_ = false; phc_ = 1.0; phs_ = 0.0;
        v1r_ = v1i_ = i1r_ = i1i_ = v2r_ = v2i_ = i2r_ = i2i_ = 0.0;
        t_ = 0.0; nstep_ = 0; peak_ = 1e-30;
    }

    void MicrostripSim::stepH()
    {
        const double cp = dt_ / kMu0;   // weight of both the curl and the CPML psi
        // Faraday differentiates E across PRIMARY edges -> per-axis 1/dxp,dyp,dzp.
        // Parallel over the outer plane: stepH only writes H and reads E, so
        // planes are independent (no data races).
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < nx_ - 1; ++i) {
            const double idxpI = idxp_[i];
            for (int j = 0; j < ny_ - 1; ++j) {
                const double idypJ = idyp_[j];
                for (int k = 0; k < nz_ - 1; ++k) {
                    const double idzpK = idzp_[k];
                    const std::size_t id = idx(i, j, k);
                    const std::size_t jx = idx(i, j + 1, k);
                    const std::size_t kz = idx(i, j, k + 1);
                    const std::size_t ix = idx(i + 1, j, k);
                    // Hx: d(Ez)/dy - d(Ey)/dz
                    const double dEzy = Ez_[jx] - Ez_[id];
                    const double dEyz = Ey_[kz] - Ey_[id];
                    psiHxy_[id] = bhY_[j] * psiHxy_[id] + ahY_[j] * (dEzy * idypJ);
                    psiHxz_[id] = bhZ_[k] * psiHxz_[id] + ahZ_[k] * (dEyz * idzpK);
                    Hx_[id] -= cp * (khY_[j] * dEzy * idypJ - khZ_[k] * dEyz * idzpK) + cp * (psiHxy_[id] - psiHxz_[id]);
                    // Hy: d(Ex)/dz - d(Ez)/dx. This Faraday loop lies in the x-z
                    // plane -- the one cut by the trace outline -- so the conformal
                    // (Dey-Mittra) weighting goes here.
                    {
                        double dExz = Ex_[kz] - Ex_[id];
                        double dEzx = Ez_[ix] - Ez_[id];
                        bool insidePec = false;
                        if (conformal_) {
                            const float a = aHy_[id];
                            if (a < 0.02f) insidePec = true;   // face essentially inside the PEC
                            else {
                                // STABLE conformal (Yu-Mittra / CFDTD-I): weight each E
                                // by the unblocked length of its own edge, but keep the
                                // FULL cell area. Strict Dey-Mittra also divides by the
                                // free area, which amplifies the update by 1/a (up to
                                // 50x here), blows past the CFL limit and diverges to
                                // NaN within a few steps -- so that term is dropped.
                                dExz = Ex_[kz]*double(fEx_[kz]) - Ex_[id]*double(fEx_[id]);
                                dEzx = Ez_[ix]*double(fEz_[ix]) - Ez_[id]*double(fEz_[id]);
                            }
                        }
                        if (insidePec) {
                            Hy_[id] = 0.0; psiHyz_[id] = 0.0; psiHyx_[id] = 0.0;
                        } else {
                            psiHyz_[id] = bhZ_[k] * psiHyz_[id] + ahZ_[k] * (dExz * idzpK);
                            psiHyx_[id] = bhX_[i] * psiHyx_[id] + ahX_[i] * (dEzx * idxpI);
                            Hy_[id] -= cp * (khZ_[k] * dExz * idzpK - khX_[i] * dEzx * idxpI)
                                     + cp * (psiHyz_[id] - psiHyx_[id]);
                        }
                    }
                    // Hz: d(Ey)/dx - d(Ex)/dy
                    const double dEyx = Ey_[ix] - Ey_[id];
                    const double dExy = Ex_[jx] - Ex_[id];
                    psiHzx_[id] = bhX_[i] * psiHzx_[id] + ahX_[i] * (dEyx * idxpI);
                    psiHzy_[id] = bhY_[j] * psiHzy_[id] + ahY_[j] * (dExy * idypJ);
                    Hz_[id] -= cp * (khX_[i] * dEyx * idxpI - khY_[j] * dExy * idypJ) + cp * (psiHzx_[id] - psiHzy_[id]);
                }
            }
        }
    }

    void MicrostripSim::stepE()
    {
        // Ampere differentiates H across DUAL edges -> per-axis 1/dxd,dyd,dzd.
        // Parallel over the outer plane: stepE only writes E and reads H.
        #pragma omp parallel for schedule(static)
        for (int i = 1; i < nx_; ++i) {
            const double idxdI = idxd_[i];
            for (int j = 1; j < ny_; ++j) {
                const double idydJ = idyd_[j];
                for (int k = 1; k < nz_; ++k) {
                    const double idzdK = idzd_[k];
                    const std::size_t id = idx(i, j, k);
                    const std::size_t im = idx(i - 1, j, k);
                    const std::size_t jm = idx(i, j - 1, k);
                    const std::size_t km = idx(i, j, k - 1);
                    // Ex: d(Hz)/dy - d(Hy)/dz
                    {
                        const double dHzy = Hz_[id] - Hz_[jm];
                        const double dHyz = Hy_[id] - Hy_[km];
                        psiExy_[id] = beY_[j] * psiExy_[id] + aeY_[j] * (dHzy * idydJ);
                        psiExz_[id] = beZ_[k] * psiExz_[id] + aeZ_[k] * (dHyz * idzdK);
                        const double curl = (keY_[j] * dHzy * idydJ - keZ_[k] * dHyz * idzdK)
                                          + (psiExy_[id] - psiExz_[id]);
                        Ex_[id] = caEx_[id] * Ex_[id] + cbEx_[id] * curl;
                    }
                    // Ey: d(Hx)/dz - d(Hz)/dx
                    {
                        const double dHxz = Hx_[id] - Hx_[km];
                        const double dHzx = Hz_[id] - Hz_[im];
                        psiEyz_[id] = beZ_[k] * psiEyz_[id] + aeZ_[k] * (dHxz * idzdK);
                        psiEyx_[id] = beX_[i] * psiEyx_[id] + aeX_[i] * (dHzx * idxdI);
                        const double curl = (keZ_[k] * dHxz * idzdK - keX_[i] * dHzx * idxdI)
                                          + (psiEyz_[id] - psiEyx_[id]);
                        Ey_[id] = caEy_[id] * Ey_[id] + cbEy_[id] * curl;
                    }
                    // Ez: d(Hy)/dx - d(Hx)/dy
                    {
                        const double dHyx = Hy_[id] - Hy_[im];
                        const double dHxy = Hx_[id] - Hx_[jm];
                        psiEzx_[id] = beX_[i] * psiEzx_[id] + aeX_[i] * (dHyx * idxdI);
                        psiEzy_[id] = beY_[j] * psiEzy_[id] + aeY_[j] * (dHxy * idydJ);
                        const double curl = (keX_[i] * dHyx * idxdI - keY_[j] * dHxy * idydJ)
                                          + (psiEzx_[id] - psiEzy_[id]);
                        Ez_[id] = caEz_[id] * Ez_[id] + cbEz_[id] * curl;
                    }
                    if (pecx_[id]) Ex_[id] = 0.0;
                    if (pecy_[id]) Ey_[id] = 0.0;
                    if (pecz_[id]) Ez_[id] = 0.0;
                }
            }
        }
    }

    void MicrostripSim::setupCpml()
    {
        // PML thickness per axis: aim for 8 cells but never let the two slabs
        // meet on a thin axis (a very short domain simply gets a thinner PML).
        auto pick = [](int n) { const int hi = (n - 2) / 2; return hi <= 0 ? 0 : std::min(8, hi); };
        npx_ = pick(nx_); npy_ = pick(ny_); npz_ = pick(nz_);
        const double eta0 = std::sqrt(kMu0 / kEps0);
        // E-node profiles (integer positions) drive the E update; H-node profiles
        // (half-cell offset) drive the H update. The y-axis gets PML only on the
        // TOP (air) side: the bottom is the PEC ground plane, so a bottom PML
        // would just absorb the guided substrate field it must not touch.
        // X, Z are uniform (dx_). The Y top-PML sits in the graded air, which is
        // uniform-coarse near the boundary -> use that local spacing (dyp_ at the
        // top) so sigma_max matches the actual cell size there.
        const double dyTop = (ny_ >= 2) ? dyp_[ny_ - 2] : dx_;
        fillCpmlProfile(nx_, npx_, npx_, 0.0, dx_,   dt_, kEps0, eta0, beX_, aeX_, keX_);
        fillCpmlProfile(nx_, npx_, npx_, 0.5, dx_,   dt_, kEps0, eta0, bhX_, ahX_, khX_);
        fillCpmlProfile(ny_, 0,    npy_, 0.0, dyTop, dt_, kEps0, eta0, beY_, aeY_, keY_);
        fillCpmlProfile(ny_, 0,    npy_, 0.5, dyTop, dt_, kEps0, eta0, bhY_, ahY_, khY_);
        fillCpmlProfile(nz_, npz_, npz_, 0.0, dx_,   dt_, kEps0, eta0, beZ_, aeZ_, keZ_);
        fillCpmlProfile(nz_, npz_, npz_, 0.5, dx_,   dt_, kEps0, eta0, bhZ_, ahZ_, khZ_);
    }

    void MicrostripSim::injectSource()
    {
        double val;
        if (pulseMode_) {
            // Broadband Gaussian pulse modulated at mid-band -> excites the whole
            // sweep band in one run (spectrum centred at pulseFc_, width ~1/tau).
            const double td = (t_ - pulseT0_) / pulseTau_;
            val = std::exp(-0.5 * td * td) * std::sin(2.0 * kPi * pulseFc_ * t_);
        } else {
            const double env = (t_ < t0_) ? 0.5 * (1.0 - std::cos(kPi * t_ / t0_)) : 1.0;
            val = env * std::sin(2.0 * kPi * fc_ * t_);
        }
        // Modal sheet: uniform vertical Ey (strip -> ground voltage) over the
        // substrate height and the strip width -> launches a clean quasi-TEM wave.
        for (int j = srcJ0_; j <= srcJ1_; ++j)
            for (int k = srcK0_; k <= srcK1_; ++k)
                Ey_[idx(srcI_, j, k)] += val;
    }

    void MicrostripSim::step(int n)
    {
        for (int s = 0; s < n; ++s) {
            stepH();
            stepE();          // CPML is folded into stepH/stepE (no separate ABC pass)
            injectSource();
            probe();          // update I, V and the transverse field energies
            t_ += dt_; ++nstep_;
        }
    }

    void MicrostripSim::syncField()
    {
        // Project the field the user asked to see (E or H) into the cell-centered
        // buffers the sampling/rendering paths read.
        const bool mag = (displayKind_ == FieldKind::Magnetic);
        const std::vector<double>& Sx = mag ? Hx_ : Ex_;
        const std::vector<double>& Sy = mag ? Hy_ : Ey_;
        const std::vector<double>& Sz = mag ? Hz_ : Ez_;
        double mx = 1e-30;
        for (int i = 0; i < nx_; ++i)
            for (int j = 0; j < ny_; ++j)
                for (int k = 0; k < nz_; ++k) {
                    const std::size_t id = idx(i, j, k);
                    cx_[id] = Sx[id]; cy_[id] = Sy[id]; cz_[id] = Sz[id];
                    const double m = std::sqrt(cx_[id]*cx_[id] + cy_[id]*cy_[id] + cz_[id]*cz_[id]);
                    if (m > mx) mx = m;
                }
        peak_ = std::max(mx, peak_ * 0.995);
        if (peak_ < 1e-30) peak_ = 1e-30;
    }

    void MicrostripSim::planeVI(int i, int kc, int kLoop0, int kLoop1, double& V, double& I) const
    {
        // Strip current: closed-loop integral of H around the trace in the
        // transverse (y-z) plane at x = i (Ampere's law -> enclosed axial
        // current, positive for a +x-directed current). Then the port voltage:
        // -integral of E.dl up through the substrate, ground -> strip. The z
        // window (kc, kLoop0, kLoop1) is per-port so each plane boxes its own
        // feed-line trace.
        // Each segment dl uses its local primary edge (z edges for the Hz sums,
        // y edges for the Hy sums; V integrates Ey along its y edges).
        double Iloop = 0.0;
        for (int k = kLoop0; k < kLoop1; ++k)
            Iloop += (Hz_[idx(i, jLoop1_, k)] - Hz_[idx(i, jLoop0_, k)]) * dzp_[k];
        for (int j = jLoop0_; j < jLoop1_; ++j)
            Iloop += (Hy_[idx(i, j, kLoop0)] - Hy_[idx(i, j, kLoop1)]) * dyp_[j];
        I = Iloop;
        double Vsum = 0.0;
        for (int j = jGnd_; j < jStrip_; ++j) Vsum -= Ey_[idx(i, j, kc)] * dyp_[j];
        V = Vsum;
    }

    bool MicrostripSim::traceZExtentAtX(double xm, int& kz0, int& kz1) const
    {
        // Scan the material scene for trace-layer copper (PEC above the ground
        // plane) crossing x = xm and return its z span; false if none (a gap).
        double zmin = 1e30, zmax = -1e30; bool found = false;
        for (const Prim& p : prims_) {
            if (p.mat != Pec || p.ymin <= 1e-9) continue;   // skip air / dielectric / ground
            if (p.kind == 0) {
                if (xm >= p.xmin && xm <= p.xmax) {
                    zmin = std::min(zmin, p.zmin); zmax = std::max(zmax, p.zmax); found = true;
                }
            } else if (p.kind == 2) {                        // disk / ring trace patch
                const double dxr = xm - p.cx;
                if (std::fabs(dxr) <= p.radius) {
                    const double hz = std::sqrt(std::max(0.0, p.radius * p.radius - dxr * dxr));
                    zmin = std::min(zmin, p.cz - hz); zmax = std::max(zmax, p.cz + hz); found = true;
                }
            }
        }
        if (!found) return false;
        kz0 = std::clamp(int(zmin / dx_),     1, nz_ - 2);
        kz1 = std::clamp(int(zmax / dx_) + 1, kz0 + 1, nz_ - 2);
        return true;
    }

    int MicrostripSim::snapToTrace(int i) const
    {
        const int lo = std::max(1, npx_ + 1), hi = nx_ - std::max(2, npx_ + 2);
        int kz0, kz1;
        i = std::clamp(i, lo, std::max(lo, hi));
        if (traceZExtentAtX((i + 0.5) * dx_, kz0, kz1)) return i;
        for (int d = 1; d < nx_; ++d) {                      // spiral out to nearest copper
            const int a = i + d, b = i - d;
            if (a <= hi && traceZExtentAtX((a + 0.5) * dx_, kz0, kz1)) return a;
            if (b >= lo && traceZExtentAtX((b + 0.5) * dx_, kz0, kz1)) return b;
        }
        return i;
    }

    void MicrostripSim::layoutPorts()
    {
        const int lo = std::max(1, npx_ + 1), hi = nx_ - std::max(2, npx_ + 2);
        auto boxFromTrace = [&](int iPlane, int& kc, int& kL0, int& kL1) {
            int kz0, kz1;
            if (traceZExtentAtX((iPlane + 0.5) * dx_, kz0, kz1)) {
                kc  = std::clamp((kz0 + kz1) / 2, 1, nz_ - 2);
                kL0 = std::clamp(kz0 - 2, 1, nz_ - 2);
                kL1 = std::clamp(kz1 + 2, kL0 + 1, nz_ - 2);
            }
        };

        // Port 1 (input reference): past the source near-field, kept on the input
        // feed (before the filter body starts).
        int i1 = srcI_ + 2 * (jStrip_ - jGnd_) + 3;
        if (feedInEndI_ > 0) i1 = std::min(i1, feedInEndI_ - 1);
        iPort1_ = snapToTrace(std::clamp(i1, lo, std::max(lo, hi)));
        boxFromTrace(iPort1_, kc1_, kLoop0_1_, kLoop1_1_);

        // Port 2 (sense / output reference): the requested plane, forced onto the
        // output feed extremity, never on the gap / filter body.
        int i2 = iP_;
        if (feedOutStartI_ < nx_) i2 = std::max(i2, feedOutStartI_ + 1);
        iP_ = snapToTrace(std::clamp(i2, lo, std::max(lo, hi)));
        boxFromTrace(iP_, kc_, kLoop0_, kLoop1_);
    }

    void MicrostripSim::probe()
    {
        double V, I;
        planeVI(iP_, kc_, kLoop0_, kLoop1_, V, I);

        // S-parameter lock-in: accumulate the V/I phasors at both reference
        // planes at the drive frequency (ratios are window-length independent,
        // so no 1/N normalisation is needed).
        if (sparamActive_) {
            double V1, I1, V2, I2;
            planeVI(iPort1_, kc1_, kLoop0_1_, kLoop1_1_, V1, I1);
            V2 = V; I2 = I;
            const double c = phc_, s = phs_;
            v1r_ += V1*c; v1i_ += V1*s; i1r_ += I1*c; i1i_ += I1*s;
            v2r_ += V2*c; v2i_ += V2*s; i2r_ += I2*c; i2i_ += I2*s;
            const double nc = c*rotc_ - s*rots_, ns = s*rotc_ + c*rots_;
            phc_ = nc; phs_ = ns;
        }

        // Broadband pulse sweep: accumulate the running DFT of the port V/I at
        // every requested frequency simultaneously (one pass over the band).
        if (pulseMode_) {
            double V1, I1;
            planeVI(iPort1_, kc1_, kLoop0_1_, kLoop1_1_, V1, I1);
            const double V2 = V, I2 = I;   // V,I already sampled at the sense plane
            const int M = int(pf_.size());
            for (int m = 0; m < M; ++m) {
                const double c = ppc_[m], s = pps_[m];
                pV1r_[m]+=V1*c; pV1i_[m]+=V1*s; pI1r_[m]+=I1*c; pI1i_[m]+=I1*s;
                pV2r_[m]+=V2*c; pV2i_[m]+=V2*s; pI2r_[m]+=I2*c; pI2i_[m]+=I2*s;
                const double nc = c*protC_[m] - s*protS_[m], ns = s*protC_[m] + c*protS_[m];
                ppc_[m] = nc; pps_[m] = ns;
            }
            // Slow-decaying peak of the port amplitude: falls once the pulse has
            // left through the CPML, so eNow_/eEver_ signals "rung down".
            const double amp = std::fabs(V1) + std::fabs(V2);
            eNow_  = std::max(amp, eNow_ * 0.999);
            eEver_ = std::max(eEver_, eNow_);
        }

        // Transverse energies per unit length at the probe plane (skip the PML
        // slabs). Wm' = 1/2 mu0 |H|^2, We' = 1/2 eps0 eps_r |E|^2, integrated
        // over the cross-section -> [J/m].
        double wm = 0.0, we = 0.0;
        for (int j = 1; j < ny_ - npy_; ++j)        // bottom = PEC ground (no PML)
            for (int k = npz_; k < nz_ - npz_; ++k) {
                const std::size_t id = idx(iP_, j, k);
                const double area = dyp_[j] * dzp_[k];   // per-cell cross-section area
                wm += (Hx_[id]*Hx_[id] + Hy_[id]*Hy_[id] + Hz_[id]*Hz_[id]) * area;
                we += (epsx_[id]*Ex_[id]*Ex_[id] + epsy_[id]*Ey_[id]*Ey_[id] + epsz_[id]*Ez_[id]*Ez_[id]) * area;
            }
        wm *= 0.5 * kMu0;
        we *= 0.5 * kEps0;

        // Peak-hold with slow decay: I and Wm' peak together (a quarter period
        // out of phase from V and We'), so tracking each maximum recovers the
        // CW amplitudes even under a standing wave.
        const double decay = 0.9999;
        iPk_  = std::max(std::fabs(I), iPk_  * decay);
        vPk_  = std::max(std::fabs(V), vPk_  * decay);
        wmPk_ = std::max(wm,           wmPk_ * decay);
        wePk_ = std::max(we,           wePk_ * decay);
    }

    double MicrostripSim::inductancePerLength() const
    {
        return (iPk_ > 1e-12) ? 2.0 * wmPk_ / (iPk_ * iPk_) : 0.0;
    }

    double MicrostripSim::capacitancePerLength() const
    {
        return (vPk_ > 1e-12) ? 2.0 * wePk_ / (vPk_ * vPk_) : 0.0;
    }

    double MicrostripSim::lineImpedance() const
    {
        const double L = inductancePerLength(), C = capacitancePerLength();
        return (L > 0.0 && C > 0.0) ? std::sqrt(L / C) : 0.0;
    }

    void MicrostripSim::setSourcePlaneX(double x)
    {
        // Keep the launch on the input feed, clear of the PML and never in the
        // filter body / gap; then rebuild both port probe boxes.
        const int lo = std::max(1, npx_ + 1), hi = nx_ - std::max(2, npx_ + 2);
        int i = std::clamp(int(x / dx_), lo, std::max(lo, hi));
        if (feedInEndI_ > 0) i = std::min(i, feedInEndI_ - 1);
        srcI_ = std::clamp(i, lo, std::max(lo, hi));
        layoutPorts();   // port 1 tracks the source; port 2 stays on the output feed
    }

    void MicrostripSim::setProbePlaneX(double x)
    {
        // Constrain the sense plane to the output feed extremity (off the gap /
        // filter body), then rebuild the port probe boxes.
        const int lo = std::max(1, npx_ + 1), hi = nx_ - std::max(2, npx_ + 2);
        int i = std::clamp(int(x / dx_), lo, std::max(lo, hi));
        if (feedOutStartI_ < nx_) i = std::max(i, feedOutStartI_ + 1);
        iP_ = std::clamp(i, lo, std::max(lo, hi));
        layoutPorts();
        iPk_ = vPk_ = wmPk_ = wePk_ = 0.0; // re-settle the readout at the new plane
    }

    std::array<double, 4> MicrostripSim::sourceMarker() const
    {
        // Sits on the copper (trace layer) at the source plane; radius = half the
        // strip width (z). srcJ1_ is the substrate top = trace layer.
        const double x = (srcI_ + 0.5) * dx_;
        const double y = ye_[std::size_t(std::clamp(srcJ1_, 0, ny_-1))];
        const double z = 0.5 * (srcK0_ + srcK1_) * dx_;
        const double r = std::max(0.5 * (srcK1_ - srcK0_) * dx_, dx_);
        return {x, y, z, r};
    }

    std::array<double, 4> MicrostripSim::senseMarker() const
    {
        // Sits on the copper at the sense plane; radius from the local trace width.
        const double x = (iP_ + 0.5) * dx_;
        const double y = ye_[std::size_t(std::clamp(jStrip_, 0, ny_-1))];
        int kz0, kz1; double z, r;
        if (traceZExtentAtX(x, kz0, kz1)) {
            z = 0.5 * (kz0 + kz1) * dx_;
            r = std::max(0.5 * (kz1 - kz0) * dx_, dx_);
        } else {                                   // no trace found -> fall back to source width
            z = 0.5 * (srcK0_ + srcK1_) * dx_;
            r = std::max(0.5 * (srcK1_ - srcK0_) * dx_, dx_);
        }
        return {x, y, z, r};
    }

    void MicrostripSim::setFrequency(double fcHz)
    {
        pulseMode_ = false;   // a CW retune leaves broadband-pulse mode
        if (fcHz > 0.0) fc_ = fcHz;
        t0_     = 3.0 / std::max(fc_, 1e6);
        spread_ = 1.0 / (2.0 * kPi * std::max(fc_, 1e6));
        reset();        // clear the field, psi and probe peaks for a clean relaunch
        setupLockin();  // retune the S-parameter lock-in to the new frequency
    }

    void MicrostripSim::sparamStartWindow()
    {
        v1r_ = v1i_ = i1r_ = i1i_ = v2r_ = v2i_ = i2r_ = i2i_ = 0.0;
        phc_ = 1.0; phs_ = 0.0;
        sparamActive_ = true;
    }

    double MicrostripSim::sparamZ0() const
    {
        // Self-consistent characteristic impedance from |V2|/|I2| at port 2. Port 2
        // now sits on the uniform OUTPUT feed (layoutPorts keeps it off the gap /
        // filter body), so this is the true Z0 of the discrete FDTD line -> the
        // wave split stays energy-consistent and passive |S| <= 1. Using a FIXED
        // analytic Z0 here instead mismatches the discrete line and manufactures
        // |S| > 1 (dips above 0 dB); the fixed feed Z0 (zRef_) is only a fallback
        // for the degenerate case where the port current collapses.
        const double v2 = v2r_*v2r_ + v2i_*v2i_, i2 = i2r_*i2r_ + i2i_*i2i_;
        if (i2 > 0.0) return std::sqrt(v2 / i2);
        return (zRef_ > 0.0) ? zRef_ : 50.0;
    }

    void MicrostripSim::s11(double& re, double& im) const
    {
        // Forward/backward split at port 1: V1+- = (V1 +- Z0*I1)/2. S11 = V1-/V1+.
        const double Z0 = sparamZ0();
        const double fpr = 0.5*(v1r_ + Z0*i1r_), fpi = 0.5*(v1i_ + Z0*i1i_); // V1+
        const double bmr = 0.5*(v1r_ - Z0*i1r_), bmi = 0.5*(v1i_ - Z0*i1i_); // V1-
        const double d = fpr*fpr + fpi*fpi + 1e-300;
        re = (bmr*fpr + bmi*fpi) / d;
        im = (bmi*fpr - bmr*fpi) / d;
    }

    void MicrostripSim::s21(double& re, double& im) const
    {
        // Transmission: forward wave at port 2 over forward (incident) at port 1.
        const double Z0 = sparamZ0();
        const double fpr = 0.5*(v1r_ + Z0*i1r_), fpi = 0.5*(v1i_ + Z0*i1i_); // V1+
        const double tpr = 0.5*(v2r_ + Z0*i2r_), tpi = 0.5*(v2i_ + Z0*i2i_); // V2+
        const double d = fpr*fpr + fpi*fpi + 1e-300;
        re = (tpr*fpr + tpi*fpi) / d;
        im = (tpi*fpr - tpr*fpi) / d;
    }

    void MicrostripSim::startPulseSweep(const std::vector<double>& freqsHz)
    {
        reset();                       // clean field, phasors and t_ = 0
        pf_ = freqsHz;
        const int M = int(pf_.size());
        double fmin = 1e30, fmax = 0.0;
        for (double f : pf_) { fmin = std::min(fmin, f); fmax = std::max(fmax, f); }
        if (M == 0 || fmax <= 0.0) { pulseMode_ = false; return; }
        // Modulated-Gaussian pulse centred at mid-band, spread wide enough that
        // its spectrum covers [fmin, fmax] with usable amplitude at the edges.
        pulseFc_  = 0.5 * (fmin + fmax);
        const double halfBand = std::max(0.5 * (fmax - fmin), 0.25 * pulseFc_);
        pulseTau_ = 1.0 / (kPi * halfBand);
        pulseT0_  = 4.0 * pulseTau_;   // launch delay: pulse ~0 at t = 0
        protC_.assign(M, 0.0); protS_.assign(M, 0.0);
        ppc_.assign(M, 1.0);   pps_.assign(M, 0.0);
        for (int m = 0; m < M; ++m) {
            const double w = 2.0 * kPi * pf_[m];
            protC_[m] = std::cos(w * dt_);
            protS_[m] = std::sin(w * dt_);
        }
        pV1r_.assign(M, 0.0); pV1i_.assign(M, 0.0); pI1r_.assign(M, 0.0); pI1i_.assign(M, 0.0);
        pV2r_.assign(M, 0.0); pV2i_.assign(M, 0.0); pI2r_.assign(M, 0.0); pI2i_.assign(M, 0.0);
        eNow_ = 0.0; eEver_ = 1e-30;
        pulseMode_ = true;
    }

    double MicrostripSim::pulseZ0(int m) const
    {
        if (m < 0 || m >= int(pf_.size())) return 0.0;
        const double v2 = std::hypot(pV2r_[m], pV2i_[m]);
        const double i2 = std::hypot(pI2r_[m], pI2i_[m]);
        return (i2 > 1e-300) ? v2 / i2 : 0.0;
    }

    void MicrostripSim::pulseS(int m, double& s11r, double& s11i,
                               double& s21r, double& s21i, double z0ext) const
    {
        s11r = s11i = s21r = s21i = 0.0;
        if (m < 0 || m >= int(pf_.size())) return;
        // Reference Z0 for the forward/backward TEM split. Prefer z0ext (the THRU
        // line's feed Z0): a reflective filter has a tiny, noisy port-2 current, so
        // its own |V2|/|I2| is garbage and blows |S| far past 1. The THRU always
        // carries a clean forward wave, so its |V2|/|I2| is the true discrete feed
        // Z0. Fall back to the local ratio, then zRef_, only if no override.
        const double v2 = std::hypot(pV2r_[m], pV2i_[m]);
        const double i2 = std::hypot(pI2r_[m], pI2i_[m]);
        const double Z0 = (z0ext > 0.0) ? z0ext
                        : (i2 > 1e-300) ? v2 / i2 : (zRef_ > 0.0 ? zRef_ : 50.0);
        const double fpr = 0.5*(pV1r_[m] + Z0*pI1r_[m]), fpi = 0.5*(pV1i_[m] + Z0*pI1i_[m]); // V1+
        const double bmr = 0.5*(pV1r_[m] - Z0*pI1r_[m]), bmi = 0.5*(pV1i_[m] - Z0*pI1i_[m]); // V1-
        const double tpr = 0.5*(pV2r_[m] + Z0*pI2r_[m]), tpi = 0.5*(pV2i_[m] + Z0*pI2i_[m]); // V2+
        const double d = fpr*fpr + fpi*fpi + 1e-300;
        s11r = (bmr*fpr + bmi*fpi) / d;  s11i = (bmi*fpr - bmr*fpi) / d;
        s21r = (tpr*fpr + tpi*fpi) / d;  s21i = (tpi*fpr - tpr*fpi) / d;
    }

    bool MicrostripSim::inside(double x, double y, double z) const
    {
        return x >= 0 && x <= domX_ && y >= 0 && y <= domY_ && z >= 0 && z <= domZ_;
    }

    double MicrostripSim::sampleCell(const std::vector<double>& A, double x, double y, double z) const
    {
        double fx = x/dx_, fz = z/dx_;              // X, Z uniform
        double fy;                                  // Y graded: fractional index
        {
            const int j = coordToIndex(ye_, y);
            const int j1 = std::min(j + 1, ny_ - 1);
            fy = (ye_[j1] > ye_[j]) ? j + (y - ye_[j]) / (ye_[j1] - ye_[j]) : double(j);
        }
        fx = std::clamp(fx, 0.0, double(nx_-1)); fy = std::clamp(fy, 0.0, double(ny_-1)); fz = std::clamp(fz, 0.0, double(nz_-1));
        const int i0 = std::min(int(fx), nx_-2), j0 = std::min(int(fy), ny_-2), k0 = std::min(int(fz), nz_-2);
        const double tx = fx-i0, ty = fy-j0, tz = fz-k0;
        auto V = [&](int i,int j,int k){ return A[idx(i,j,k)]; };
        const double c00 = V(i0,j0,k0)*(1-tx)+V(i0+1,j0,k0)*tx, c10 = V(i0,j0+1,k0)*(1-tx)+V(i0+1,j0+1,k0)*tx;
        const double c01 = V(i0,j0,k0+1)*(1-tx)+V(i0+1,j0,k0+1)*tx, c11 = V(i0,j0+1,k0+1)*(1-tx)+V(i0+1,j0+1,k0+1)*tx;
        return (c00*(1-ty)+c10*ty)*(1-tz) + (c01*(1-ty)+c11*ty)*tz;
    }

    std::array<double,3> MicrostripSim::fieldVector(double x, double y, double z, double) const
    {
        if (cx_.empty()) return {0,0,0};
        return {sampleCell(cx_,x,y,z), sampleCell(cy_,x,y,z), sampleCell(cz_,x,y,z)};
    }

    std::pair<double,double> MicrostripSim::transverseField(double x, double y, double z, double) const
    {
        const auto v = fieldVector(x,y,z,0.0); return {v[0], v[1]};
    }

    double MicrostripSim::surfaceCurrent(double x, double y, double z) const
    {
        // |Js| = |n x H| = |H tangential| to a horizontal conductor face = the
        // in-plane magnetic field (Hx, Hz). Sampled from the live H buffers, so
        // it is valid even while the display shows E.
        if (Hx_.empty()) return 0.0;
        const double hx = sampleCell(Hx_, x, y, z);
        const double hz = sampleCell(Hz_, x, y, z);
        return std::sqrt(hx*hx + hz*hz);
    }

    std::vector<Particle> MicrostripSim::sampleGrid(int nx, int ny, int nz, bool cutawayOn, float minIntensity, double) const
    {
        std::vector<Particle> out;
        if (nx<=1||ny<=1||nz<=1||cx_.empty()) return out;
        const double peak = peak_;
        const float hW=float(domX_)*0.5f, hH=float(domY_)*0.5f, hD=float(domZ_)*0.5f;
        out.reserve(std::size_t(nx)*ny*nz);
        for (int i=0;i<nx;++i){ const double ux=domX_*(i+0.5)/nx;
        for (int j=0;j<ny;++j){ const double uy=domY_*(j+0.5)/ny;
        for (int k=0;k<nz;++k){ const double uz=domZ_*(k+0.5)/nz;
            const auto g=fieldVector(ux,uy,uz,0.0);
            const double mag=std::sqrt(g[0]*g[0]+g[1]*g[1]+g[2]*g[2]);
            const float t=float(std::min(1.0,mag/peak));
            if(t<minIntensity) continue;
            const float cx=float(ux)-hW, cy=float(uy)-hH, cz=float(uz)-hD;
            if(cutawayOn && cx>0.0f && cy>0.0f && cz>0.0f) continue;
            float r,gg,bb; fireColormap(t,r,gg,bb);
            Particle p; p.x=cx;p.y=cy;p.z=cz;p.r=r;p.g=gg;p.b=bb;p.intensity=t;
            out.push_back(p);
        }}}
        return out;
    }

} // namespace waveguide
