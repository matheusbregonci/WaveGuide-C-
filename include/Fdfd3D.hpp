#pragma once

// Phase 3 of the frequency-domain solver: a 3D FDFD (finite-difference
// frequency-domain) core on a uniform Yee grid. It assembles the grad-div
// regularized curl-curl operator
//
//     L(E) = curl(curl E) - grad(div E) - k0^2 * eps_c * E
//
// (validated in Phase 2 to reproduce a rectangular-cavity resonance) into a
// sparse complex matrix and factors it with a direct solver. Per-cell complex
// permittivity + a PEC mask describe the microstrip (ground plane, dielectric
// substrate, copper trace). A driven source J gives the phasor field E at one
// frequency. PML and ports are added in a later step.
//
// Eigen is kept out of this header (pimpl) so the rest of the app need not see it.

#include <complex>
#include <vector>

namespace waveguide
{
    class Fdfd3D
    {
    public:
        using cd = std::complex<double>;

        Fdfd3D(int nx, int ny, int nz, double dx);
        ~Fdfd3D();
        Fdfd3D(const Fdfd3D&)            = delete;
        Fdfd3D& operator=(const Fdfd3D&) = delete;

        int    nx() const { return nx_; }
        int    ny() const { return ny_; }
        int    nz() const { return nz_; }
        double dx() const { return dx_; }
        int    cell(int i, int j, int k) const { return (k * ny_ + j) * nx_ + i; }

        // ---- material setup (defaults: air, no conductor) ----
        void setEps(int i, int j, int k, cd epsr) { eps_[cell(i,j,k)] = epsr; }
        void setPEC(int i, int j, int k)          { pec_[cell(i,j,k)] = 1;    }
        // PML (complex coordinate stretch) thickness in cells: nx on both x ends,
        // nz on both z ends, nyTop on the +y (air) side only (the -y face is the
        // ground plane / PEC). 0 disables a face.
        void setPML(int nx, int nyTop, int nz) { pmlX_=nx; pmlYtop_=nyTop; pmlZ_=nz; factored_=false; }
        // Optional grad-div stabilizer (default off: pure curl-curl - k0^2 eps).
        void setGradDiv(bool on) { gradDiv_=on; factored_=false; }
        // Linear solver: false = direct SparseLU (small grids), true = iterative
        // BiCGSTAB + incomplete-LU preconditioner (scales to large grids).
        void setIterative(bool on) { iterative_=on; factored_=false; }
        int    lastIters() const { return lastIters_; }   // iterative solve stats
        double lastError() const { return lastError_; }
        void reset();   // clear materials + source + factorization

        // ---- assemble + factor the operator at frequency f [Hz] ----
        bool factor(double fHz);
        bool factored() const { return factored_; }

        // ---- source (right-hand side) ----
        void clearSource();
        void addSource(int i, int j, int k, int comp, cd amp); // comp: 0=Ex,1=Ey,2=Ez

        // ---- solve L E = b for the current source ----
        bool solve();

        // ---- field access ----
        cd     E(int i, int j, int k, int comp) const;
        // Magnetic field H = (j / (w mu0)) * curl(E), from the solved E (plain
        // curl, valid in the physical region where S-parameters are measured).
        cd     H(int i, int j, int k, int comp) const;
        double omega() const { return w_; }   // angular frequency of the last factor()
        double normE() const;

        // Matrix-free operator apply (exposed for assembly / tests). k0sq = (w/c)^2.
        void applyL(double k0sq, const std::vector<cd>& e, std::vector<cd>& v) const;
        // Multiply the assembled sparse matrix by a vector (validation only).
        void matVec(const std::vector<cd>& e, std::vector<cd>& v) const;

    private:
        int    nx_, ny_, nz_, N_, M_;
        double dx_;
        std::vector<cd>            eps_;   // per-cell relative permittivity (complex)
        std::vector<unsigned char> pec_;   // per-cell PEC mask (all 3 comps forced 0)
        std::vector<cd>            b_, E_;  // source and solution, size M_ = 3*N_
        int  pmlX_ = 0, pmlYtop_ = 0, pmlZ_ = 0;   // PML thickness (cells)
        bool gradDiv_ = false;                       // optional grad-div stabilizer
        bool iterative_ = false;                      // solver: false=SparseLU, true=BiCGSTAB
        int  lastIters_ = 0;                          // iterative solve iteration count
        double lastError_ = 0.0;                      // iterative solve residual
        double w_ = 0.0;                             // angular frequency of last factor()
        bool factored_ = false;

        cd g(const std::vector<cd>& A, int i, int j, int k) const {
            if (i<0||i>=nx_||j<0||j>=ny_||k<0||k>=nz_) return cd(0,0);
            return A[cell(i,j,k)];
        }

        struct Impl;        // holds the Eigen sparse matrix + factorization
        Impl* impl_ = nullptr;
    };
}
