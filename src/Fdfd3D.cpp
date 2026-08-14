#include "Fdfd3D.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>

namespace waveguide
{
    namespace { constexpr double kC0 = 299792458.0; }

    using SpMat = Eigen::SparseMatrix<std::complex<double>>;

    struct Fdfd3D::Impl
    {
        SpMat A;
        Eigen::SparseLU<SpMat> lu;                                          // direct
        Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<std::complex<double>>> it; // iterative
    };

    Fdfd3D::Fdfd3D(int nx, int ny, int nz, double dx)
        : nx_(nx), ny_(ny), nz_(nz), N_(nx*ny*nz), M_(3*nx*ny*nz), dx_(dx)
    {
        eps_.assign(N_, cd(1.0, 0.0));
        pec_.assign(N_, 0);
        b_.assign(M_, cd(0,0));
        E_.assign(M_, cd(0,0));
        impl_ = new Impl();
    }
    Fdfd3D::~Fdfd3D() { delete impl_; }

    void Fdfd3D::reset()
    {
        std::fill(eps_.begin(), eps_.end(), cd(1,0));
        std::fill(pec_.begin(), pec_.end(), (unsigned char)0);
        std::fill(b_.begin(), b_.end(), cd(0,0));
        std::fill(E_.begin(), E_.end(), cd(0,0));
        factored_ = false;
    }

    // L(E) = curl(curl E) - grad(div E) - k0^2 eps_c E, with PEC cells forced to
    // an identity row (E = 0 inside the conductor). Matrix-free.
    void Fdfd3D::applyL(double k0sq, const std::vector<cd>& e, std::vector<cd>& v) const
    {
        const int N = N_; const double h = 1.0 / dx_;
        std::vector<cd> Ex(e.begin(), e.begin()+N), Ey(e.begin()+N, e.begin()+2*N), Ez(e.begin()+2*N, e.end());
        std::vector<cd> Hx(N), Hy(N), Hz(N), Dv(N);
        for (int k=0;k<nz_;++k) for (int j=0;j<ny_;++j) for (int i=0;i<nx_;++i) { int c=cell(i,j,k);
            Hx[c]=h*((g(Ez,i,j+1,k)-g(Ez,i,j,k))-(g(Ey,i,j,k+1)-g(Ey,i,j,k)));
            Hy[c]=h*((g(Ex,i,j,k+1)-g(Ex,i,j,k))-(g(Ez,i+1,j,k)-g(Ez,i,j,k)));
            Hz[c]=h*((g(Ey,i+1,j,k)-g(Ey,i,j,k))-(g(Ex,i,j+1,k)-g(Ex,i,j,k)));
            Dv[c]=0.5*h*((g(Ex,i+1,j,k)-g(Ex,i-1,j,k))+(g(Ey,i,j+1,k)-g(Ey,i,j-1,k))+(g(Ez,i,j,k+1)-g(Ez,i,j,k-1)));
        }
        for (int k=0;k<nz_;++k) for (int j=0;j<ny_;++j) for (int i=0;i<nx_;++i) { int c=cell(i,j,k);
            if (pec_[c]) { v[c]=e[c]; v[N+c]=e[N+c]; v[2*N+c]=e[2*N+c]; continue; }  // identity
            cd cx=h*((g(Hz,i,j,k)-g(Hz,i,j-1,k))-(g(Hy,i,j,k)-g(Hy,i,j,k-1)));
            cd cy=h*((g(Hx,i,j,k)-g(Hx,i,j,k-1))-(g(Hz,i,j,k)-g(Hz,i-1,j,k)));
            cd cz=h*((g(Hy,i,j,k)-g(Hy,i-1,j,k))-(g(Hx,i,j,k)-g(Hx,i,j-1,k)));
            cd gx=0.5*h*(g(Dv,i+1,j,k)-g(Dv,i-1,j,k));
            cd gy=0.5*h*(g(Dv,i,j+1,k)-g(Dv,i,j-1,k));
            cd gz=0.5*h*(g(Dv,i,j,k+1)-g(Dv,i,j,k-1));
            const cd m = k0sq * eps_[c];
            v[c]=cx-gx-m*Ex[c]; v[N+c]=cy-gy-m*Ey[c]; v[2*N+c]=cz-gz-m*Ez[c];
        }
    }

    bool Fdfd3D::factor(double fHz)
    {
        const double w = 2.0*3.14159265358979323846*fHz;
        w_ = w;
        const double k0sq = (w*w)/(kC0*kC0);
        const double h = 1.0 / dx_;
        // ---- PML: complex coordinate stretch 1/s, s = 1 + sigma/(j w eps0) ----
        // Each spatial derivative is multiplied by 1/s at its node (E-node for the
        // backward/central diffs, H-node = +1/2 cell for the forward diffs).
        constexpr double kEps0 = 8.8541878128e-12;
        const double eta0 = std::sqrt((4.0e-7*3.14159265358979323846) / kEps0);
        const double a = w * kEps0;
        // Complex-coordinate stretch s = kappa + sigma/(j w eps0), polynomially
        // graded over the PML depth (kappa helps absorb the near/evanescent field).
        auto invS = [&](double pos,int n,int nlo,int nhi)->cd {
            double rho = 0.0;
            if (nlo>0 && pos < nlo)                 rho = (nlo - pos) / nlo;
            else if (nhi>0 && pos > (n-1) - nhi)    rho = (pos - ((n-1) - nhi)) / nhi;
            rho = rho<0?0:(rho>1?1:rho);
            if (rho<=0) return cd(1,0);
            const double m=3.0, rm=std::pow(rho,m);
            const double sigma = 0.8*(m+1.0)/(eta0*dx_) * rm;
            const double kappa = 1.0 + 4.0*rm;                  // kappa: 1 -> 5
            return cd(1,0)/cd(kappa, -sigma/a);                 // 1/(kappa + sigma/(j a))
        };
        std::vector<cd> ixE(nx_),ixH(nx_),iyE(ny_),iyH(ny_),izE(nz_),izH(nz_);
        for (int i=0;i<nx_;++i){ ixE[i]=invS(i,nx_,pmlX_,pmlX_);   ixH[i]=invS(i+0.5,nx_,pmlX_,pmlX_); }
        for (int j=0;j<ny_;++j){ iyE[j]=invS(j,ny_,0,pmlYtop_);    iyH[j]=invS(j+0.5,ny_,0,pmlYtop_); } // PML on +y only
        for (int k=0;k<nz_;++k){ izE[k]=invS(k,nz_,pmlZ_,pmlZ_);   izH[k]=invS(k+0.5,nz_,pmlZ_,pmlZ_); }

        // Build the elementary operators as sparse matrices (each with its own
        // g()=0 boundary clipping AND its PML stretch), then compose
        // L = Cb*Cf - G*D - k0^2 eps. This reproduces the pure operator exactly
        // when the PML is off, and absorbs at the graded faces when it is on.
        auto inR = [&](int i,int j,int k){ return i>=0&&i<nx_&&j>=0&&j<ny_&&k>=0&&k<nz_; };
        std::vector<Eigen::Triplet<cd>> tCf, tCb, tD, tG;
        auto add = [&](std::vector<Eigen::Triplet<cd>>& T,int row,int comp,int i,int j,int k,cd v){
            if (inR(i,j,k)) T.emplace_back(row, comp*N_+cell(i,j,k), v); };
        auto addS = [&](std::vector<Eigen::Triplet<cd>>& T,int row,int i,int j,int k,cd v){ // scalar col
            if (inR(i,j,k)) T.emplace_back(row, cell(i,j,k), v); };
        for (int k=0;k<nz_;++k) for (int j=0;j<ny_;++j) for (int i=0;i<nx_;++i) {
            const int c=cell(i,j,k);
            const cd hxE=h*ixE[i], hyE=h*iyE[j], hzE=h*izE[k];   // backward/central: E-node
            const cd hxH=h*ixH[i], hyH=h*iyH[j], hzH=h*izH[k];   // forward: H-node
            const cd hhxE=0.5*hxE, hhyE=0.5*hyE, hhzE=0.5*hzE;   // central 1/(2dx) at E-node
            // forward curl E->H : d/dy,d/dz,... at the H-node
            add(tCf,0*N_+c,2,i,j+1,k,+hyH); add(tCf,0*N_+c,2,i,j,k,-hyH); add(tCf,0*N_+c,1,i,j,k+1,-hzH); add(tCf,0*N_+c,1,i,j,k,+hzH);
            add(tCf,1*N_+c,0,i,j,k+1,+hzH); add(tCf,1*N_+c,0,i,j,k,-hzH); add(tCf,1*N_+c,2,i+1,j,k,-hxH); add(tCf,1*N_+c,2,i,j,k,+hxH);
            add(tCf,2*N_+c,1,i+1,j,k,+hxH); add(tCf,2*N_+c,1,i,j,k,-hxH); add(tCf,2*N_+c,0,i,j+1,k,-hyH); add(tCf,2*N_+c,0,i,j,k,+hyH);
            // backward curl H->out : d/dy,d/dz,... at the E-node
            add(tCb,0*N_+c,2,i,j,k,+hyE); add(tCb,0*N_+c,2,i,j-1,k,-hyE); add(tCb,0*N_+c,1,i,j,k,-hzE); add(tCb,0*N_+c,1,i,j,k-1,+hzE);
            add(tCb,1*N_+c,0,i,j,k,+hzE); add(tCb,1*N_+c,0,i,j,k-1,-hzE); add(tCb,1*N_+c,2,i,j,k,-hxE); add(tCb,1*N_+c,2,i-1,j,k,+hxE);
            add(tCb,2*N_+c,1,i,j,k,+hxE); add(tCb,2*N_+c,1,i-1,j,k,-hxE); add(tCb,2*N_+c,0,i,j,k,-hyE); add(tCb,2*N_+c,0,i,j-1,k,+hyE);
            // divergence E->scalar (central, E-node)
            add(tD,c,0,i+1,j,k,+hhxE); add(tD,c,0,i-1,j,k,-hhxE);
            add(tD,c,1,i,j+1,k,+hhyE); add(tD,c,1,i,j-1,k,-hhyE);
            add(tD,c,2,i,j,k+1,+hhzE); add(tD,c,2,i,j,k-1,-hhzE);
            // gradient scalar->E (central, E-node)
            addS(tG,0*N_+c,i+1,j,k,+hhxE); addS(tG,0*N_+c,i-1,j,k,-hhxE);
            addS(tG,1*N_+c,i,j+1,k,+hhyE); addS(tG,1*N_+c,i,j-1,k,-hhyE);
            addS(tG,2*N_+c,i,j,k+1,+hhzE); addS(tG,2*N_+c,i,j,k-1,-hhzE);
        }
        SpMat Cf(M_,M_),Cb(M_,M_),Dm(N_,M_),Gm(M_,N_);
        Cf.setFromTriplets(tCf.begin(),tCf.end());
        Cb.setFromTriplets(tCb.begin(),tCb.end());
        Dm.setFromTriplets(tD.begin(),tD.end());
        Gm.setFromTriplets(tG.begin(),tG.end());
        // Pure curl-curl - k0^2 eps is already invertible for a driven problem (the
        // mass term lifts the gradient null space), and it avoids the checkerboard
        // that the wide central grad-div introduces. The grad-div is only kept as
        // an optional stabilizer for source-heavy / near-resonant cases.
        SpMat L = gradDiv_ ? SpMat((Cb*Cf) - (Gm*Dm)) : SpMat(Cb*Cf);

        // Add -k0^2 eps on the diagonal, then impose PEC rows as identity.
        std::vector<Eigen::Triplet<cd>> tA; tA.reserve(L.nonZeros() + M_);
        for (int col=0; col<L.outerSize(); ++col)
            for (SpMat::InnerIterator it(L,col); it; ++it) {
                const int row=int(it.row()); const int rc=row%N_;
                if (pec_[rc]) continue;   // PEC rows overwritten below
                tA.emplace_back(row, int(it.col()), it.value());
            }
        for (int comp=0;comp<3;++comp) for (int c=0;c<N_;++c) {
            const int row=comp*N_+c;
            if (pec_[c]) tA.emplace_back(row,row,cd(1,0));
            else         tA.emplace_back(row,row,-k0sq*eps_[c]);
        }
        impl_->A.resize(M_,M_);
        impl_->A.setFromTriplets(tA.begin(),tA.end());   // sums the mass into the diagonal
        impl_->A.makeCompressed();
        if (iterative_) {
            impl_->it.preconditioner().setDroptol(1e-3);
            impl_->it.preconditioner().setFillfactor(20);
            impl_->it.setMaxIterations(4000);
            impl_->it.setTolerance(1e-7);
            impl_->it.compute(impl_->A);
            factored_ = (impl_->it.info()==Eigen::Success);
        } else {
            impl_->lu.compute(impl_->A);
            factored_ = (impl_->lu.info()==Eigen::Success);
        }
        return factored_;
    }

    void Fdfd3D::clearSource() { std::fill(b_.begin(), b_.end(), cd(0,0)); }

    void Fdfd3D::addSource(int i, int j, int k, int comp, cd amp)
    {
        const int c = cell(i,j,k);
        if (pec_[c]) return;                 // no source inside a conductor
        b_[comp*N_ + c] += amp;
    }

    bool Fdfd3D::solve()
    {
        if (!factored_) return false;
        Eigen::Map<const Eigen::VectorXcd> bb(b_.data(), M_);
        Eigen::VectorXcd x;
        if (iterative_) {
            x = impl_->it.solve(bb);
            lastIters_ = int(impl_->it.iterations());
            lastError_ = impl_->it.error();
            if (impl_->it.info() != Eigen::Success && lastError_ > 1e-3) return false;
        } else {
            x = impl_->lu.solve(bb);
            if (impl_->lu.info() != Eigen::Success) return false;
        }
        for (int m=0; m<M_; ++m) E_[m] = x[m];
        return true;
    }

    void Fdfd3D::matVec(const std::vector<cd>& e, std::vector<cd>& v) const
    {
        Eigen::Map<const Eigen::VectorXcd> ev(e.data(), M_);
        Eigen::VectorXcd r = impl_->A * ev;
        v.resize(M_);
        for (int m=0;m<M_;++m) v[m]=r[m];
    }

    Fdfd3D::cd Fdfd3D::E(int i, int j, int k, int comp) const { return E_[comp*N_ + cell(i,j,k)]; }

    Fdfd3D::cd Fdfd3D::H(int i, int j, int k, int comp) const
    {
        // H = (j / (w mu0)) * (curl E), forward-difference Yee curl (plain, no PML).
        constexpr double kMu0 = 4.0e-7*3.14159265358979323846;
        const double h = 1.0/dx_;
        auto Eb = [&](int c, int a, int b, int d) -> cd {
            if (a<0||a>=nx_||b<0||b>=ny_||d<0||d>=nz_) return cd(0,0);
            return E_[c*N_ + cell(a,b,d)];
        };
        cd curl(0,0);
        if      (comp==0) curl = h*((Eb(2,i,j+1,k)-Eb(2,i,j,k)) - (Eb(1,i,j,k+1)-Eb(1,i,j,k)));
        else if (comp==1) curl = h*((Eb(0,i,j,k+1)-Eb(0,i,j,k)) - (Eb(2,i+1,j,k)-Eb(2,i,j,k)));
        else              curl = h*((Eb(1,i+1,j,k)-Eb(1,i,j,k)) - (Eb(0,i,j+1,k)-Eb(0,i,j,k)));
        return cd(0,1) / (w_*kMu0) * curl;
    }

    double Fdfd3D::normE() const { double s=0; for (auto& z : E_) s += std::norm(z); return std::sqrt(s); }
}
