#include "Geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace waveguide
{

    // CSG tree node. Leaves are primitives; internal nodes are boolean ops.
    struct Geometry::Node
    {
        enum class Kind { Box, CylinderZ, Union, Diff, Intersect } kind;

        // Box: center + half sizes. CylinderZ: axis (px,py), z range, radius.
        double cx = 0, cy = 0, cz = 0, hx = 0, hy = 0, hz = 0;
        double px = 0, py = 0, z0 = 0, z1 = 0, r = 0;

        std::shared_ptr<const Node> a, b; // children (ops)
    };

    namespace
    {
        bool nodeInside(const Geometry::Node& n, double x, double y, double z)
        {
            using K = Geometry::Node::Kind;
            switch (n.kind)
            {
            case K::Box:
                return std::abs(x - n.cx) <= n.hx &&
                       std::abs(y - n.cy) <= n.hy &&
                       std::abs(z - n.cz) <= n.hz;
            case K::CylinderZ:
            {
                const double dx = x - n.px, dy = y - n.py;
                return dx * dx + dy * dy <= n.r * n.r && z >= n.z0 && z <= n.z1;
            }
            case K::Union:
                return nodeInside(*n.a, x, y, z) || nodeInside(*n.b, x, y, z);
            case K::Diff:
                return nodeInside(*n.a, x, y, z) && !nodeInside(*n.b, x, y, z);
            case K::Intersect:
                return nodeInside(*n.a, x, y, z) && nodeInside(*n.b, x, y, z);
            }
            return false;
        }

        Aabb merge(const Aabb& p, const Aabb& q)
        {
            return {std::min(p.xmin, q.xmin), std::min(p.ymin, q.ymin),
                    std::min(p.zmin, q.zmin), std::max(p.xmax, q.xmax),
                    std::max(p.ymax, q.ymax), std::max(p.zmax, q.zmax)};
        }
        Aabb clip(const Aabb& p, const Aabb& q)
        {
            return {std::max(p.xmin, q.xmin), std::max(p.ymin, q.ymin),
                    std::max(p.zmin, q.zmin), std::min(p.xmax, q.xmax),
                    std::min(p.ymax, q.ymax), std::min(p.zmax, q.zmax)};
        }

        Aabb nodeBounds(const Geometry::Node& n)
        {
            using K = Geometry::Node::Kind;
            switch (n.kind)
            {
            case K::Box:
                return {n.cx - n.hx, n.cy - n.hy, n.cz - n.hz,
                        n.cx + n.hx, n.cy + n.hy, n.cz + n.hz};
            case K::CylinderZ:
                return {n.px - n.r, n.py - n.r, n.z0,
                        n.px + n.r, n.py + n.r, n.z1};
            case K::Union:
                return merge(nodeBounds(*n.a), nodeBounds(*n.b));
            case K::Diff:
                return nodeBounds(*n.a); // subtracting can only shrink `a`
            case K::Intersect:
                return clip(nodeBounds(*n.a), nodeBounds(*n.b));
            }
            return {0, 0, 0, 0, 0, 0};
        }

        std::shared_ptr<Geometry::Node> makeOp(Geometry::Node::Kind k,
                                               std::shared_ptr<const Geometry::Node> a,
                                               std::shared_ptr<const Geometry::Node> b)
        {
            auto n = std::make_shared<Geometry::Node>();
            n->kind = k;
            n->a = std::move(a);
            n->b = std::move(b);
            return n;
        }
    } // namespace

    std::size_t VoxelMask::solidCount() const
    {
        std::size_t c = 0;
        for (std::uint8_t v : occ) c += v ? 1u : 0u;
        return c;
    }

    Geometry Geometry::box(double cx, double cy, double cz,
                           double sx, double sy, double sz)
    {
        auto n = std::make_shared<Node>();
        n->kind = Node::Kind::Box;
        n->cx = cx; n->cy = cy; n->cz = cz;
        n->hx = 0.5 * sx; n->hy = 0.5 * sy; n->hz = 0.5 * sz;
        return Geometry(n);
    }

    Geometry Geometry::cylinderZ(double cx, double cy,
                                 double z0, double z1, double radius)
    {
        auto n = std::make_shared<Node>();
        n->kind = Node::Kind::CylinderZ;
        n->px = cx; n->py = cy;
        n->z0 = std::min(z0, z1); n->z1 = std::max(z0, z1);
        n->r = radius;
        return Geometry(n);
    }

    Geometry Geometry::unite(const Geometry& other) const
    {
        if (empty()) return other;
        if (other.empty()) return *this;
        return Geometry(makeOp(Node::Kind::Union, root_, other.root_));
    }
    Geometry Geometry::subtract(const Geometry& other) const
    {
        if (empty() || other.empty()) return *this;
        return Geometry(makeOp(Node::Kind::Diff, root_, other.root_));
    }
    Geometry Geometry::intersect(const Geometry& other) const
    {
        if (empty() || other.empty()) return Geometry();
        return Geometry(makeOp(Node::Kind::Intersect, root_, other.root_));
    }

    bool Geometry::inside(double x, double y, double z) const
    {
        return root_ && nodeInside(*root_, x, y, z);
    }

    Aabb Geometry::bounds() const
    {
        if (!root_) return {0, 0, 0, 0, 0, 0};
        return nodeBounds(*root_);
    }

    VoxelMask Geometry::voxelize(int nx, int ny, int nz) const
    {
        return voxelize(nx, ny, nz, bounds());
    }

    VoxelMask Geometry::voxelize(int nx, int ny, int nz, const Aabb& region) const
    {
        VoxelMask m;
        m.nx = std::max(1, nx); m.ny = std::max(1, ny); m.nz = std::max(1, nz);
        m.box = region;
        m.occ.assign(std::size_t(m.nx) * m.ny * m.nz, 0);
        if (!root_ || !region.valid()) return m;

        const double dx = region.sizeX() / m.nx;
        const double dy = region.sizeY() / m.ny;
        const double dz = region.sizeZ() / m.nz;
        for (int k = 0; k < m.nz; ++k)
        {
            const double z = region.zmin + (k + 0.5) * dz;
            for (int j = 0; j < m.ny; ++j)
            {
                const double y = region.ymin + (j + 0.5) * dy;
                for (int i = 0; i < m.nx; ++i)
                {
                    const double x = region.xmin + (i + 0.5) * dx;
                    if (nodeInside(*root_, x, y, z))
                        m.occ[(std::size_t(k) * m.ny + j) * m.nx + i] = 1;
                }
            }
        }
        return m;
    }

} // namespace waveguide
