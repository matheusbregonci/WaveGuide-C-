#pragma once

#include <cstdint>
#include <memory>
#include <vector>

// Constructive Solid Geometry (CSG) for arbitrary waveguide / cavity shapes.
//
// A shape is a tree of primitives (axis-aligned box, z-aligned cylinder)
// combined with boolean operators (union / difference / intersection). This is
// the geometry substrate for the future numerical solver (Phase 3): a T-junction
// is just box(main) ∪ box(stub); an iris is box − box; etc.
//
// The two things the rest of the program needs from a geometry are:
//   - inside(x,y,z): is a MODEL-space point solid (inside the metal walls)?
//   - bounds():      the axis-aligned bounding box, for grids and framing.
// voxelize() rasterizes the shape into a boolean mask for a grid-based solver.
//
// All coordinates are in meters, in the same model frame the field models use.

namespace waveguide
{

    struct Aabb
    {
        double xmin, ymin, zmin;
        double xmax, ymax, zmax;
        double sizeX() const { return xmax - xmin; }
        double sizeY() const { return ymax - ymin; }
        double sizeZ() const { return zmax - zmin; }
        bool valid() const { return xmax >= xmin && ymax >= ymin && zmax >= zmin; }
    };

    // Boolean occupancy grid over an AABB (1 = solid). Ready to hand to an
    // FDFD/FDTD/FEM assembler as the material mask.
    struct VoxelMask
    {
        int nx = 0, ny = 0, nz = 0;
        Aabb box{};
        std::vector<std::uint8_t> occ; // size nx*ny*nz, row-major (x fastest)

        std::uint8_t at(int i, int j, int k) const
        {
            if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) return 0;
            return occ[(std::size_t(k) * ny + j) * nx + i];
        }
        std::size_t solidCount() const;
    };

    class Geometry
    {
    public:
        Geometry() = default;

        // Axis-aligned box: center (cx,cy,cz) and FULL sizes (sx,sy,sz).
        static Geometry box(double cx, double cy, double cz,
                            double sx, double sy, double sz);
        // Circular cylinder along z: axis at (cx,cy), spanning z in [z0,z1].
        static Geometry cylinderZ(double cx, double cy,
                                  double z0, double z1, double radius);

        // Boolean composition (value semantics — returns a new Geometry).
        Geometry unite(const Geometry& other) const;
        Geometry subtract(const Geometry& other) const;
        Geometry intersect(const Geometry& other) const;

        bool empty() const { return root_ == nullptr; }

        // Is the point solid (inside the shape)?
        bool inside(double x, double y, double z) const;

        // Axis-aligned bounding box of the shape.
        Aabb bounds() const;

        // Rasterize into an nx*ny*nz occupancy grid over `region` (defaults to
        // bounds()). A cell is solid if its center is inside().
        VoxelMask voxelize(int nx, int ny, int nz) const;
        VoxelMask voxelize(int nx, int ny, int nz, const Aabb& region) const;

        // Opaque CSG node; its definition lives in Geometry.cpp.
        struct Node;

    private:
        explicit Geometry(std::shared_ptr<const Node> root) : root_(std::move(root)) {}
        std::shared_ptr<const Node> root_;
    };

} // namespace waveguide
