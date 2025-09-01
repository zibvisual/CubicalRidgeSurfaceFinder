#pragma once

#include "Vec.hpp"
#include "Dims.hpp"

class CornerBBox
{
public:
    CornerBBox(VecFloat min_corner, VecFloat max_corner) : m_min_corner(min_corner), m_max_corner(max_corner){}

    static CornerBBox fromCornerAndVoxelsize(VecFloat min, VecFloat voxelsize, Dims dims){
        return CornerBBox(min, min + static_cast<VecFloat>(dims) * voxelsize);
    }

    static CornerBBox fromCenterAndVoxelsize(VecFloat origin, VecFloat voxelsize, Dims dims){
        const auto min = origin - (voxelsize * 0.5f);
        return fromCornerAndVoxelsize(min, voxelsize, dims);
    }

    static CornerBBox fromCornerBbox(VecFloat min_corner, VecFloat max_corner){
        return CornerBBox(min_corner, max_corner);
    }

    static CornerBBox fromCenterBBox(VecFloat origin, VecFloat max_center, Dims dims){
        const auto voxelsize = VecFloat(
            dims[0] == 1 ? 0 : (max_center[0] - origin[0]) / static_cast<float>(dims[0] - 1),
            dims[1] == 1 ? 0 : (max_center[1] - origin[1]) / static_cast<float>(dims[1] - 1),
            dims[2] == 1 ? 0 : (max_center[2] - origin[2]) / static_cast<float>(dims[2] - 1)
        );
        return fromCenterAndVoxelsize(origin, voxelsize, dims);
    }

    /**
     * We assume the bounding box to include the voxels of each signal (and all voxels are equally sized (per dimension))
     */
    VecFloat voxelsize(Dims dims) const {
        return VecFloat(
            (m_max_corner - m_min_corner) / static_cast<VecFloat>(dims)
        );
    }

    VecFloat min_corner() const {
        return m_min_corner;
    }

    VecFloat max_corner() const {
        return m_max_corner;
    }

    VecFloat origin(Dims dims) const {
        return m_min_corner - (voxelsize(dims) * 0.5f);
    }

    bool operator==(const CornerBBox& rhs) const
    {
        return m_min_corner == rhs.m_min_corner && m_max_corner == rhs.m_max_corner;
    }

    bool equals(const CornerBBox& rhs, const float epsilon = 0.0001) const {
        return m_min_corner.equals(rhs.m_min_corner, epsilon) && m_max_corner.equals(rhs.m_max_corner, epsilon);
    }

private:
    VecFloat m_min_corner;
    VecFloat m_max_corner;
};

inline std::ostream& operator<<(std::ostream& os, const CornerBBox& bbox)
{
    os << "from " << bbox.min_corner() << " to " << bbox.max_corner();
    return os;
}

class CenterBBox
{
public:
    CenterBBox(VecFloat min_center, VecFloat max_center) : m_min_center(min_center), m_max_center(max_center){}

    static CenterBBox fromCornerAndVoxelsize(VecFloat min, VecFloat voxelsize, Dims dims){
        return CenterBBox(min + voxelsize * 0.5f, min + static_cast<VecFloat>(dims) * voxelsize - voxelsize * 0.5f);
    }

    static CenterBBox fromCenterAndVoxelsize(VecFloat origin, VecFloat voxelsize, Dims dims){
        return CenterBBox(origin, origin + (static_cast<VecFloat>(dims) - 1.0f) * voxelsize);
    }

    // static CenterBBox fromCornerBBox(VecFloat min_center, VecFloat max_center){
    //     return CenterBBox(min_center, max_center);
    // }

    static CenterBBox fromCenterBBox(VecFloat min_center, VecFloat max_center){
        return CenterBBox(min_center, max_center);
    }

    /**
     * We assume the bounding box to include the voxels of each signal (and all voxels are equally sized (per dimension))
     */
    VecFloat voxelsize(Dims dims) const {
        return VecFloat(
            (dims[0] == 1) ? 0.f : (m_max_center[0] - m_min_center[0]) / (dims[0] - 1),
            (dims[1] == 1) ? 0.f : (m_max_center[1] - m_min_center[1]) / (dims[1] - 1),
            (dims[2] == 1) ? 0.f : (m_max_center[2] - m_min_center[2]) / (dims[2] - 1)
        );
    }

    VecFloat min_center() const {
        return m_min_center;
    }

    VecFloat max_center() const {
        return m_max_center;
    }

    VecFloat origin(Dims dims) const {
        return m_min_center;
    }

    bool operator==(const CenterBBox& rhs) const
    {
        return m_min_center == rhs.m_min_center && m_max_center == rhs.m_max_center;
    }

    bool equals(const CenterBBox& rhs, const float epsilon = 0.0001) const {
        return m_min_center.equals(rhs.m_min_center, epsilon) && m_max_center.equals(rhs.m_max_center, epsilon);
    }

private:
    VecFloat m_min_center;
    VecFloat m_max_center;
};

inline std::ostream& operator<<(std::ostream& os, const CenterBBox& bbox)
{
    os << "from " << bbox.min_center() << " to " << bbox.max_center();
    return os;
}