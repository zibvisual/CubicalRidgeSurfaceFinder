#pragma once

#include "Dims.hpp"
#include "Vec.hpp"
#include "BBox.hpp"

/**
 * Lattice in which the m_origin represents the exact location of the first "measured" signal.
 */
class Lattice {
public:
    Lattice() : m_dims(Dims(0,0,0)), m_origin(VecFloat(0.f,0.f,0.f)), m_voxelsize(VecFloat(1.f, 1.f, 1.f)){}
    Lattice(Dims dims) : m_dims(dims), m_origin(VecFloat(0.f,0.f,0.f)), m_voxelsize(VecFloat(1.f,1.f,1.f)){}
    Lattice(Dims dims, VecFloat origin, VecFloat voxelsize) : m_dims(dims), m_origin(origin), m_voxelsize(voxelsize){}
    Lattice(Dims dims, CornerBBox bbox) : m_dims(dims), m_origin(bbox.origin(dims)), m_voxelsize(bbox.voxelsize(dims)){}

    Dims dims() const {
        return m_dims;
    }

    std::size_t size() const {
        return m_dims.size();
    }

    VecFloat voxelsize() const {
        return m_voxelsize;
    }

    void setVoxelSize(VecFloat voxelsize) {
        m_voxelsize = voxelsize;
    }

    VecFloat origin() const {
        return m_origin;
    }

    void setOrigin(VecFloat orig) {
        m_origin = orig;
    }

    CornerBBox cornerbox() const {
        return CornerBBox::fromCenterAndVoxelsize(m_origin, m_voxelsize, m_dims);
    }

    CenterBBox centerbox() const {
        return CenterBBox::fromCenterAndVoxelsize(m_origin, m_voxelsize, m_dims);
    }

    void setBBox(CornerBBox bbox){
        m_voxelsize = bbox.voxelsize(m_dims);
        m_origin = bbox.origin(m_dims);
    }

    void setBBox(CenterBBox bbox){
        m_voxelsize = bbox.voxelsize(m_dims);
        m_origin = bbox.origin(m_dims);
    }

    VecFloat worldPosition(VecSize gridLocation) const {
        return static_cast<VecFloat>(gridLocation) * m_voxelsize + m_origin;
    }

    VecFloat cornerPosition(VecSize gridLocation) const {
        return (static_cast<VecFloat>(gridLocation) - 0.5f) * m_voxelsize + m_origin;
    }

    VecFloat gridPosition(VecFloat worldPosition) const {
        return (worldPosition - m_origin) / m_voxelsize;
    }

    VecSize gridLocation(VecFloat worldPosition) const {
        // 0.5f and floor act as round()
        return static_cast<VecSize>((gridPosition(worldPosition) + 0.5f).floor().max(0.0)).clamp(m_dims);
    }

    // These are tempporary functions and should be chaned later on. Instead of std::size_t, FieldLocations should be created.
    std::size_t c_index(VecSize gridLocation) const {
        return gridLocation.z() * m_dims.x() * m_dims.y() + gridLocation.y() * m_dims.x() + gridLocation.x();
    }  

    std::size_t c_index(VecFloat worldPosition) const {
        return c_index(gridLocation(worldPosition));
    }

    static VecSize gridLocationFromCIndex(std::size_t index, Dims dims){
        const uint64_t z = index / (dims.y() * dims.x());
        const uint64_t rest = index - z * (dims.y() * dims.x());
        const uint64_t y = rest / dims.x();
        const uint64_t x = rest - y * dims.x();
        return VecSize(x, y, z);
    }

private:
    Dims m_dims;
    VecFloat m_origin;
    VecFloat m_voxelsize;
};