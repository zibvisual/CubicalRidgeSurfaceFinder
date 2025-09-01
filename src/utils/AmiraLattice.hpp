#pragma once

#include "Dims.hpp"
#include "Vec.hpp"
#include "BBox.hpp"

class Lattice {
public:
    Lattice() : m_dims(Dims(0,0,0)), m_voxelsize(VecFloat(1.f, 1.f, 1.f)), m_origin(VecFloat(0.f,0.f,0.f)){}
    Lattice(Dims dims) : m_dims(dims), m_voxelsize(VecFloat(1.f,1.f,1.f)), m_origin(VecFloat(0.f,0.f,0.f)){}
    Lattice(Dims dims, VecFloat voxelsize, VecFloat origin) : m_dims(dims), m_voxelsize(voxelsize), m_origin(origin){}
    Lattice(Dims dims, CornerBBox bbox) : m_dims(dims), m_voxelsize(bbox.voxelsize(dims)), m_origin(bbox.origin()){}

    Dims dims() const {
        return m_dims;
    }

    std::size_t size() const {
        return m_dims.size();
    }

    VecFloat voxelsize() const {
        return m_voxelsize;
    }
    VecFloat origin() const {
        return m_origin;
    }

    CornerBBox bbox() const {
        return CornerBBox::fromCenterAndVoxelSize(m_origin, m_voxelsize, m_dims);
    }

    VecFloat worldPosition(VecSize gridLocation) const {
        auto offsets = VecFloat(
            (gridLocation[0] < 1 ? 0.0f : (gridLocation[0] >= m_dims[0] ? 1.0f : 0.5f)),
            (gridLocation[1] < 1 ? 0.0f : (gridLocation[1] >= m_dims[1] ? 1.0f : 0.5f)),
            (gridLocation[2] < 1 ? 0.0f : (gridLocation[2] >= m_dims[2] ? 1.0f : 0.5f))
        );
        return (static_cast<VecFloat>(gridLocation) - offsets) * m_voxelsize + m_origin;   
    }

    VecFloat gridPosition(VecFloat worldPosition) const {
        return (worldPosition - m_origin + (m_voxelsize * 0.5f)) / m_voxelsize;
    }

    VecSize gridLocation(VecFloat worldPosition) const {
        return static_cast<VecSize>(gridPosition(worldPosition).floor().max(0.0)).clamp(m_dims);
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
    VecFloat m_voxelsize;
    VecFloat m_origin;
};