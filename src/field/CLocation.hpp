#pragma once

#include <utils/Vec.hpp>
#include <utils/Dims.hpp>

class CLocation {
public:
    CLocation(VecSize location, Dims dims)
        : m_index(dims.contains(location) ? location.z() * dims.x() * dims.y() + location.y() * dims.x() + location.x() : std::numeric_limits<uint64_t>::max())
    {}

    CLocation(VecInt location, Dims dims)
    : m_index(dims.contains(location) ? static_cast<uint64_t>(location.z()) * dims.x() * dims.y() + static_cast<uint64_t>(location.y()) * dims.x() + static_cast<uint64_t>(location.x()) : std::numeric_limits<uint64_t>::max())
    {}

    CLocation(uint64_t index)
        : m_index(index)
    {}

    // Empty CLocation (out of bounds)
    CLocation()
        : m_index(std::numeric_limits<uint64_t>::max())
    {}

    bool operator==(const CLocation& rhs) const 
    {
        return m_index == rhs.m_index;
    }

    // THIS IS DANGEROUS (for empty field locations!)
    VecSize location_unsafe(Dims dims) const
    {
        const uint64_t z = m_index / (dims.y() * dims.x());
        const uint64_t rest = m_index - z * (dims.y() * dims.x());
        const uint64_t y = rest / dims.x();
        const uint64_t x = rest - y * dims.x();
        return VecSize(x, y, z);
    }

    std::optional<VecSize> location(Dims dims) const {
        if (dims.size() <= m_index){
            return std::optional<VecSize>();
        }
        return location_unsafe(dims);
    }

    uint64_t index() const
    {
        return m_index;
    }

    bool empty() const
    {
        return m_index == std::numeric_limits<uint64_t>::max();
    }

    bool valid() const
    {
        return m_index != std::numeric_limits<uint64_t>::max();
    }

private:
    uint64_t m_index;
};