#pragma once

#include <utils/Vec.hpp>
#include <utils/Dims.hpp>

class BLocation {
    public:

    BLocation(VecSize location)
    : m_index((location.z() << 40) | (location.y() << 20) | location.x())
    {}

    BLocation(VecSize location, Dims dims)
        : m_index(dims.contains(location) ? (location.z() << 40) | (location.y() << 20) | location.x() : (static_cast<uint64_t>(0b1111) << 60))
    {}

    BLocation(VecInt location)
    : m_index((static_cast<uint64_t>(location.z()) << 40) | (static_cast<uint64_t>(location.y()) << 20) | static_cast<uint64_t>(location.x()))
    {}

    BLocation(VecInt location, Dims dims)
    : m_index(dims.contains(location) ? (static_cast<uint64_t>(location.z()) << 40) | (static_cast<uint64_t>(location.y()) << 20) | static_cast<uint64_t>(location.x()) : (static_cast<uint64_t>(0b1111) << 60))
    {}

    BLocation(uint64_t index)
        : m_index(index)
    {}

    // Empty BLocation (out of bounds)
    BLocation()
        : m_index((static_cast<uint64_t>(0b1111) << 60))
    {}

    bool operator==(const BLocation& rhs) const 
    {
        return m_index == rhs.m_index;
    }

    VecSize location_unsafe() const
    {
        return VecSize(
            m_index & 0xFFFFF,
            (m_index >> 20) & 0xFFFFF,
            (m_index >> 40) & 0xFFFFF
        );
    }

    std::optional<VecSize> location() const {
        if (empty()){
            return std::optional<VecSize>();
        }
        return location_unsafe();
    }

    uint64_t index() const
    {
        return m_index;
    }

    bool empty() const
    {
        return (m_index >> 60) & 0b1111 == 0b1111;
    }

    bool valid() const
    {
        return (m_index >> 60) & 0b1111 != 0b1111;
    }

private:
    uint64_t m_index;
};