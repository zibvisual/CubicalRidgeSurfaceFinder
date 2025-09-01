#pragma once

#include <cstddef>
#include "Vec.hpp"

class Dims
{
public:
    Dims(std::size_t x) : m_x(x), m_y(x), m_z(x) {}
    Dims(std::size_t x, std::size_t y, std::size_t z) : m_x(x), m_y(y), m_z(z) {}

    bool contains(std::size_t x, std::size_t y, std::size_t z) const
    {
        return x < m_x && y < m_y && z < m_z;
    }

    bool contains(int x, int y, int z) const
    {
        return x >= 0 && static_cast<std::size_t>(x) < m_x && y >= 0 && static_cast<std::size_t>(y) < m_y && z >= 0 && static_cast<std::size_t>(z) < m_z;
    }

    bool contains(VecSize loc) const
    {
        return contains(loc.x(), loc.y(), loc.z());
    }

    bool contains(VecInt loc) const
    {
        return contains(loc.x(), loc.y(), loc.z());
    }

    std::size_t x() const
    {
        return m_x;
    }

    std::size_t y() const
    {
        return m_y;
    }

    std::size_t z() const
    {
        return m_z;
    }

    std::size_t size() const {
        return m_x * m_y * m_z;
    }

    std::size_t& operator[](int index){
        if(index == 0){
            return m_x;
        }else if (index == 1){
            return m_y;
        }else if (index == 2){
            return m_z;
        }
        throw std::out_of_range("Invalid indexing");
    }

    const std::size_t& operator[](int index) const {
        if(index == 0){
            return m_x;
        }else if (index == 1){
            return m_y;
        }else if (index == 2){
            return m_z;
        }
        throw std::out_of_range("Invalid indexing");
    }

    bool operator==(const Dims& rhs) const {
        return m_x == rhs.m_x && m_y == rhs.m_y && m_z == rhs.m_z;
    }

    // VecFloat operator*(VecFloat const& v){
    //     return VecFloat(
    //         v.x() * static_cast<float>(m_x),
    //         v.y() * static_cast<float>(m_y),
    //         v.z() * static_cast<float>(m_z)
    //     );
    // }

    Dims extend(std::size_t size) const {
        return Dims(m_x + size, m_y + size, m_z + size);
    }

    explicit operator VecSize() const {
        return VecSize(
            static_cast<std::size_t>(m_x),
            static_cast<std::size_t>(m_y),
            static_cast<std::size_t>(m_z)
        );
    }
    explicit operator VecInt() const {
        return VecInt(
            static_cast<int>(m_x),
            static_cast<int>(m_y),
            static_cast<int>(m_z)
        );

    }
    explicit operator VecFloat() const {
        return VecFloat(
            static_cast<float>(m_x),
            static_cast<float>(m_y),
            static_cast<float>(m_z)
        );
    }

private:
    std::size_t m_x;
    std::size_t m_y;
    std::size_t m_z;
};