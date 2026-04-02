#pragma once

#include <cstddef>
#include <optional>
#include "Vec.hpp"

class Dims
{
public:
    Dims() : m_x(0), m_y(0), m_z(0) {}
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

    bool empty() const
    {
        return m_x == 0 && m_y == 0 && m_z == 0;
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

    size_t c_index_unsafe(int x, int y, int z) const {
        return static_cast<size_t>(z) * m_x * m_y + static_cast<size_t>(y) * m_x + static_cast<size_t>(x);
    }

    std::optional<size_t> c_index(int x, int y, int z) const {
        return contains(x,y,z) ? std::optional<size_t>(c_index_unsafe(x,y,z)) : std::nullopt;
    }

    size_t c_index_unsafe(size_t x, size_t y, size_t z) const {
        return z * m_x * m_y + y * m_x + x;
    }

    std::optional<size_t> c_index(size_t x, size_t y, size_t z) const {
        return contains(x,y,z) ? std::optional<size_t>(c_index_unsafe(x,y,z)) : std::nullopt;
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

inline std::ostream& operator<<(std::ostream& os, const Dims& dims)
{
    os << "["<< dims.x() << ", " << dims.y() << ", " << dims.z() << "]";
    return os;
}