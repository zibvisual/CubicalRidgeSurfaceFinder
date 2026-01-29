#pragma once

#include <cstddef>
#include <cmath>
#include <ostream>
#include <optional>
#include <array>

class Dims;
class VecFloat;
class VecSize;

class VecInt
{
public:
    using value_type = int;

    static const VecInt LEFT;
    static const VecInt RIGHT;
    static const VecInt UP;
    static const VecInt DOWN;
    static const VecInt FORWARD;
    static const VecInt BACKWARD;

    VecInt() : m_values({0,0,0}) {}
    VecInt(int x) : m_values({x,x,x}) {}
    VecInt(int x, int y, int z) : m_values({x,y,z}) {}

    inline int x() const {return m_values[0];}
    inline int y() const {return m_values[1];}
    inline int z() const {return m_values[2];}

    const VecInt operator+(const VecInt &other) const;
    const VecInt operator+(const int &other) const;
    const VecInt operator-() const;
    const VecInt operator-(const VecInt &other) const;
    const VecInt operator-(const int &other) const;
    const VecFloat operator*(const float &scalar) const;
    const VecInt operator*(const int &scalar) const;

    const VecInt abs() const;

    inline int& operator[](int index) {return m_values[index];}
    inline const int& operator[](int index) const {return m_values[index];}

    bool operator==(const VecInt& other) const;
    bool operator!=(const VecInt& other) const;

    // fallible conversion
    std::optional<VecSize> tryVecSize() const;

    explicit operator VecSize() const;

private:
    std::array<int, 3> m_values;
};

// const VecInt VecInt::LEFT = {-1, 0, 0};
// const VecInt VecInt::RIGHT = {1, 0, 0};
// const VecInt VecInt::UP = {0, 1, 0};
// const VecInt VecInt::DOWN = {0, -1, 0};
// const VecInt VecInt::FORWARD = {0, 0, -1};
// const VecInt VecInt::BACKWARD = {0, 0, 1};

class VecSize
{
public:
    using value_type = std::size_t;

    VecSize() : m_values({0,0,0}) {}
    VecSize(std::size_t x) : m_values({x,x,x}) {}
    VecSize(std::size_t x, std::size_t y, std::size_t z) : m_values({x,y,z}) {}

    inline std::size_t x() const {return m_values[0];}
    inline std::size_t y() const {return m_values[1];}
    inline std::size_t z() const {return m_values[2];}

    VecSize min(std::size_t min) const;
    VecSize max(std::size_t max) const;
    VecSize clamp(std::size_t min, std::size_t max) const;

    VecSize min(VecSize min) const;
    VecSize max(VecSize max) const;
    VecSize clamp(VecSize min, VecSize max) const;

    VecSize clamp(Dims dims) const;

    inline std::size_t& operator[](int index) {return m_values[index];}
    inline const std::size_t& operator[](int index) const {return m_values[index];}

    bool operator==(const VecSize& rhs) const;

    explicit operator VecInt() const;
    explicit operator VecFloat() const;

private:
    std::array<std::size_t, 3> m_values;
};

class VecFloat
{
public:
    using value_type = float;

    VecFloat() : m_values({0,0,0}) {}
    VecFloat(float x) : m_values({x,x,x}) {}
    VecFloat(float x, float y, float z) : m_values({x,y,z}) {}

    inline float x() const {return m_values[0];}
    inline float y() const {return m_values[1];}
    inline float z() const {return m_values[2];}

    VecFloat floor() const;
    VecFloat min(float min) const;
    VecFloat max(float max) const;
    VecFloat clamp(float min, float max) const;

    inline float& operator[](int index) {return m_values[index];}
    inline const float& operator[](int index) const {return m_values[index];}

    bool operator==(const VecFloat& rhs) const;
    bool equals(const VecFloat& rhs, const float epsilon = 0.0001) const;

    VecFloat operator+(VecFloat const& other) const;
    const VecFloat operator-(const VecFloat &other) const;
    const VecFloat operator*(const float &other) const;
    const VecFloat operator*(const VecFloat &other) const;
    const VecFloat operator/(const VecFloat &other) const;

    explicit operator VecInt() const;
    explicit operator VecSize() const;

    bool lexicographic_order_exact(const VecFloat& other) const;
    // This is a valid compare function for sorting (and the translation rule for equiv expression holds true)
    bool lexicographic_order_epsilon(const VecFloat& other, const float epsilon = 0.0001) const;

private:
    std::array<float, 3> m_values;
};


inline std::ostream& operator<<(std::ostream& os, const VecSize& vec)
{
    os << "("<< vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
    return os;
}
inline std::ostream& operator<<(std::ostream& os, const VecInt& vec)
{
    os << "("<< vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
    return os;
}
inline std::ostream& operator<<(std::ostream& os, const VecFloat& vec)
{
    os << "("<< vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
    return os;
}

template<>
struct std::hash<VecFloat>
{
    std::size_t operator()(const VecFloat& vec) const
    {
        return ((std::hash<float>()(vec.x())
        ^ (std::hash<float>()(vec.y()) << 1)) >> 1)
        ^ (std::hash<float>()(vec.z()) << 1);
    }
};