#include "Vec.hpp"
#include <stdexcept>
#include <algorithm>

#include "Dims.hpp"

const VecInt VecInt::operator+(const VecInt &other) const
{
    return VecInt(
        m_values[0] + other.m_values[0],
        m_values[1] + other.m_values[1],
        m_values[2] + other.m_values[2]);
}

const VecInt VecInt::operator+(const int &other) const
{
    return VecInt(
        m_values[0] + other,
        m_values[1] + other,
        m_values[2] + other);
}

const VecInt VecInt::operator-() const
{
    return VecInt(
        -m_values[0],
        -m_values[1],
        -m_values[2]);
}

const VecInt VecInt::operator-(const VecInt &other) const
{
    return VecInt(
        m_values[0] - other.m_values[0],
        m_values[1] - other.m_values[1],
        m_values[2] - other.m_values[2]);
}

const VecInt VecInt::operator-(const int &other) const
{
    return VecInt(
        m_values[0] - other,
        m_values[1] - other,
        m_values[2] - other);
}

const VecFloat VecInt::operator*(const float &scalar) const
{
    return VecFloat(
        m_values[0] * scalar,
        m_values[1] * scalar,
        m_values[2] * scalar);
}

const VecInt VecInt::operator*(const int &scalar) const
{
    return VecInt(
        m_values[0] * scalar,
        m_values[1] * scalar,
        m_values[2] * scalar);
}

const VecInt VecInt::abs() const {
    return VecInt(
        std::abs(m_values[0]),
        std::abs(m_values[1]),
        std::abs(m_values[2])
    );
}

bool VecInt::operator==(const VecInt& other) const {
    return m_values[0] == other.m_values[0] && m_values[1] == other.m_values[1] && m_values[2] == other.m_values[2];
}
bool VecInt::operator!=(const VecInt& other) const {
    return m_values[0] != other.m_values[0] || m_values[1] != other.m_values[1] || m_values[2] != other.m_values[2];
}

std::optional<VecSize> VecInt::tryVecSize() const {
    if (m_values[0] < 0 || m_values[1] < 0 || m_values[2] < 0){
        return std::optional<VecSize>();
    }
    return static_cast<VecSize>(*this);
}

VecInt::operator VecSize() const {
    return VecSize(
        static_cast<std::size_t>(m_values[0]),
        static_cast<std::size_t>(m_values[1]),
        static_cast<std::size_t>(m_values[2])
    );
}

const VecInt VecInt::LEFT = {-1, 0, 0};
const VecInt VecInt::RIGHT = {1, 0, 0};
const VecInt VecInt::UP = {0, 1, 0};
const VecInt VecInt::DOWN = {0, -1, 0};
const VecInt VecInt::FORWARD = {0, 0, -1};
const VecInt VecInt::BACKWARD = {0, 0, 1};

VecSize VecSize::min(std::size_t min) const{
    return VecSize(
        std::min(m_values[0], min),
        std::min(m_values[1], min),
        std::min(m_values[2], min)
    );
}
VecSize VecSize::max(std::size_t max) const{
    return VecSize(
        std::max(m_values[0], max),
        std::max(m_values[1], max),
        std::max(m_values[2], max)
    );
}
VecSize VecSize::clamp(std::size_t min, std::size_t max) const{
    return VecSize(
        std::clamp(m_values[0], min, max),
        std::clamp(m_values[1], min, max),
        std::clamp(m_values[2], min, max)
    );
}

VecSize VecSize::min(VecSize min) const{
    return VecSize(
        std::min(m_values[0], min.m_values[0]),
        std::min(m_values[1], min.m_values[1]),
        std::min(m_values[2], min.m_values[2])
    );
}
VecSize VecSize::max(VecSize max) const{
    return VecSize(
        std::max(m_values[0], max.m_values[0]),
        std::max(m_values[1], max.m_values[1]),
        std::max(m_values[2], max.m_values[2])
    );
}
VecSize VecSize::clamp(VecSize min, VecSize max) const{
    return VecSize(
        std::clamp(m_values[0], min.m_values[0], max.m_values[0]),
        std::clamp(m_values[1], min.m_values[1], max.m_values[1]),
        std::clamp(m_values[2], min.m_values[2], max.m_values[2])
    );
}

VecSize VecSize::clamp(Dims dims) const{
    return VecSize(
        std::min(m_values[0], dims[0]),
        std::min(m_values[1], dims[1]),
        std::min(m_values[2], dims[2])
    );
}

bool VecSize::operator==(const VecSize& rhs) const {
    return m_values[0] == rhs.m_values[0] && m_values[1] == rhs.m_values[1] && m_values[2] == rhs.m_values[2];
}

VecSize::operator VecInt() const {
    return VecInt(
        static_cast<int>(m_values[0]),
        static_cast<int>(m_values[1]),
        static_cast<int>(m_values[2])
    );
}

VecSize::operator VecFloat() const {
    return VecFloat(
        static_cast<float>(m_values[0]),
        static_cast<float>(m_values[1]),
        static_cast<float>(m_values[2])
    );
}

VecFloat VecFloat::floor() const
{
    return VecFloat(
        std::floor(m_values[0]),
        std::floor(m_values[1]),
        std::floor(m_values[2])  
    );
}

VecFloat VecFloat::min(float min) const
{
    return VecFloat(
        std::min(m_values[0], min),
        std::min(m_values[1], min),
        std::min(m_values[2], min)
    );
}

VecFloat VecFloat::max(float max) const
{
    return VecFloat(
        std::max(m_values[0], max),
        std::max(m_values[1], max),
        std::max(m_values[2], max)
    );
}

VecFloat VecFloat::clamp(float min, float max) const
{
    return VecFloat(
        std::clamp(m_values[0], min, max),
        std::clamp(m_values[1], min, max),
        std::clamp(m_values[2], min, max)
    );
}

bool VecFloat::operator==(const VecFloat& rhs) const {
    return m_values[0] == rhs.m_values[0] && m_values[1] == rhs.m_values[1] && m_values[2] == rhs.m_values[2];
}

bool VecFloat::equals(const VecFloat& rhs, const float epsilon /*= 0.0001*/) const {
    return abs(m_values[0] - rhs.m_values[0]) + abs(m_values[1] - rhs.m_values[1]) + abs(m_values[2] - rhs.m_values[2]) < epsilon;
}

VecFloat VecFloat::operator+(VecFloat const& other) const {
    return VecFloat(
        other.x() + m_values[0],
        other.y() + m_values[1],
        other.z() + m_values[2]
    );
}

const VecFloat VecFloat::operator-(const VecFloat &other) const
{
    return VecFloat(
        this->m_values[0] - other.m_values[0],
        this->m_values[1] - other.m_values[1],
        this->m_values[2] - other.m_values[2]);
}


const VecFloat VecFloat::operator*(const float &other) const
{
    return VecFloat(
        this->m_values[0] * other,
        this->m_values[1] * other,
        this->m_values[2] * other);
}

const VecFloat VecFloat::operator*(const VecFloat &other) const
{
    return VecFloat(
        this->m_values[0] * other.m_values[0],
        this->m_values[1] * other.m_values[1],
        this->m_values[2] * other.m_values[2]);
}

const VecFloat VecFloat::operator/(const VecFloat &other) const
{
    return VecFloat(
        this->m_values[0] / other.m_values[0],
        this->m_values[1] / other.m_values[1],
        this->m_values[2] / other.m_values[2]);
}

VecFloat::operator VecInt() const {
    return VecInt(
        static_cast<int>(m_values[0]),
        static_cast<int>(m_values[1]),
        static_cast<int>(m_values[2])
    );
}

VecFloat::operator VecSize() const {
    return VecSize(
        static_cast<std::size_t>(m_values[0]),
        static_cast<std::size_t>(m_values[1]),
        static_cast<std::size_t>(m_values[2])
    );
}

bool VecFloat::lexicographic_order_exact(const VecFloat& other) const {
    if(m_values[0] < other.m_values[0]){
        return true;
    }else if(m_values[0] > other.m_values[0]){
        return false;
    }

    if(m_values[1] < other.m_values[1]){
        return true;
    }else if (m_values[1] > other.m_values[1]){
        return false;
    }

    if(m_values[2] < other.m_values[2]){
        return true;
    }else if (m_values[2] > other.m_values[2]){
        return false;
    }

    return false;
}

// This is a valid compare function for sorting (and the translation rule for equiv expression holds true)
bool VecFloat::lexicographic_order_epsilon(const VecFloat& other, const float epsilon /*= 0.0001*/) const {
    const float reciprocal = 1.0/epsilon;

    const float ax = round(m_values[0] * reciprocal);
    const float bx = round(other.m_values[0] * reciprocal);
    if(ax < bx){
        return true;
    }else if(ax > bx){
        return false;
    }

    const float ay = round(m_values[1] * reciprocal);
    const float by = round(other.m_values[1] * reciprocal);
    if(ay < by){
        return true;
    }else if(ay > by){
        return false;
    }

    const float az = round(m_values[2] * reciprocal);
    const float bz = round(other.m_values[2] * reciprocal);
    if(az < bz){
        return true;
    }else if(az > bz){
        return false;
    }

    return false;
}