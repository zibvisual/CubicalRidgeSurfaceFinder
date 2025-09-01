#pragma once

#include <cstdint>
#include <optional>

#include <utils/Dims.hpp>

class Direction
{
public:
  enum Value : uint8_t
  {
    // x-axis (+1,-1)
    RIGHT = 0,
    LEFT = 3,
    // y-axis (+1,-1)
    UP = 1,
    DOWN = 4,
    // z-axis (+1,-1)
    FORWARD = 2,
    BACKWARD = 5,
  };

  Direction() = default;
  constexpr Direction(Value dir) : value(dir) { }

  // Allow switch and comparisons.
  constexpr operator Value() const { return value; }

  // Prevent usage: if(fruit)
  explicit operator bool() const = delete;

//   constexpr bool operator==(Direction a) const { return value == a.value; }
//   constexpr bool operator!=(Direction a) const { return value != a.value; }
//   constexpr bool operator==(Direction::Value a) const { return value == a; }
//   constexpr bool operator!=(Direction::Value a) const { return value != a; }

  Direction flip() const{
    return static_cast<Direction::Value>((value + 3) % 6);
  }

  bool sign() const {
    return value < 3;
  }

  // returns 0 for x, 1 for y and 2 for z
  uint8_t abs() const {
    return value % 3;
  }

  int64_t c_offset(Dims dims) const {
    switch (value)
        {
            case Direction::LEFT:
                return -1;
            case Direction::RIGHT:
                return 1;
            case Direction::DOWN:
                return -dims[0];
            case Direction::UP:
                return dims[0];
            case Direction::BACKWARD:
                return -dims[0] * dims[1];
            case Direction::FORWARD:
                return dims[0] * dims[1];
            default:
                return 0;
        }
  }

  VecInt c_translation() const {
    switch (value)
    {
        case Direction::LEFT:
            return VecInt(-1, 0, 0);
        case Direction::RIGHT:
            return VecInt(1, 0, 0);
        case Direction::DOWN:
            return VecInt(0, -1, 0);
        case Direction::UP:
            return VecInt(0, 1, 0);
        case Direction::BACKWARD:
            return VecInt(0, 0, -1);
        case Direction::FORWARD:
            return VecInt(0, 0, 1);
        default:
            return VecInt(0, 0, 0);
    }
  }

  static std::optional<Direction> from_c_offset(int64_t offset, Dims dims){
        if (offset == -1)
            return Direction::LEFT;
        else if (offset == 1)
            return Direction::RIGHT;
        else if (offset == -dims[0])
            return Direction::DOWN;
        else if (offset == dims[0])
            return Direction::UP;
        else if (offset == -dims[0] * dims[1])
            return Direction::BACKWARD;
        else if (offset == dims[0] * dims[1])
            return Direction::FORWARD;
        else
            return std::optional<Direction>();
  }

private:
  Value value;
};