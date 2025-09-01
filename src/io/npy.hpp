/*
   Copyright 2017-2023 Leon Merten Lohse
   Modified 2025 by Nicolas Klenert

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef NPY_HPP_
#define NPY_HPP_

#include <algorithm>
#include <array>
#include <complex>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>
#include <filesystem>
#include "common.hpp"

namespace npy
{

/* Compile-time test for byte order.
   If your compiler does not define these per default, you may want to define
   one of these constants manually.
   Defaults to little endian order. */
#if defined(__BYTE_ORDER) && __BYTE_ORDER == __BIG_ENDIAN || defined(__BIG_ENDIAN__) || defined(__ARMEB__) || \
    defined(__THUMBEB__) || defined(__AARCH64EB__) || defined(_MIBSEB) || defined(__MIBSEB) || defined(__MIBSEB__)
  const bool const_big_endian = true;
#else
  const bool const_big_endian = false;
#endif

  const size_t magic_string_length = 6;
  const std::array<char, magic_string_length> magic_string = {'\x93', 'N', 'U', 'M', 'P', 'Y'};

  const char little_endian_char = '<';
  const char big_endian_char = '>';
  const char no_endian_char = '|';

  constexpr std::array<char, 3> endian_chars = {little_endian_char, big_endian_char, no_endian_char};
  constexpr std::array<char, 4> numtype_chars = {'f', 'i', 'u', 'c'};

  constexpr char host_endian_char = (const_big_endian ? big_endian_char : little_endian_char);

  /* npy array length */
  using ndarray_len_t = unsigned long int;

  class shape_t
  {
  public:
    shape_t() : m_data(){}
    shape_t(ndarray_len_t x, ndarray_len_t y, ndarray_len_t z) : m_data({z,y,x}) {}  

    inline void push_back(ndarray_len_t item)
    {
      m_data.push_back(item);
    }

    inline ndarray_len_t size() const
    {
      ndarray_len_t size = 1;
      for (ndarray_len_t i : m_data)
        size *= i;

      return size;
    }

    inline ndarray_len_t ndims() const
    {
      return static_cast<ndarray_len_t>(m_data.size());
    }

    inline const std::vector<ndarray_len_t> &raw() const
    {
      return m_data;
    }

    inline ndarray_len_t &operator[](int index)
    {
      return m_data[index];
    }
  private:
    std::vector<ndarray_len_t> m_data;
  };

  inline std::ostream &operator<<(std::ostream &os, const shape_t &shape)
  {
    os << "[";
    for (auto item : shape.raw())
    {
      os << item << ", ";
    }
    os << "]";
    return os;
  }

  using version_t = std::pair<char, char>;

  enum class stype_t
  {
    invalid,
    int8,
    int16,
    int32,
    int64,
    uint8,
    uint16,
    uint32,
    uint64,
    f32,
    f64,
  };

  inline std::ostream &operator<<(std::ostream &os, const stype_t &stype)
  {
    switch (stype)
    {
    case stype_t::int8:
      os << "<int8>";
      break;
    case stype_t::int16:
      os << "<int16>";
      break;
    case stype_t::int32:
      os << "<int32>";
      break;
    case stype_t::int64:
      os << "<int64>";
      break;
    case stype_t::uint8:
      os << "<uint8>";
      break;
    case stype_t::uint16:
      os << "<uint16>";
      break;
    case stype_t::uint32:
      os << "<uint32>";
      break;
    case stype_t::uint64:
      os << "<uint64>";
      break;
    case stype_t::f32:
      os << "<f32>";
      break;
    case stype_t::f64:
      os << "<f64>";
      break;
    default:
      os << "<invalid type>";
    }
    return os;
  }

  template <typename T>
  struct stype_t_helper
  {
    static const stype_t val = stype_t::invalid;
  };

  template <>
  struct stype_t_helper<int8_t>
  {
    static const stype_t val = stype_t::int8;
  };

  template <>
  struct stype_t_helper<int16_t>
  {
    static const stype_t val = stype_t::int16;
  };

  template <>
  struct stype_t_helper<int32_t>
  {
    static const stype_t val = stype_t::int32;
  };

  template <>
  struct stype_t_helper<int64_t>
  {
    static const stype_t val = stype_t::int64;
  };

  template <>
  struct stype_t_helper<uint8_t>
  {
    static const stype_t val = stype_t::uint8;
  };

  template <>
  struct stype_t_helper<uint16_t>
  {
    static const stype_t val = stype_t::uint16;
  };

  template <>
  struct stype_t_helper<uint32_t>
  {
    static const stype_t val = stype_t::uint32;
  };

  template <>
  struct stype_t_helper<uint64_t>
  {
    static const stype_t val = stype_t::uint64;
  };

  template <>
  struct stype_t_helper<float>
  {
    static const stype_t val = stype_t::f32;
  };

  template <>
  struct stype_t_helper<double>
  {
    static const stype_t val = stype_t::f64;
  };

  /**
   * Check if casting does not (potentially) reduce precision
   */
  inline bool upcast(stype_t from, stype_t to)
  {
    switch (from)
    {
    case stype_t::int8:
      return to == stype_t::int8 || to == stype_t::int16 || to == stype_t::int32 || to == stype_t::int64;
    case stype_t::int16:
      return to == stype_t::int16 || to == stype_t::int32 || to == stype_t::int64;
    case stype_t::int32:
      return to == stype_t::int32 || to == stype_t::int64;
    case stype_t::int64:
      return to == stype_t::int64;
    case stype_t::uint8:
      return to == stype_t::int16 || to == stype_t::int32 || to == stype_t::int64 || to == stype_t::uint8 || to == stype_t::uint16 || to == stype_t::uint32 || to == stype_t::uint64;
    case stype_t::uint16:
      return to == stype_t::int32 || to == stype_t::int64 || to == stype_t::uint16 || to == stype_t::uint32 || to == stype_t::uint64;
    case stype_t::uint32:
      return to == stype_t::int64 || to == stype_t::uint32 || to == stype_t::uint64;
    case stype_t::uint64:
      return to == stype_t::uint64;
    case stype_t::f32:
      return to == stype_t::f32 || to == stype_t::f64;
    case stype_t::f64:
      return to == stype_t::f64;
    default:
      return false;
    }
  }

  inline stype_t from_raw(char kind, unsigned int itemsize)
  {
    switch (kind)
    {
    case 'f':
      switch (itemsize)
      {
      case 4:
        return stype_t::f32;
      case 8:
        return stype_t::f64;
      }
      break;
    case 'i':
      switch (itemsize)
      {
      case 1:
        return stype_t::int8;
      case 2:
        return stype_t::int16;
      case 4:
        return stype_t::int32;
      case 8:
        return stype_t::int64;
      }
      break;
    case 'u':
      switch (itemsize)
      {
      case 1:
        return stype_t::uint8;
      case 2:
        return stype_t::uint16;
      case 4:
        return stype_t::uint32;
      case 8:
        return stype_t::uint64;
      }
      break;
    }
    return stype_t::invalid;
  }

  inline char stype_kind(stype_t stype)
  {
    switch (stype)
    {
    case stype_t::int8:
      return 'i';
    case stype_t::int16:
      return 'i';
    case stype_t::int32:
      return 'i';
    case stype_t::int64:
      return 'i';
    case stype_t::uint8:
      return 'u';
    case stype_t::uint16:
      return 'u';
    case stype_t::uint32:
      return 'u';
    case stype_t::uint64:
      return 'u';
    case stype_t::f32:
      return 'f';
    case stype_t::f64:
      return 'f';
    default:
      throw std::runtime_error("invalid stype");
    }
  }

  inline unsigned int stype_itemsize(stype_t stype)
  {
    switch (stype)
    {
    case stype_t::int8:
      return 1;
    case stype_t::int16:
      return 2;
    case stype_t::int32:
      return 4;
    case stype_t::int64:
      return 8;
    case stype_t::uint8:
      return 1;
    case stype_t::uint16:
      return 2;
    case stype_t::uint32:
      return 4;
    case stype_t::uint64:
      return 8;
    case stype_t::f32:
      return 4;
    case stype_t::f64:
      return 8;
    default:
      throw std::runtime_error("invalid stype");
    }
  }

  inline char stype_byteorder(stype_t stype, bool big_endian)
  {
    if (stype == stype_t::int8 || stype == stype_t::uint8)
    {
      return no_endian_char;
    }
    return big_endian ? big_endian_char : little_endian_char;
  }

  struct dtype_t
  {
    stype_t stype;
    bool big_endian;

    dtype_t(stype_t stype, bool big_endian = const_big_endian) : stype(stype), big_endian(big_endian) {}

    inline std::string str() const
    {
      std::stringstream ss;
      ss << stype_byteorder(stype, big_endian) << stype_kind(stype) << stype_itemsize(stype);
      return ss.str();
    }
  };

  // struct dtype_t
  // {
  //   char byteorder;
  //   char kind;
  //   unsigned int itemsize;

  //   dtype_t(char byteorder, char kind, unsigned int itemsize)
  //       : byteorder(byteorder), kind(kind), itemsize(itemsize)
  //   {
  //   }

  //   inline std::string str() const
  //   {
  //     std::stringstream ss;
  //     ss << byteorder << kind << itemsize;
  //     return ss.str();
  //   }

  //   inline std::tuple<const char, const char, const unsigned int> tie() const
  //   {
  //     return std::tie(byteorder, kind, itemsize);
  //   }

  //   inline stype_t stype_t() const
  //   {
  //     return from_raw(kind, itemsize);
  //   }
  // };

  struct header_t
  {
    dtype_t dtype;
    bool fortran_order;
    shape_t shape;
  };

  inline void write_magic(std::ostream &ostream, version_t version)
  {
    ostream.write(magic_string.data(), magic_string_length);
    ostream.put(version.first);
    ostream.put(version.second);
  }

  inline version_t read_magic(std::istream &istream)
  {
    std::array<char, magic_string_length + 2> buf{};
    istream.read(buf.data(), sizeof(buf));

    if (!istream)
    {
      throw std::runtime_error("io error: failed reading file");
    }

    if (!std::equal(magic_string.begin(), magic_string.end(), buf.begin()))
      throw std::runtime_error("this file does not have a valid npy format.");

    version_t version;
    version.first = buf[magic_string_length];
    version.second = buf[magic_string_length + 1];

    return version;
  }

  // helpers
  inline bool is_digits(const std::string &str) { return std::all_of(str.begin(), str.end(), ::isdigit); }

  template <typename T, size_t N>
  inline bool in_array(T val, const std::array<T, N> &arr)
  {
    return std::find(std::begin(arr), std::end(arr), val) != std::end(arr);
  }

  inline dtype_t parse_descr(std::string typestring)
  {
    if (typestring.length() < 3)
    {
      throw std::runtime_error("invalid length");
    }

    char byteorder_c = typestring.at(0);
    char kind_c = typestring.at(1);
    std::string itemsize_s = typestring.substr(2);

    if (!in_array(byteorder_c, endian_chars))
    {
      throw std::runtime_error("invalid byteorder");
    }

    if (!in_array(kind_c, numtype_chars))
    {
      throw std::runtime_error("invalid kind");
    }

    if (!is_digits(itemsize_s))
    {
      throw std::runtime_error("invalid itemsize");
    }
    unsigned int itemsize = std::stoul(itemsize_s);

    return dtype_t(from_raw(kind_c, itemsize), byteorder_c == big_endian_char);
  }

  namespace pyparse
  {

    /**
      Removes leading and trailing whitespaces
      */
    inline std::string trim(const std::string &str)
    {
      const std::string whitespace = " \t";
      auto begin = str.find_first_not_of(whitespace);

      if (begin == std::string::npos)
        return "";

      auto end = str.find_last_not_of(whitespace);

      return str.substr(begin, end - begin + 1);
    }

    inline std::string get_value_from_map(const std::string &mapstr)
    {
      size_t sep_pos = mapstr.find_first_of(":");
      if (sep_pos == std::string::npos)
        return "";

      std::string tmp = mapstr.substr(sep_pos + 1);
      return trim(tmp);
    }

    /**
       Parses the string representation of a Python dict

       The keys need to be known and may not appear anywhere else in the data.
     */
    inline std::unordered_map<std::string, std::string> parse_dict(std::string in, const std::vector<std::string> &keys)
    {
      std::unordered_map<std::string, std::string> map;

      if (keys.size() == 0)
        return map;

      in = trim(in);

      // unwrap dictionary
      if ((in.front() == '{') && (in.back() == '}'))
        in = in.substr(1, in.length() - 2);
      else
        throw std::runtime_error("Not a Python dictionary.");

      std::vector<std::pair<size_t, std::string>> positions;

      for (auto const &value : keys)
      {
        size_t pos = in.find("'" + value + "'");

        if (pos == std::string::npos)
          throw std::runtime_error("Missing '" + value + "' key.");

        std::pair<size_t, std::string> position_pair{pos, value};
        positions.push_back(position_pair);
      }

      // sort by position in dict
      std::sort(positions.begin(), positions.end());

      for (size_t i = 0; i < positions.size(); ++i)
      {
        std::string raw_value;
        size_t begin{positions[i].first};
        size_t end{std::string::npos};

        std::string key = positions[i].second;

        if (i + 1 < positions.size())
          end = positions[i + 1].first;

        raw_value = in.substr(begin, end - begin);

        raw_value = trim(raw_value);

        if (raw_value.back() == ',')
          raw_value.pop_back();

        map[key] = get_value_from_map(raw_value);
      }

      return map;
    }

    /**
      Parses the string representation of a Python boolean
      */
    inline bool parse_bool(const std::string &in)
    {
      if (in == "True")
        return true;
      if (in == "False")
        return false;

      throw std::runtime_error("Invalid python boolan.");
    }

    /**
      Parses the string representation of a Python str
      */
    inline std::string parse_str(const std::string &in)
    {
      if ((in.front() == '\'') && (in.back() == '\''))
        return in.substr(1, in.length() - 2);

      throw std::runtime_error("Invalid python string.");
    }

    /**
      Parses the string represenatation of a Python tuple into a vector of its items
     */
    inline std::vector<std::string> parse_tuple(std::string in)
    {
      std::vector<std::string> v;
      const char seperator = ',';

      in = trim(in);

      if ((in.front() == '(') && (in.back() == ')'))
        in = in.substr(1, in.length() - 2);
      else
        throw std::runtime_error("Invalid Python tuple.");

      std::istringstream iss(in);

      for (std::string token; std::getline(iss, token, seperator);)
      {
        v.push_back(token);
      }

      return v;
    }

    template <typename T>
    inline std::string write_tuple(const std::vector<T> &v)
    {
      if (v.size() == 0)
        return "()";

      std::ostringstream ss;
      ss.imbue(std::locale("C"));

      if (v.size() == 1)
      {
        ss << "(" << v.front() << ",)";
      }
      else
      {
        const std::string delimiter = ", ";
        // v.size() > 1
        ss << "(";
        std::copy(v.begin(), v.end() - 1, std::ostream_iterator<T>(ss, delimiter.c_str()));
        ss << v.back();
        ss << ")";
      }

      return ss.str();
    }

    inline std::string write_boolean(bool b)
    {
      if (b)
        return "True";
      else
        return "False";
    }

  } // namespace pyparse

  inline header_t parse_header(std::string header)
  {
    /*
       The first 6 bytes are a magic string: exactly "x93NUMPY".
       The next 1 byte is an unsigned byte: the major version number of the file
       format, e.g. x01. The next 1 byte is an unsigned byte: the minor version
       number of the file format, e.g. x00. Note: the version of the file format
       is not tied to the version of the numpy package. The next 2 bytes form a
       little-endian unsigned short int: the length of the header data HEADER_LEN.
       The next HEADER_LEN bytes form the header data describing the array's
       format. It is an ASCII string which contains a Python literal expression of
       a dictionary. It is terminated by a newline ('n') and padded with spaces
       ('x20') to make the total length of the magic string + 4 + HEADER_LEN be
       evenly divisible by 16 for alignment purposes. The dictionary contains
       three keys:

       "descr" : dtype.descr
       An object that can be passed as an argument to the numpy.dtype()
       constructor to create the array's dtype. "fortran_order" : bool Whether the
       array data is Fortran-contiguous or not. Since Fortran-contiguous arrays
       are a common form of non-C-contiguity, we allow them to be written directly
       to disk for efficiency. "shape" : tuple of int The shape of the array. For
       repeatability and readability, this dictionary is formatted using
       pprint.pformat() so the keys are in alphabetic order.
     */

    // remove trailing newline
    if (header.back() != '\n')
      throw std::runtime_error("invalid header");
    header.pop_back();

    // parse the dictionary
    std::vector<std::string> keys{"descr", "fortran_order", "shape"};
    auto dict_map = npy::pyparse::parse_dict(header, keys);

    if (dict_map.size() == 0)
      throw std::runtime_error("invalid dictionary in header");

    std::string descr_s = dict_map["descr"];
    std::string fortran_s = dict_map["fortran_order"];
    std::string shape_s = dict_map["shape"];

    std::string descr = npy::pyparse::parse_str(descr_s);
    dtype_t dtype = parse_descr(descr);

    // convert literal Python bool to C++ bool
    bool fortran_order = npy::pyparse::parse_bool(fortran_s);

    // parse the shape tuple
    auto shape_v = npy::pyparse::parse_tuple(shape_s);

    shape_t shape;
    for (auto item : shape_v)
    {
      auto dim = static_cast<ndarray_len_t>(std::stoul(item));
      shape.push_back(dim);
    }

    return {dtype, fortran_order, shape};
  }

  inline std::string write_header_dict(const std::string &descr, bool fortran_order, const shape_t &shape)
  {
    std::string s_fortran_order = npy::pyparse::write_boolean(fortran_order);
    std::string shape_s = npy::pyparse::write_tuple(shape.raw());

    return "{'descr': '" + descr + "', 'fortran_order': " + s_fortran_order + ", 'shape': " + shape_s + ", }";
  }

  inline void write_header(std::ostream &out, const header_t &header)
  {
    std::string header_dict = write_header_dict(header.dtype.str(), header.fortran_order, header.shape);

    size_t length = magic_string_length + 2 + 2 + header_dict.length() + 1;

    version_t version{1, 0};
    if (length >= 255 * 255)
    {
      length = magic_string_length + 2 + 4 + header_dict.length() + 1;
      version = {2, 0};
    }
    size_t padding_len = 16 - length % 16;
    std::string padding(padding_len, ' ');

    // write magic
    write_magic(out, version);

    // write header length
    if (version == version_t{1, 0})
    {
      auto header_len = static_cast<uint16_t>(header_dict.length() + padding.length() + 1);

      std::array<uint8_t, 2> header_len_le16{static_cast<uint8_t>((header_len >> 0) & 0xff),
                                             static_cast<uint8_t>((header_len >> 8) & 0xff)};
      out.write(reinterpret_cast<char *>(header_len_le16.data()), 2);
    }
    else
    {
      auto header_len = static_cast<uint32_t>(header_dict.length() + padding.length() + 1);

      std::array<uint8_t, 4> header_len_le32{
          static_cast<uint8_t>((header_len >> 0) & 0xff), static_cast<uint8_t>((header_len >> 8) & 0xff),
          static_cast<uint8_t>((header_len >> 16) & 0xff), static_cast<uint8_t>((header_len >> 24) & 0xff)};
      out.write(reinterpret_cast<char *>(header_len_le32.data()), 4);
    }

    out << header_dict << padding << '\n';
  }

  inline std::string read_header(std::istream &istream)
  {
    // check magic bytes an version number
    version_t version = read_magic(istream);

    uint32_t header_length = 0;
    if (version == version_t{1, 0})
    {
      std::array<uint8_t, 2> header_len_le16{};
      istream.read(reinterpret_cast<char *>(header_len_le16.data()), 2);
      header_length = (header_len_le16[0] << 0) | (header_len_le16[1] << 8);

      if ((magic_string_length + 2 + 2 + header_length) % 16 != 0)
      {
        // TODO(llohse): display warning
      }
    }
    else if (version == version_t{2, 0})
    {
      std::array<uint8_t, 4> header_len_le32{};
      istream.read(reinterpret_cast<char *>(header_len_le32.data()), 4);

      header_length =
          (header_len_le32[0] << 0) | (header_len_le32[1] << 8) | (header_len_le32[2] << 16) | (header_len_le32[3] << 24);

      if ((magic_string_length + 2 + 4 + header_length) % 16 != 0)
      {
        // TODO(llohse): display warning
      }
    }
    else
    {
      throw std::runtime_error("unsupported file format version");
    }

    auto buf_v = std::vector<char>(header_length);
    istream.read(buf_v.data(), header_length);
    std::string header(buf_v.data(), header_length);

    return header;
  }

  template <typename T>
  struct npy_data
  {
    std::vector<T> data = {};
    shape_t shape = {};
    bool fortran_order = false;
    stype_t orig_data = stype_t_helper<T>::val;
  };

  template <typename T>
  struct npy_data_ptr
  {
    const T *data_ptr = nullptr;
    shape_t shape = {};
    bool fortran_order = false;
  };

  inline stype_t probe_npy_scalartype(std::istream &in)
  {
    std::string header_s = read_header(in);
    // parse header
    header_t header = parse_header(header_s);
    return header.dtype.stype;
  }

  inline stype_t probe_npy_scalartype(const std::string &filename)
  {
    auto stream = filename_to_ifstream(filename, "npy");
    return probe_npy_scalartype(stream);
  }

  template <typename S, typename T>
  inline void read_bytes_into_vec(std::istream &in, shape_t shape, std::vector<T>& target, bool transpose = false)
  {
    // create a buffer vector
    const auto size = shape.size();
    const auto dims = shape.ndims();
    auto buf = std::vector<S>();
    buf.resize(size);
    in.read(reinterpret_cast<char *>(buf.data()), sizeof(S) * size);
    // manual conversion (and optional transpose)
    if(transpose){
      auto iter = std::vector<ndarray_len_t>(dims, 0);
      auto i = 0;
      while(true){
        // calculate index
        std::size_t index = 0;
        std::size_t mult = 1;
        for(auto axis = 0; axis < dims; ++axis){
          index += iter[axis] * mult;
          mult *= shape[axis];
        }
        target[i++] = static_cast<T>(buf[index]);
        // increment (from right) or break out
        auto axis = 0;
        while(axis < dims){
          const auto invaxis = dims - axis - 1;
          iter[invaxis] += 1;
          if(iter[invaxis] == shape[invaxis]){
            iter[invaxis] = 0;
            ++axis;
          }else{
            break;
          }
        }
        if(axis >= dims){
          break;
        }
      }
    }else{
      for(auto i = 0; i < size; ++i){
        target[i] = static_cast<T>(buf[i]);
      }
    }
  }

  template <typename T>
  inline void read_bytes_into_vec_dyn(std::istream &in, stype_t stype, shape_t shape, std::vector<T>& target, bool transpose = false)
  {
    target.resize(shape.size());
    // read the data
    if(stype == stype_t_helper<T>::val && !transpose){
      // data types from source and target are the same, we only need to copy the bytes (if we do not need to transpose)
      in.read(reinterpret_cast<char *>(target.data()), sizeof(T) * shape.size());  
    }else{
      // warn if its not an upcast
      if (!upcast(stype, stype_t_helper<T>::val)){
        std::cout << "WARNING: input image was converted to type with lower precision" << std::endl;
      }
      // data types are different, use a buffer (dependent on src type)
      switch (stype)
      {
      case stype_t::int8:
        return read_bytes_into_vec<int8_t,T>(in, shape, target, transpose);
      case stype_t::int16:
        return read_bytes_into_vec<int16_t,T>(in, shape, target, transpose);
      case stype_t::int32:
        return read_bytes_into_vec<int32_t,T>(in, shape, target, transpose);
      case stype_t::int64:
        return read_bytes_into_vec<int64_t,T>(in, shape, target, transpose);
      case stype_t::uint8:
        return read_bytes_into_vec<uint8_t,T>(in, shape, target, transpose);
      case stype_t::uint16:
        return read_bytes_into_vec<uint16_t,T>(in, shape, target, transpose);
      case stype_t::uint32:
        return read_bytes_into_vec<uint32_t,T>(in, shape, target, transpose);
      case stype_t::uint64:
        return read_bytes_into_vec<uint64_t,T>(in, shape, target, transpose);
      case stype_t::f32:
        return read_bytes_into_vec<float,T>(in, shape, target, transpose);
      case stype_t::f64:
        return read_bytes_into_vec<double,T>(in, shape, target, transpose);
      default:
        throw std::invalid_argument("Given type is not valid");
      }
    }
  }

  /**
   * Always return a c-contiguous array
   */
  template <typename T>
  inline npy_data<T> read_npy(std::istream &in)
  {
    // check deconstructor of header_s
    std::string header_s = read_header(in);
    // parse header
    header_t header = parse_header(header_s);

    npy_data<T> data;

    data.shape = header.shape;
    data.fortran_order = false;
    data.orig_data = header.dtype.stype;

    // read the data
    read_bytes_into_vec_dyn(in, header.dtype.stype, data.shape, data.data, header.fortran_order);
    return data;
  }

  template <typename T>
  inline npy_data<T> read_npy(const std::string &filename)
  {
    auto stream = filename_to_ifstream(filename, "npy");
    return read_npy<T>(stream);
  }

  template <typename T>
  inline void write_npy(std::ostream &out, const npy_data<T> &data)
  {
    //  static_assert(has_typestring<T>::value, "T type not
    //  understood");
    const dtype_t dtype = dtype_t(stype_t_helper<T>::val);

    header_t header{dtype, data.fortran_order, data.shape};
    write_header(out, header);

    auto size = static_cast<size_t>(data.shape.size());

    out.write(reinterpret_cast<const char *>(data.data.data()), sizeof(T) * size);
  }

  template <typename T>
  inline void write_npy(const std::string &filename, const npy_data<T> &data)
  {
    auto stream = filename_to_ofstream(filename, "npy");
    write_npy<T>(stream, data);
  }

  template <typename T>
  inline void write_npy(std::ostream &out, const npy_data_ptr<T> &data_ptr)
  {
    // const dtype_t dtype = dtype_map.at(std::type_index(typeid(T)));
    const dtype_t dtype = dtype_t(stype_t_helper<T>::val);

    header_t header{dtype, data_ptr.fortran_order, data_ptr.shape};
    write_header(out, header);

    auto size = static_cast<size_t>(data_ptr.shape.size());

    out.write(reinterpret_cast<const char *>(data_ptr.data_ptr), sizeof(T) * size);
  }

  template <typename T>
  inline void write_npy(const std::string &filename, const npy_data_ptr<T> &data_ptr)
  {
    auto stream = filename_to_ofstream(filename, "npy");
    write_npy<T>(stream, data_ptr);
  }
} // namespace npy

#endif // NPY_HPP_
