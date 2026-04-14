#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <filesystem>
#include <ios>
#include <bit>
#include <algorithm>

inline std::filesystem::path add_suffix(const std::filesystem::path& path, const std::string& suffix = ""){
    auto filename = path.stem();
    filename += suffix;
    filename += path.extension();
    // allow relative image paths from the working directory
    auto result = std::filesystem::current_path() / path.parent_path() / filename;
    return result;
}

inline std::ifstream path_to_ifstream(std::filesystem::path& input, const std::string& default_extension, const bool enforce_extenstion = false){
    // if no extension is set, use the default one
    if(!input.has_extension() || enforce_extenstion){
        input.replace_extension(default_extension);
    }
    std::ifstream stream(input, std::ifstream::binary);
    if (!stream)
    {
        throw std::ios_base::failure("io error: failed to open file \"" + input.string() + "\"");
    }
    return stream;
}

inline std::ofstream path_to_ofstream(std::filesystem::path& output, const std::string& extension){
    // enforce the correct extension
    output.replace_extension(extension);
    // create necessary directories
    if(!output.parent_path().empty()){
        std::filesystem::create_directories(output.parent_path());
    }
    std::ofstream stream(output, std::ofstream::binary);
    if (!stream)
    {
        throw std::ios_base::failure("io error: failed to open file \"" + output.string() + "\"");
    }else{
        // std::cout << "Write into file " << output << std::endl;
    }
    return stream;
}

// Scalar types enum
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

inline std::string to_c_string(const stype_t& stype)
{
    switch (stype)
    {
    case stype_t::int8:
        return "signed char";
    case stype_t::int16:
        return "short";
    case stype_t::int32:
        return "int";
    case stype_t::int64:
        return "long long";
    case stype_t::uint8:
        return "unsigned char";
    case stype_t::uint16:
        return "unsigned short";
    case stype_t::uint32:
        return "unsigned int";
    case stype_t::uint64:
        return "unsigned long long";
    case stype_t::f32:
        return "float";
    case stype_t::f64:
        return "double";
    default:
        return "<invalid type>";
    }
}

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
        return to == stype_t::int8 || to == stype_t::int16 || to == stype_t::int32 || to == stype_t::int64 || to == stype_t::f32 || to == stype_t::f64;
    case stype_t::int16:
        return to == stype_t::int16 || to == stype_t::int32 || to == stype_t::int64 || to == stype_t::f32 || to == stype_t::f64;
    case stype_t::int32:
        return to == stype_t::int32 || to == stype_t::int64 || to == stype_t::f64;
    case stype_t::int64:
        return to == stype_t::int64;
    case stype_t::uint8:
        return to == stype_t::int16 || to == stype_t::int32 || to == stype_t::int64 || to == stype_t::uint8 || to == stype_t::uint16 || to == stype_t::uint32 || to == stype_t::uint64 || to == stype_t::f32 || to == stype_t::f64;
    case stype_t::uint16:
        return to == stype_t::int32 || to == stype_t::int64 || to == stype_t::uint16 || to == stype_t::uint32 || to == stype_t::uint64 || to == stype_t::f32 || to == stype_t::f64;
    case stype_t::uint32:
        return to == stype_t::int64 || to == stype_t::uint32 || to == stype_t::uint64 || to == stype_t::f64;
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

// template <typename T>
// void
// swap_endian(T* data, std::size_t size)
// {
//     uint8_t buf[sizeof(T)];
//     for (std::size_t i = 0; i < size; ++i)
//     {
//         auto* ptr_start = reinterpret_cast<uint8_t*>(&data[i]);
//         std::reverse_copy(ptr_start, ptr_start + sizeof(S), &buf[0]);
//         auto* ptr = reinterpret_cast<T*>(&buf[0]);
//         data[i] = *ptr;
//     }
// }

template <typename S, typename T>
void
convert(S* from, T* to, std::size_t size)
{
    for (std::size_t i = 0; i < size; ++i)
    {
        to[i] = static_cast<T>(from[i]);
    }
}

/**
 * Swaps endian of S
 */
template <typename S, typename T>
void
convert_and_swap_endian(S* from, T* to, std::size_t size)
{
    uint8_t buf[sizeof(S)];
    for (std::size_t i = 0; i < size; ++i)
    {
        auto* ptr_start = reinterpret_cast<uint8_t*>(&from[i]);
        // TODO: we could change from by swapping bytes instead of copying
        std::reverse_copy(ptr_start, ptr_start + sizeof(S), &buf[0]);
        auto* ptr = reinterpret_cast<S*>(&buf[0]);
        to[i] = static_cast<T>(*ptr);
    }
}

template <class S, class T>
void
read(std::ifstream& inputFile, T* data, std::size_t size, std::endian endian)
{
    if(stype_t_helper<S>::val == stype_t_helper<T>::val && std::endian::native == endian){
        // same endianess and same value -> just copy the bytes
        inputFile.read(reinterpret_cast<char*>(data), size * sizeof(T));
    }
    else{
        // we need to buffer input
        S buf[1024];
        auto counter = 0;
        while(counter < size){
            auto len = std::min(size - counter, static_cast<std::size_t>(1024));
            inputFile.read(reinterpret_cast<char*>(buf), len * sizeof(S));
            std::endian::native == endian ? convert(buf, data, len) : convert_and_swap_endian(buf, data, len);
        }
    }
}



// template <typename T>
// inline void read_bytes_into_vec_dyn(std::istream &in, stype_t stype, Dims dims, std::vector<T>& target)
// {
//     target.resize(dims.size());
//     // read the data
//     if(stype == stype_t_helper<T>::val){
//         // data types from source and target are the same, we only need to copy the bytes
//         in.read(reinterpret_cast<char *>(target.data()), sizeof(T) * dims.size()); // TODO: swap byte order if necessary --> see nrrd
//     }else{
//         // warn if its not an upcast
//         if (!upcast(stype, stype_t_helper<T>::val)){
//             std::cout << "WARNING: input image was converted to type with lower precision (" << stype << " to " << stype_t_helper<T>::val << ")" << std::endl;
//         }
//         // data types are different, use a buffer (dependent on src type)
//         // for now just error
//         std::cerr << "Error: Data type could no be converted" << std::endl;
//         // switch (stype)
//         // {
//         // case stype_t::int8:
//         //     return read_bytes_into_vec<int8_t,T>(in, dims, target);
//         // case stype_t::int16:
//         //     return read_bytes_into_vec<int16_t,T>(in, dims, target);
//         // case stype_t::int32:
//         //     return read_bytes_into_vec<int32_t,T>(in, dims, target);
//         // case stype_t::int64:
//         //     return read_bytes_into_vec<int64_t,T>(in, dims, target);
//         // case stype_t::uint8:
//         //     return read_bytes_into_vec<uint8_t,T>(in, dims, target);
//         // case stype_t::uint16:
//         //     return read_bytes_into_vec<uint16_t,T>(in, dims, target);
//         // case stype_t::uint32:
//         //     return read_bytes_into_vec<uint32_t,T>(in, dims, target);
//         // case stype_t::uint64:
//         //     return read_bytes_into_vec<uint64_t,T>(in, dims, target);
//         // case stype_t::f32:
//         //     return read_bytes_into_vec<float,T>(in, dims, target);
//         // case stype_t::f64:
//         //     return read_bytes_into_vec<double,T>(in, dims, target);
//         // default:
//         //     throw std::invalid_argument("Given type is not valid");
//         // }
//     }
// }

template <typename T>
void
parseToken(char* token, T& val)
{
    // should never be called
    val = static_cast<T>(0);
}