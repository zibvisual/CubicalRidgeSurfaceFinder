#pragma once

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdint>
#include <optional>
#include <bit>

#include <utils/Vec.hpp>
#include <utils/BBox.hpp>
#include <utils/Dims.hpp>

namespace nrrd {

template <typename T>
struct nrrd_data
{
    Lattice lattice = Lattice();
    std::vector<T> data = {};
};

template <class T>
int fillScalarField(std::ifstream& inputFile, T* data, std::size_t size);

inline bool is_nrrd_file(std::ifstream& inputFile)
{
    std::string magic_bytes = "NRRD00";
    std::string line;

    return std::getline(inputFile, line) && line.starts_with(magic_bytes);
}

template <typename T>
std::optional<nrrd_data<T>>
read_nrrd(std::ifstream& inputFile, SpatialArguments command_args = SpatialArguments())
{
    if(!is_nrrd_file(inputFile)){
        return {};
    }

    const auto supported_encodings = std::array<std::string, 4>({"raw", "txt", "text", "ascii"});

    const std::string propertyDimension = "dimension: ";
    const std::string propertyType = "type: ";
    const std::string propertySizes = "sizes: ";
    const std::string propertySpacings = "spacings: ";
    const std::string propertyEncoding = "encoding: ";
    const std::string propertyMin = "axis mins: ";
    const std::string propertyMax = "axis maxs: ";
    const std::string propertyEndian = "endian: ";

    std::string line;
    std::string type;
    std::string str_encoding;
    std::endian endian = std::endian::native;
    Dims sizes;
    std::optional<VecFloat> spacings;
    std::optional<VecFloat> spaceMins;
    std::optional<VecFloat> spaceMaxs;
    std::optional<bool> centered;

    int encoding = 0;

    while (std::getline(inputFile, line))
    {
        std::stringstream line_stream(line);
        std::string property;
        // header is finished if we reached an empty line
        if(!(line_stream >> property)){
            break;
        }

        if (property == "dimension:")
        {
            int dim;
            if(!(line_stream >> dim)){
                throw std::invalid_argument("Property \"dimension\" must contain a value!");
            }
            if(dim != 3){
                throw std::invalid_argument("Only 3 dimensional images are supported!");
            }
        }
        else if (property == "type:")
        {
            // get the rest of the line
            type = line.substr(propertyType.length());
        }
        else if (property == "sizes:")
        {
            std::size_t x, y, z;
            if(!(line_stream >> x >> y >> z)){
                throw std::invalid_argument("Property \"sizes\" must have 3 integers");
            }
            sizes = Dims(x,y,z);
        }
        else if (property == "spacings:")
        {
            float x, y, z;
            if(!(line_stream >> x >> y >> z)){
                throw std::invalid_argument("Property \"spacings\" must have 3 floats");
            }
            spacings = VecFloat(x,y,z);
        }
        else if (property == "encoding:")
        {
            std::string str_encoding = line.substr(propertyEncoding.length());
            if (std::find(supported_encodings.begin(), supported_encodings.end(), str_encoding) == supported_encodings.end())
            {
                throw std::invalid_argument("Only raw and ascii encoding is supported! Property \"encoding\" was " + str_encoding);
                // std::cerr << "Error: Only raw and ascii encoding is supported!" << std::endl;
                // return {};
            }
            encoding = (str_encoding == "raw");
        }
        else if (line.starts_with(propertyMax) || property == "axismaxs:")
        {
            if(!property.ends_with(":")){
                // consume another token
                line_stream >> property;
            }
            float x, y, z;
            if(!(line_stream >> x >> y >> z)){
                throw std::invalid_argument("Property \"axis maxs\" must have 3 floats");
            }
            spaceMaxs = VecFloat(x,y,z);
        }
        else if (line.starts_with(propertyMin) || property == "axismins:")
        {
            if(!property.ends_with(":")){
                // consume another token
                line_stream >> property;
            }
            float x, y, z;
            if(!(line_stream >> x >> y >> z)){
                throw std::invalid_argument("Property \"axis mins\" must have 3 floats");
            }
            spaceMins = VecFloat(x,y,z);
        }
        else if (property == "endian:")
        {
            std::string endian_str;
            if(!(line_stream >> endian_str)){
                throw std::invalid_argument("Property \"endian\" must contain a value!");
            }
            if(endian_str == "big"){
                endian = std::endian::big;
            }else if (endian_str == "little"){
                endian = std::endian::little;
            }else{
                throw std::invalid_argument("Property \"endian\" can only be either \"big\" or \"little\", but we got \"!" + endian_str + "\"");
            }
        }
        else if (property == "centers:" || property == "centerings:"){
            std::string c1, c2, c3;
            if(!(line_stream >> c1 >> c2 >> c3)){
                throw std::invalid_argument("Property \""+ property +"\" must contain 3 values!");
            }
            if(c1 != c2 || c2 != c3 || c1 != c3){
                throw std::invalid_argument("Image centers of all axis are not the same. We only support either all centerings to be \"cell\" or \"node\"");
            }
            if(c1 == "cell"){
                centered = true;
            }else if(c1 == "node"){
                centered = false;
            }else{
                //if(c1 == "???" || c1 == "none")
                throw std::invalid_argument("All centerings must either be all \"cell\" or all \"node\"");
            }
        }
    }

    // check dims
    if (sizes.empty()){
        std::cerr << "Nrrd file must specifiy dimension of data" << std::endl;
        return {};
    }

    // create spatial args
    SpatialArguments args = SpatialArguments();
    if (spaceMins && spaceMaxs)
    {
        args.bbox = GenericBBox(spaceMins.value(), spaceMaxs.value());
    }
    args.origin = spaceMins;
    args.voxelsize = spacings;
    args.bbox_signal = centered;
    args.origin_signal = centered;

    // create lattice
    Lattice lattice = args.update(command_args).toInformation().toLattice(sizes);

    auto size = sizes.size();
    auto data = std::vector<T>(size);

    if(!encoding){
        // ascii
        fillScalarField(inputFile, data.data(), size);
        return nrrd_data{lattice, data};
    }

    // binary/raw
    if (type == "signed char" || type == "int8" || type == "int8_t")
    {
        read<int8_t, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "uchar" || type == "unsigned char" || type == "uint8" || type == "uint8_t")
    {
        read<uint8_t, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "short" || type == "short int" || type == "signed short" || type == "signed short int" || type == "int16" || type == "int16_t")
    {
        read<int16_t, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "ushort" || type == "unsigned short" || type == "unsigned short int" || type == "uint16" || type == "uint16_t")
    {
        read<uint16_t, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "int" || type == "signed int" || type == "int32" || type == "int32_t")
    {
        read<int32_t, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "uint" || type == "unsigned int" || type == "uint32" || type == "uint32_t")
    {
        read<uint32_t, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "longlong" || type == "long long" || type == "long long int" || type == "signed long long" || type == "signed long long int" || type == "int64" || type == "int64_t")
    {
        read<int64_t, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "ulonglong" || type == "unsigned long long" || type == "unsigned long long int" || type == "uint64" || type == "uint64_t")
    {
        read<uint64_t, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "float")
    {
        read<float, T>(inputFile, data.data(), size, endian);
    }
    else if (type == "double")
    {
        read<double, T>(inputFile, data.data(), size, endian);
    }
    else
    {
        std::cerr << "Datatype " << type << " is not supported" << std::endl;
        return {};
    }

    return nrrd_data{lattice, data};
}

template <class T>
int
fillScalarField(std::ifstream& inputFile, T* data, std::size_t size)
{

    // ascii
    std::size_t counter = 0;
    T val;
    std::string line;
    while (std::getline(inputFile, line))
    {
        char* str = const_cast<char*>(line.c_str());
        char* tokens = strtok(str, " ");
        while (tokens != NULL)
        {
            parseToken<T>(tokens, val);
            data[counter] = val;
            tokens = strtok(NULL, " ");
            counter++;
            if (counter >= size)
                break;
        }
        if (counter >= size)
            break;
    }
    
    return 1;
}

} // end namespace nrrd