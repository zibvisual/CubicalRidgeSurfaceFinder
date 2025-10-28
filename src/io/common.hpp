#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <filesystem>
#include <ios>

inline std::ifstream filename_to_ifstream(const std::string &filename, const std::string &extension, const bool replace_extension = false){
    // allow relative image paths
    auto input = std::filesystem::current_path() / filename;
    if (!input.has_extension() || replace_extension)
    {
        input = input.replace_extension(extension);
    }
    std::ifstream stream(input, std::ifstream::binary);
    if (!stream)
    {
        throw std::ios_base::failure("io error: failed to open file " + input.string() + ".");
    }
    return stream;
}

inline std::ofstream filename_to_ofstream(const std::string &filename, const std::string &extension, const bool replace_extension = false){
    // allow relative image paths
    auto output = std::filesystem::current_path() / filename;
    if (!output.has_extension() || replace_extension)
    {
        output = output.replace_extension(extension);
    }
    std::ofstream stream(output, std::ofstream::binary);
    if (!stream)
    {
        throw std::ios_base::failure("io error: failed to open file " + output.string() + ".");
    }
    return stream;
}