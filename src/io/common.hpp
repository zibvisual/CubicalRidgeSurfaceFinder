#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <filesystem>
#include <ios>

inline std::filesystem::path filename_to_path(const std::filesystem::path& filename, const std::string& suffix = ""){
    auto comp = filename.stem();
    comp += suffix;
    // allow relative image paths
    return std::filesystem::current_path() / comp;
}

inline void path_extension(std::filesystem::path& path, const std::string& extension, const bool replace_extension = false){
    if(!path.has_extension() || replace_extension){
        path.replace_extension(extension);
    }else{
        path += "." + extension;
    }
}

inline std::ifstream path_to_ifstream(std::filesystem::path& input, const std::string& extension, const bool replace_extension = false){
    path_extension(input, extension, replace_extension);
    std::ifstream stream(input, std::ifstream::binary);
    if (!stream)
    {
        throw std::ios_base::failure("io error: failed to open file " + input.string() + ".");
    }
    return stream;
}

inline std::ofstream path_to_ofstream(std::filesystem::path& output, const std::string& extension, const bool replace_extension = false){
    path_extension(output, extension, replace_extension);
    std::filesystem::create_directories(output.parent_path());
    std::ofstream stream(output, std::ofstream::binary);
    if (!stream)
    {
        throw std::ios_base::failure("io error: failed to open file " + output.string() + ".");
    }else{
        std::cout << "Write into file " << output << std::endl;
    }
    return stream;
}