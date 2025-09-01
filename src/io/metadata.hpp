/**
 * Class to read and write metafiles, which can be appended to e.g. numpy files to have more information.
 * This is used to include spatial information to numpy arrays.
 *
 * bbox: [0,1]
 * signal-bbox: true/false
 * /// OR
 * origin: [0.5]
 * voxelsize: [0.5]
 */

#pragma once

#include <vector>
#include <array>
#include <utils/Vec.hpp>
#include <utils/Lattice.hpp>
#include "common.hpp"
#include <iostream>
#include <sstream>
#include <string>

/**
 * Class to read and write simple metadata files.
 */
struct metadata_t
{
    std::vector<std::string> keywords;
    std::vector<std::vector<std::string>> data;
};

inline metadata_t read_metadata(std::istream &input)
{
    std::vector<std::string> keywords;
    std::vector<std::vector<std::string>> datas;

    std::string line;
    while (std::getline(input, line))
    {
        // use strtok as we allow some other symbols as well
        char* head = strtok(line.data(), " :=");
        if(head == nullptr || head[0] == '#'){
            // ignore empty lines or comments
            continue;
        }
        // first token is a keyword
        std::string keyword = head;
        keywords.push_back(keyword);

        std::vector<std::string> data;
        head = strtok(nullptr, " :=,;{}[]");
        while(head != nullptr){
            if(head[0] == '#'){
                break;
            }
            std::string datapoint = head;
            data.push_back(datapoint);
            head = strtok(nullptr, " ,;{}[]");
        }
        datas.push_back(data);
    }
    return metadata_t{keywords, datas};
}

inline metadata_t read_metadata(const std::string &filename)
{
    auto stream = filename_to_ifstream(filename, "txt");
    return read_metadata(stream);
}

inline void write_metadata(std::ostream &output, metadata_t data)
{
    std::size_t size = std::min(data.keywords.size(), data.data.size());
    for(std::size_t i = 0; i < size; ++i){
        output << data.keywords[i] << ": ";
        if(data.data[i].size() == 1){
            output << data.data[i][0] << std::endl;
        }else if (data.data[i].size() > 1){
            output << "[ ";
            for(std::size_t j = 0; j < data.data[i].size(); ++j){
                output << data.data[i][j] << " ";
            }
            output << "]" << std::endl;
        }
    }
}

inline void write_metadata(const std::string &filename, metadata_t data)
{
    auto stream = filename_to_ofstream(filename, "txt");
    write_metadata(stream, data);
}

inline std::ostream& operator<<(std::ostream& os, const metadata_t& data)
{
    std::size_t size = std::min(data.keywords.size(), data.data.size());
    for(std::size_t i = 0; i < size; ++i){
        os << data.keywords[i] << ": ";
        if(data.data[i].size() == 1){
            os << data.data[i][0] << std::endl;
        }else if (data.data[i].size() > 1){
            os << "[ ";
            for(std::size_t j = 0; j < data.data[i].size(); ++j){
                os << data.data[i][j] << " ";
            }
            os << "]" << std::endl;
        }
    }
    return os;
}