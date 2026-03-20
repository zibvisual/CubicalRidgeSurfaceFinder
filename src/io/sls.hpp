

#pragma once

#include <vector>
#include <array>
#include <utils/Vec.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "common.hpp"
#include <graph/LineSet.hpp>

// /**
//  * Reads the graph data.
//  */
// inline graph_data_t read_graph_data(std::istream& input){
//     std::unordered_map<std::string, VecFloat> points;
//     std::vector<std::array<std::string, 2>> edges;

//    std::string line;
//    while(std::getline(input, line)){
//       std::stringstream line_stream(line);
//       std::string token;
//       // ignore empty lines
//       if(!(line_stream >> token)){
//          continue;
//       }
//       // ignore comments
//       if(token.starts_with("#")){
//          continue;
//       }

//       // check if the token exists in the mapping --> we assume an edge
//       if (points.contains(token)){
//         // edge
//         std::string to;
//         if(!(line_stream >> to) || !points.contains(to)){
//             // ERROR: no valid second id for edge given
//             continue;
//         }
//         // TODO: read edge meta values
//         edges.push_back({token, to});
//       }else{
//         // vertex
//         float x, y, z;
//         if(!(line_stream >> x >> y >> z)){
//            // ERROR: 3 float points are MANDATORY
//            continue;
//         }
//         // TODO: read vertex meta values
//         points[token] = VecFloat(x,y,z);
//       }
//    }
//    return graph_data_t{points, edges};
// }

// inline graph_data_t read_graph_data(const std::string &filename)
// {
//   auto stream = filename_to_ifstream(filename, "sls");
//   return read_graph_data(stream);
// }

/**
 * Write seeding line set file format (see https://www.csc.kth.se/~weinkauf/notes/seedinglinesets.html).
 * 
 * The template parameter should be a std::tuple!
 */
template<typename T>
inline void write_line_data(std::ostream& output, const LineSet<T>& lineset){
    // header
    output << "SeedingLS" << std::endl << "Version 1" << std::endl << std::endl;
    output << "NumOfSeedingLines " << lineset.size() << std::endl;
    output << "NumOfDataPerPoint " << lineset.data_size() << std::endl << std::endl;

    // data
    for(const auto& line : lineset.lines()){
        output << "LinePoints " << line.size() << std::endl;
        for(const auto& point : line){
            output << point.point.x() << " " << point.point.y() << " " << point.point.z() << std::endl;
            // std::apply([&output](auto&&... args){((output << args << std::endl), ...);}, point.data);
            point.print_data(output);
        }
        output << std::endl;
    }
}

template<typename T>
inline void write_line_data(const std::string &filename, const LineSet<T>& data)
{
  auto stream = filename_to_ofstream(filename, "sls");
  write_line_data(stream, data);
}