

 #pragma once

#include <vector>
#include <array>
#include <utils/Vec.hpp>
#include "common.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

/**
 * Class to read and write trivial graph file format (with included spatial information)
 */
struct compact_graph_data_t {
   std::vector<VecFloat> points;
   std::vector<std::array<std::size_t, 2>> edges;
};

struct graph_data_t {
    std::unordered_map<std::string, VecFloat> points;
    std::vector<std::array<std::string, 2>> edges;
};

/**
 * Reads the graph data.
 */
inline graph_data_t read_graph_data(std::istream& input){
    std::unordered_map<std::string, VecFloat> points;
    std::vector<std::array<std::string, 2>> edges;

   std::string line;
   while(std::getline(input, line)){
      std::stringstream line_stream(line);
      std::string token;
      // ignore empty lines
      if(!(line_stream >> token)){
         continue;
      }
      // ignore comments
      if(token.starts_with("#")){
         continue;
      }

      // check if the token exists in the mapping --> we assume an edge
      if (points.contains(token)){
        // edge
        std::string to;
        if(!(line_stream >> to) || !points.contains(to)){
            // ERROR: no valid second id for edge given
            continue;
        }
        // TODO: read edge meta values
        edges.push_back({token, to});
      }else{
        // vertex
        float x, y, z;
        if(!(line_stream >> x >> y >> z)){
           // ERROR: 3 float points are MANDATORY
           continue;
        }
        // TODO: read vertex meta values
        points[token] = VecFloat(x,y,z);
      }
   }
   return graph_data_t{points, edges};
}

inline graph_data_t read_graph_data(const std::string &filename)
{
  auto stream = filename_to_ifstream(filename, "tgf");
  return read_graph_data(stream);
}

/**
 * Reads the graph data but does not care about the ids of the vertices or edges.
 */
inline compact_graph_data_t read_graph_data_rename(std::istream& input){
    std::vector<VecFloat> points;
    std::vector<std::array<std::size_t, 2>> edges;
    std::unordered_map<std::string, std::size_t> mapping;

   std::string line;
   while(std::getline(input, line)){
      std::stringstream line_stream(line);
      std::string token;
      // ignore empty lines
      if(!(line_stream >> token)){
         continue;
      }
      // ignore comments
      if(token.starts_with("#")){
         continue;
      }

      // check if the token exists in the mapping --> we assume an edge
      if (mapping.contains(token)){
        // edge
        auto from = mapping[token];
        if(!(line_stream >> token) || !mapping.contains(token)){
            // ERROR: no valid second id for edge given
            continue;
        }
        auto to = mapping[token];
        // TODO: read edge values
        edges.push_back({from, to});
      }else{
        // vertex
        float x, y, z;
        if(!(line_stream >> x >> y >> z)){
           // ERROR: 3 float points are MANDATORY
           continue;
        }
        // TODO: read vertex values
        mapping[token] = points.size();
        points.push_back(VecFloat(x,y,z));
      }
   }
   return compact_graph_data_t{points, edges};
}

inline compact_graph_data_t read_graph_data_rename(const std::string &filename)
{
  auto stream = filename_to_ifstream(filename, "tgf");
  return read_graph_data_rename(stream);
}

inline void write_graph_data(std::ostream& output, graph_data_t data){
   // write all points first
   for(auto point : data.points){
      output << point.first << " " << point.second.x() << " " << point.second.y() << " " << point.second.z() << std::endl;  
   }

   // add a # line to make it valid tgf file format
   output << "#" << std::endl;

   // write edges
   for(auto edge : data.edges){
    output << edge[0] << " " << edge[1] << std::endl;
   }
}

inline void write_graph_data(const std::string &filename, graph_data_t data)
{
  auto stream = filename_to_ofstream(filename, "tgf");
  write_graph_data(stream, data);
}