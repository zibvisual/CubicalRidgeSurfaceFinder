#pragma once

#include <vector>
#include <array>
#include <utils/Vec.hpp>
#include "common.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <exception>

/** 
 * Class to read and write simple vertexsets. Not all features are supported.
 */
template <class T>
struct vertex_data {
   std::vector<VecFloat> points;
   std::vector<T> data;
};

// template <class T>
// Just do simple float datapoint for now. Later on, one can use tuple structs (and constexpr) to allow for generic reads
inline vertex_data<float> read_vertexset(std::istream& input, float defaults){
   std::vector<VecFloat> points;
   std::vector<float> data;

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

      float x, y, z;
      if(token == "v"){
        if(!(line_stream >> x >> y >> z)){
            // ERROR: 3 float points are MANDATORY
        }
      }else{
        try{
            x = std::stof(token);
        } catch(const std::exception &err){
            // ERROR: neither v nor a float
        }
        if(!(line_stream >> y >> z)){
            // ERROR: 3 float points are MANDATORY
        }
      }
      points.push_back(VecFloat(x,y,z));

      // get data if there is data
      float d;
      if(!(line_stream >> d)){
        d = defaults;
      }
      data.push_back(d);
   }
   return vertex_data<float>{points, data};
}

inline vertex_data<float> read_vertexset(const std::string &filename, float defaults)
{
  auto stream = filename_to_ifstream(filename, "obj");
  return read_vertexset(stream, defaults);
}

inline void write_vertexset(std::ostream& output, vertex_data<float> data){
    // write all points
    for(std::size_t i = 0; i < data.points.size(); ++i){
        auto point = data.points[i];
        output << "v " << point.x() << " " << point.y() << " " << point.z();
        output << " " << data.data[i] << std::endl;
    }
}

inline void write_vertexset(const std::string &filename, vertex_data<float> data)
{
  auto stream = filename_to_ofstream(filename, "obj");
  write_vertexset(stream, data);
}