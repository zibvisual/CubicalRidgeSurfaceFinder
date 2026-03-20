#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <span>

#include "common.hpp"
#include <utils/Vec.hpp>
#include <surface/Triangle.hpp>

// wavefront:
// # point:
// v float float float
// # normal:
// vn float float float
// # face
// f index index index
// # patches
// g groupname/groupnumber

/** 
 * Class to read and write simple wavefront object files. Not all features are supported.
 * We especially only allow for triangles and no other polygon.
 * 
 * IMPORTANT: triangles and patch indices start with 0 and not 1 as in wavefront file.
 */
struct wavefront_data_t {
   std::vector<VecFloat> points;
   std::vector<std::array<std::size_t, 3>> triangles;
};

inline wavefront_data_t read_wavefront(std::istream& input){
   std::vector<VecFloat> points;
   std::vector<std::array<std::size_t, 3>> triangles;

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

      if(token == "v") {
         float x, y, z;
         if(!(line_stream >> x >> y >> z)){
            // ERROR: 3 float points are MANDATORY
         }
         points.push_back(VecFloat(x,y,z));
      }else if (token == "f") {
         std::size_t p1, p2, p3;
         if(!(line_stream >> p1 >> p2 >> p3)){
            // ERROR: 3 indices (for triangles) are MANDATORY
         }
         // TODO: check if we have more than just triangles... (more than 3 indices)
         triangles.push_back({p1-1,p2-1,p3-1});
      }
   }
   return wavefront_data_t{points, triangles};
}

inline wavefront_data_t read_wavefront(std::filesystem::path input)
{
  auto stream = path_to_ifstream(input, "obj");
  return read_wavefront(stream);
}

inline void write_wavefront(std::ostream& output, wavefront_data_t data){
   // write all points first
   for(const auto& point : data.points){
      output << "v " << point.x() << " " << point.y() << " " << point.z() << std::endl;  
   }
   
   // write triangles
   for(const auto& triangle: data.triangles){
      output << "f " << (triangle[0]+1) << " " << (triangle[1]+1) << " " << (triangle[2]+1) << std::endl;
   }
}

inline void write_wavefront(std::filesystem::path output, wavefront_data_t data)
{
   auto stream = path_to_ofstream(output, "obj");
   write_wavefront(stream, data);
}

inline void write_wavefront(std::ostream& output, std::span<const VecFloat> points, std::span<const surface::SimpleTriangle> triangles)
{
   // write all points first
   for(auto point : points){
      output << "v " << point.x() << " " << point.y() << " " << point.z() << std::endl;  
   }
   
   // write triangles
   for(std::size_t i = 0; i < triangles.size(); ++i){
      output << "f " << (triangles[i][0]+1) << " " << (triangles[i][1]+1) << " " << (triangles[i][2]+1) << std::endl;
   }
}

inline void write_wavefront(std::filesystem::path output, std::span<const VecFloat> points, std::span<const surface::SimpleTriangle> triangles)
{
   auto stream = path_to_ofstream(output, "obj");
   write_wavefront(stream, points, triangles);
}