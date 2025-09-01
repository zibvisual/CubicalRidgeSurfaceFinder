#pragma once

#include <vector>
#include <array>
#include <utils/Vec.hpp>
#include "common.hpp"
#include <iostream>
#include <sstream>
#include <string>

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
 * 
 * IMPORTANT: triangles and patch indices start with 0 and not 1 as in wavefront file.
 */
struct wavefront_data_t {
   std::vector<VecFloat> points;
   std::vector<std::array<std::size_t, 3>> triangles;
   std::vector<std::size_t> patches;
};

inline wavefront_data_t read_wavefront(std::istream& input){
   std::vector<VecFloat> points;
   std::vector<std::array<std::size_t, 3>> triangles;
   std::vector<std::size_t> patches;

   std::string line;
   std::size_t last_patch_triangle = 0;
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
      }else if (token == "g" || token == "o"){
         // create a new patch (ignore empty patches)
         if(last_patch_triangle == triangles.size()){
            continue;
         }
         patches.push_back(triangles.size());
         last_patch_triangle = triangles.size();
      }
   }
   return wavefront_data_t{points, triangles, patches};
}

inline wavefront_data_t read_wavefront(const std::string &filename)
{
  auto stream = filename_to_ifstream(filename, "obj");
  return read_wavefront(stream);
}

inline void write_wavefront(std::ostream& output, wavefront_data_t data){
   // write all points first
   for(auto point : data.points){
      output << "v " << point.x() << " " << point.y() << " " << point.z() << std::endl;  
   }
   // write triangles (with patches interspersed)
   std::size_t patch_counter = 0;

   output << "g " << (patch_counter + 1) << std::endl;
   for(std::size_t i = 0; i < data.triangles.size(); ++i){
      if(patch_counter < data.patches.size() && data.patches[patch_counter] == i){
         ++patch_counter;
         output << "g " << (patch_counter + 1) << std::endl;
      }
      output << "f " << (data.triangles[i][0]+1) << " " << (data.triangles[i][1]+1) << " " << (data.triangles[i][2]+1) << std::endl;
   }
}

inline void write_wavefront(const std::string &filename, wavefront_data_t data)
{
  auto stream = filename_to_ofstream(filename, "obj");
  write_wavefront(stream, data);
}