#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>

#include "Dims.hpp"
#include "UnionFind.hpp"

/**
 * Generated sample points to use for seed points in the Ridge Surface Finder. Threshold must be positive.
 */
std::vector<uint32_t> sample(const float* data, Dims dims, float threshold = 0.0f){
    auto uf = UnionFind(dims.size());

    // generate the components
    size_t counter = 0;
    for(size_t z = 0; z < dims.z(); ++z){
        for(size_t y = 0; y < dims.y(); ++y){
            for(size_t x = 0; x < dims.x(); ++x){
                if(data[counter] <= threshold){
                    ++counter;
                    continue;
                }
                auto right = dims.c_index(x+1, y, z);
                auto up = dims.c_index(x, y+1, z);
                auto prev = dims.c_index(x, y, z+1);
                if(right.has_value() && data[right.value()] > threshold){
                    uf.join(counter, right.value());
                }
                if(up.has_value() && data[up.value()] > threshold){
                    uf.join(counter, up.value());
                }
                if(prev.has_value() && data[prev.value()] > threshold){
                    uf.join(counter, prev.value());
                }
                ++counter;
            }
        }
    }

    // track maxima
    std::unordered_map<uint32_t, std::pair<float, uint32_t>> maxima;
    for(size_t i = 0; i < dims.size(); ++i){
        if(data[i] > threshold){
            auto id = uf.find(i);
            if(data[i] > maxima[id].first){
                maxima[id] = std::pair(data[i], i);
            }
        }
    }

    // extract to result
    std::vector<uint32_t> result;
    result.reserve(maxima.size());
    for(auto pair : maxima){
        result.push_back(pair.second.second);
    }
    return result;
}