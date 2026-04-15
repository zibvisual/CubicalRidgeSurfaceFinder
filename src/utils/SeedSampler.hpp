#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>

#include "Dims.hpp"
#include "UnionFind.hpp"

#include <fastmarching/FastMarching.hpp>
#include <fastmarching/FastMarchingBitSetObserver.hpp>

/**
 * Generated sample points to use for seed points in the Ridge Surface Finder. Threshold must be positive.
 */
std::vector<uint32_t> sample(const float* data, Dims dims, float threshold = 0.01f){
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

std::vector<uint32_t> greedyFmSampling(progressbar::Reporter& report, const RawFieldView<float>& data, float distance, float min, float max, std::size_t maxPoints = std::numeric_limits<std::size_t>::max())
{
    // sort by magnitude of the data
    std::priority_queue<std::pair<float, uint64_t> > queue;
    // create FM with a observer which sets the bitfield
    auto fm = fastmarching::FastMarching<fastmarching::ObserverBitSet<mutil::HashMapMappingView<BigHashMap<uint64_t, float>>>>(report);
    fm.observer().resize(data.dims().size());
    fm.data(&data);
    fm.thresholds(min, max);
    fm.maxDistance(distance);

    // push queue with all values
    for(std::size_t i = 0; i < data.dims().size(); ++i){
        float val = data.data()[i];
        if(val > min){
            queue.push({val, i});
        }
    }
    const float queue_max = static_cast<float>(queue.size());
    const float points_max = static_cast<float>(maxPoints);

    // greedy algo
    report.start("Searching for seed points");
    auto points = std::vector<uint32_t>();
    while (!queue.empty() && points.size() < maxPoints)
    {
        std::pair<float, uint64_t> entry = queue.top();
        queue.pop();
        auto index = entry.second;

        // near another seed point?
        if(fm.observer().bitset()[index]){
            continue;
        }
        
        // set point and run fm
        points.push_back(index);
        fm.setStartVoxel(index);
        fm.march();

        if(points.size() % 128 == 0){
            report.update(fmt::format("Queue: ({:d})", queue.size()));
            // report.update(std::max((queue_max - queue.size()) / queue_max, points.size() / points_max));
        }
    }
    report.end();

    return points;
}

// Two things to improve:
// - because of fm, most border points are not covered --> march a bit further with an inverse FM (which likes borders the most) OR with simple euclidean (a few more voxels in all directions)
// - we only need to greedily take the local maximas. When turning off voxels, we can check the border, what new paths are possible.
std::vector<uint32_t> greedySampling(progressbar::Reporter& report, const RawFieldView<float>& data, float distance, float min, float max, std::size_t maxPoints = std::numeric_limits<std::size_t>::max())
{
    const auto dims = data.getDims();

    // create FM with a observer which sets the bitfield
    auto fm = fastmarching::FastMarching<fastmarching::ObserverBitSet<mutil::HashMapMappingView<BigHashMap<uint64_t, float>>>>(report);
    fm.observer().resize(data.dims().size());
    fm.data(&data);
    fm.thresholds(min, max);
    fm.maxDistance(distance);

    // take biggest maxima first
    std::priority_queue<std::pair<float, uint64_t> > queue;

    report.start("Searching for seed points");
    // initialize the queue
    for(auto gridPos : data.gridPositions())
    {
        // ugly, should be improved
        const auto index = data.createLocation(gridPos).index();
        const auto center = data.get_unchecked(index);
        if(center < min){
            continue;
        }

        if (gridPos[0] != 0 && center < data.get_unchecked(index - 1))
        {
            continue;
        }
        if (gridPos[0] + 1 != dims[0] && center < data.get_unchecked(index + 1))
        {
            continue;
        }
        if (gridPos[1] != 0 && center < data.get_unchecked(index - dims[0]))
        {
            continue;
        }
        if (gridPos[1] + 1 != dims[1] && center < data.get_unchecked(index + dims[0]))
        {
            continue;
        }
        if (gridPos[2] != 0 && center < data.get_unchecked(index - dims[0] * dims[1]))
        {
            continue;
        }
        if (gridPos[2] + 1 != dims[2]  && center < data.get_unchecked(index + dims[0] * dims[1]))
        {
            continue;
        }

        // gridPos is a locale maxima
        queue.push({center, index});
    }

    const float queue_max = static_cast<float>(queue.size());
    const float points_max = static_cast<float>(maxPoints);

    // greedy algo
    auto points = std::vector<uint32_t>();
    while (!queue.empty() && points.size() < maxPoints)
    {
        std::pair<float, uint64_t> entry = queue.top();
        queue.pop();
        auto index = entry.second;

        // skip if covered by another seed point
        if(fm.observer().bitset()[index]){
            continue;
        }
        
        // set point and run fm
        points.push_back(index);
        fm.setStartVoxel(index);
        fm.march();

        // TODO: further expand bits by some neighborhood (done with border)

        // check border for new locale maxima
        for(auto border : fm.borderVoxels())
        {
            // skip if covered by another seed point
            if(fm.observer().bitset()[border]){
                continue;
            }

            // check for neigbhorhodd
            const auto gridPos = data.lattice().gridLocationFromCIndex(border, data.lattice().dims());
            const auto center = data.get_unchecked(border);
            if(center < min){
                continue;
            }
    
            if (gridPos[0] != 0 && center < data.get_unchecked(border - 1))
            {
                continue;
            }
            if (gridPos[0] + 1 != dims[0] && center < data.get_unchecked(border + 1))
            {
                continue;
            }
            if (gridPos[1] != 0 && center < data.get_unchecked(border - dims[0]))
            {
                continue;
            }
            if (gridPos[1] + 1 != dims[1] && center < data.get_unchecked(border + dims[0]))
            {
                continue;
            }
            if (gridPos[2] != 0 && center < data.get_unchecked(border - dims[0] * dims[1]))
            {
                continue;
            }
            if (gridPos[2] + 1 != dims[2]  && center < data.get_unchecked(border + dims[0] * dims[1]))
            {
                continue;
            }

            // locale maxima found
            queue.push({center, border});
        }

        if(points.size() % 128 == 0){
            report.update(fmt::format("Queue: ({:d})", queue.size()));
            // report.update(std::max((queue_max - queue.size()) / queue_max, points.size() / points_max));
        }
    }
    report.end();

    return points;
}