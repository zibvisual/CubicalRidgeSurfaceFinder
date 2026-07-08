#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <queue>

#include "Dims.hpp"
#include "UnionFind.hpp"

#include <fastmarching/FastMarching.hpp>
#include <fastmarching/FastMarchingBitSetObserver.hpp>

#include <field/GridPositionBoxIterator.hpp>

namespace {
    bool inside_box(VecSize pos, VecSize min, VecSize max){
        return pos[0] >= min[0] && pos[0] < max[0] && pos[1] >= min[1] && pos[1] < max[1] && pos[2] >= min[2] && pos[2] < max[2] ;
    }
}


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

std::vector<uint32_t> greedyFmSampling(progressbar::Reporter& report, const RawFieldView<float>& data, float distance, float min, float max, float threshold, std::size_t margin, std::size_t maxPoints = std::numeric_limits<std::size_t>::max())
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
    VecSize maxBox = (VecSize)((VecInt) data.dims() - VecInt(margin));
    for(auto gridPos : GridPositionBoxIterator(VecSize(margin), maxBox)){
        const auto i = data.createLocation(gridPos).index();
        float val = data.data()[i];
        if(val > threshold){
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

/**  Greedy Sampling which only considers local maximas to reduce queue size.
 * 
 */
std::vector<uint32_t> greedySmartFmSampling(progressbar::Reporter& report, const RawFieldView<float>& data, float distance, float min, float max, float threshold, std::size_t margin, std::size_t maxPoints = std::numeric_limits<std::size_t>::max())
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
    const VecSize maxBox = (VecSize)((VecInt) data.dims() - VecInt(margin));
    const VecSize minBox = VecSize(margin);
    for(auto gridPos : GridPositionBoxIterator(minBox, maxBox))
    {
        // ugly, should be improved
        const auto index = data.createLocation(gridPos).index();
        const auto center = data.get_unchecked(index);
        if(center <= threshold){
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

        // check border for new locale maxima
        for(auto border : fm.borderVoxels())
        {
            // skip if covered by another seed point
            if(fm.observer().bitset()[border]){
                continue;
            }

            // skip if near image border
            const auto gridPos = data.lattice().gridLocationFromCIndex(border, data.lattice().dims());
            if(!inside_box(gridPos, minBox, maxBox)){
                continue;
            }

            // check for neigbhorhood
            const auto center = data.get_unchecked(border);
            if(center <= threshold){
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

    // TODO: save bitset as image if wanted
    auto bits = RawField<uint8_t>(dims);
    for(auto index = 0; index < dims.size(); ++index){
        bits.data()[index] = (uint8_t) fm.observer().bitset()[index];
    }
    bits.save("./debug_tmp/exclusion_zone.nrrd");

    return points;
}

/**
 * Greedy sampling which also grows further out than FM itself.
 */
std::vector<uint32_t> greedyShellSampling(progressbar::Reporter& report, const RawFieldView<float>& data, float distance, float min, float max, float threshold, std::size_t margin, std::size_t shell = 5, std::size_t maxPoints = std::numeric_limits<std::size_t>::max())
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
    VecSize maxBox = (VecSize)((VecInt) data.dims() - VecInt(margin));
    VecSize minBox = VecSize(margin);
    for(auto gridPos : GridPositionBoxIterator(minBox, maxBox))
    {
        // ugly, should be improved
        const auto index = data.createLocation(gridPos).index();
        const auto center = data.get_unchecked(index);
        if(center <= threshold){
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

        // grow shell-value voxels further out
        auto buffer1 = std::vector<uint64_t>();
        auto buffer2 = std::vector<uint64_t>();

        auto* current_layer = &buffer1;
        auto* next_layer = &buffer2;
        for(auto border : fm.borderVoxels())
        { 
            // skip if covered by another seed point
            if(fm.observer().bitset()[border]){
                continue;
            }

            // skip if near image border
            const auto gridPos = data.lattice().gridLocationFromCIndex(border, data.lattice().dims());
            if(!inside_box(gridPos, minBox, maxBox)){
                continue;
            }

            next_layer->push_back(border);
            fm.observer().bitset()[border] = true;
        }

        std::size_t shell_size = 0;
        while(shell_size < shell){
            auto* tmp = current_layer;
            current_layer = next_layer;
            next_layer = tmp;
            next_layer->clear();

            // grow one layer
            while(!current_layer->empty()){
                auto voxel = current_layer->back();
                current_layer->pop_back();

                // check for neighborhood
                const auto gridPos = data.lattice().gridLocationFromCIndex(voxel, data.lattice().dims());

                const auto left = voxel - 1;
                const auto right = voxel + 1;
                const auto up = voxel - dims[0];
                const auto down = voxel + dims[0];
                const auto above = voxel - dims[0] * dims[1];
                const auto below = voxel + dims[0] * dims[1];

                if (gridPos[0] != 0 && data.get_unchecked(left) > threshold && !fm.observer().bitset()[left])
                {

                    next_layer->push_back(left);
                    fm.observer().bitset()[left] = true;
                }
                if (gridPos[0] + 1 != dims[0] && data.get_unchecked(right) > threshold && !fm.observer().bitset()[right])
                {
                    next_layer->push_back(right);
                    fm.observer().bitset()[right] = true;
                }
                if (gridPos[1] != 0 && data.get_unchecked(up) > threshold && !fm.observer().bitset()[up])
                {
                    next_layer->push_back(up);
                    fm.observer().bitset()[up] = true;
                }
                if (gridPos[1] + 1 != dims[1] && data.get_unchecked(down) > threshold && !fm.observer().bitset()[down])
                {
                    next_layer->push_back(down);
                    fm.observer().bitset()[down] = true;
                }
                if (gridPos[2] != 0 && data.get_unchecked(above) > threshold && !fm.observer().bitset()[above])
                {
                    next_layer->push_back(above);
                    fm.observer().bitset()[above] = true;
                }
                if (gridPos[2] + 1 != dims[2] && data.get_unchecked(below) > threshold && !fm.observer().bitset()[below])
                {
                    next_layer->push_back(below);
                    fm.observer().bitset()[below] = true;
                }
            }

            ++shell_size;
        }

        // check border for new locale maxima
        for(auto border : *next_layer)
        {
            // skip if covered by another seed point
            if(fm.observer().bitset()[border]){
                continue;
            }

            // check for neigbhorhood
            const auto gridPos = data.lattice().gridLocationFromCIndex(border, data.lattice().dims());
            const auto center = data.get_unchecked(border);
            if(center <= threshold){
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
    
    // TODO: save bitset as image if wanted
    // auto bits = RawField<uint8_t>(dims);
    // for(auto index = 0; index < dims.size(); ++index){
    //     bits.data()[index] = (uint8_t) fm.observer().bitset()[index];
    // }
    // bits.save("./debug_tmp/exclusion_zone.nrrd");

    return points;
}

/**
 * Greedy sampling which also grows further out than FM itself. May not fully cover the data (even if unklikely).
 */
std::vector<uint32_t> greedyWaterflowSampling(progressbar::Reporter& report, const RawFieldView<float>& data, float distance, float min, float max, float threshold, std::size_t margin, std::size_t maxPoints = std::numeric_limits<std::size_t>::max())
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
    VecSize maxBox = (VecSize)((VecInt) data.dims() - VecInt(margin));
    VecSize minBox = VecSize(margin);
    for(auto gridPos : GridPositionBoxIterator(minBox, maxBox))
    {
        // ugly, should be improved
        const auto index = data.createLocation(gridPos).index();
        const auto center = data.get_unchecked(index);
        if(center <= threshold){
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
        // Idea: allow moving "down" from the potential until everything > min is covered (border is enough)

        // grow downwards
        // we use a heap so we now while filling up which voxels are border voxels
        auto cmp = [](std::pair<uint64_t, float> left, std::pair<uint64_t, float> right) { return left.second < right.second; };
        auto heap = std::priority_queue<std::pair<uint64_t, float>, std::vector<std::pair<uint64_t, float>>, decltype(cmp)>();
        auto new_border = std::vector<uint64_t>();
        for(auto border : fm.borderVoxels())
        { 
            // skip if covered by another seed point
            if(fm.observer().bitset()[border]){
                continue;
            }

            // skip if near image border
            const auto gridPos = data.lattice().gridLocationFromCIndex(border, data.lattice().dims());
            if(!inside_box(gridPos, minBox, maxBox)){
                continue;
            }

            heap.push(std::make_pair(border, data.get_unchecked(border)));
            fm.observer().bitset()[border] = true;
        }

        while(!heap.empty()){
            auto pair = heap.top();
            auto voxel = pair.first;
            auto center = pair.second;
            heap.pop();

            // check for neighborhood
            const auto gridPos = data.lattice().gridLocationFromCIndex(voxel, data.lattice().dims());

            const auto left = voxel - 1;
            const auto right = voxel + 1;
            const auto up = voxel - dims[0];
            const auto down = voxel + dims[0];
            const auto above = voxel - dims[0] * dims[1];
            const auto below = voxel + dims[0] * dims[1];

            if (gridPos[0] != 0 && data.get_unchecked(left) > threshold && !fm.observer().bitset()[left])
            {
                auto val = data.get_unchecked(left);
                if(center > val){
                    heap.push(std::make_pair(left, val));
                    fm.observer().bitset()[left] = true;
                }else{
                    new_border.push_back(left);
                }
            }
            if (gridPos[0] + 1 != dims[0] && data.get_unchecked(right) > threshold && !fm.observer().bitset()[right])
            {
                auto val = data.get_unchecked(right);
                if(center > val){
                    heap.push(std::make_pair(right, val));
                    fm.observer().bitset()[right] = true;
                }else{
                    new_border.push_back(right);
                }
            }
            if (gridPos[1] != 0 && data.get_unchecked(up) > threshold && !fm.observer().bitset()[up])
            {
                auto val = data.get_unchecked(up);
                if(center > val){
                    heap.push(std::make_pair(up, val));
                    fm.observer().bitset()[up] = true;
                }else{
                    new_border.push_back(up);
                }
            }
            if (gridPos[1] + 1 != dims[1] && data.get_unchecked(down) > threshold && !fm.observer().bitset()[down])
            {
                auto val = data.get_unchecked(down);
                if(center > val){
                    heap.push(std::make_pair(down, val));
                    fm.observer().bitset()[down] = true;
                }else{
                    new_border.push_back(down);
                }
            }
            if (gridPos[2] != 0 && data.get_unchecked(above) > threshold && !fm.observer().bitset()[above])
            {
                auto val = data.get_unchecked(above);
                if(center > val){
                    heap.push(std::make_pair(above, val));
                    fm.observer().bitset()[above] = true;
                }else{
                    new_border.push_back(above);
                }
            }
            if (gridPos[2] + 1 != dims[2] && data.get_unchecked(below) > threshold && !fm.observer().bitset()[below])
            {
                auto val = data.get_unchecked(below);
                if(center > val){
                    heap.push(std::make_pair(below, val));
                    fm.observer().bitset()[below] = true;
                }else{
                    new_border.push_back(below);
                }
            }
        }

        // check border for new locale maxima
        for(auto border : new_border)
        {
            // skip if covered by another seed point
            if(fm.observer().bitset()[border]){
                continue;
            }

            // check for neigbhorhood
            const auto gridPos = data.lattice().gridLocationFromCIndex(border, data.lattice().dims());
            const auto center = data.get_unchecked(border);
            if(center <= threshold){
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

    // TODO: save bitset as image if wanted
    // auto bits = RawField<uint8_t>(dims);
    // for(auto index = 0; index < dims.size(); ++index){
    //     bits.data()[index] = (uint8_t) fm.observer().bitset()[index];
    // }
    // bits.save("./debug_tmp/exclusion_zone.nrrd");

    return points;
}