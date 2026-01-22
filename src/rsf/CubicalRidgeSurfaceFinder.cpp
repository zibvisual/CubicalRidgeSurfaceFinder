#include "CubicalRidgeSurfaceFinder.h"
#include <utils/Vec.hpp>
#include <utils/Dims.hpp>

#include "Flooding.hpp"
#include "Amanatides.hpp"
#include "FaceGraph.hpp"
#include <utils/StructureTensor.hpp>
// #include <hxspatialgraph/internal/HxSpatialGraph.h>

#include <Eigen/Eigenvalues>

#include <stack>
#include <limits>
#include <unordered_set>

// to measure time elapsed
#include <chrono>

#include <fmt/base.h>
#include <fmt/format.h>

// PERF: To improve peformance, one might use the fact that we calculate face-surfaces (which always have at most 4 neighbor faces).
// PERF: Converting these Face-Surfaces to HalfEdgeMeshes can be faster than usual, as face-surfaces and their neighborhood can be constructed in such a way
// PERF: that the calculation of the smallest circles (that is the faces of the halfedgemeshes itself, as faces in the face-surface become points) are fast!
// PERF: That is: arriving at a face from a specific direction, one can just go always "to the right" or "to the left" direction neighbor to furhter go along
// PERF: the smallest circle. Doing this for all circles results in the mesh itself.
namespace
{
    inline uint64_t
    label_to_seed(uint64_t label)
    {
        return (label >> 1) - 1;
    }

    inline uint64_t
    seed_to_label(uint64_t seed, bool side)
    {
        return static_cast<uint64_t>((seed << 1) + 2 + static_cast<uint64_t>(side));
    }

    inline uint64_t
    seed_to_label(uint64_t seed, uint64_t label)
    {
        return static_cast<uint64_t>((seed << 1) + 2 + (label & 1));
    }

    inline uint64_t
    label_flip(uint64_t label)
    {
        return label ^ 1;
    }

} // namespace

namespace ridgesurface
{
    CubicalRidgeSurfaceFinder::CubicalRidgeSurfaceFinder(progressbar::Progressbar& progress)
        : m_progressbar(progress)
        , m_seeds()
        , m_seeds_update()
        , m_graph()
        , m_label_view(0)
        , ts_sphere_coloring(nullptr)
        , ts_sphere_distance(nullptr)
        , ts_integralcurve1(nullptr)
        , ts_integralcurve2(nullptr)
        , m_fm(progress)
        , m_cum_label()
        , m_probability(nullptr)
        , m_cum_time()
        , m_cum_distance()
        , m_surface_update()
        , m_patch_orientation()
        , m_surface_writer()
        , m_seedpoints_candidates()
        , m_lattice(Lattice(Dims(0,0,0)))
    {
        m_cum_label.setUndefinedValue(0);
        m_cum_time.setUndefinedValue(INFINITY);
        m_cum_distance.setUndefinedValue(INFINITY);
    }

    CubicalRidgeSurfaceFinder::~CubicalRidgeSurfaceFinder()
    {
    }

    void
    CubicalRidgeSurfaceFinder::setInput(RawField<float>* probability)
    {
        m_probability = probability;
        if(probability){
            m_lattice = probability->lattice();
        }else{
            m_lattice = Lattice();
        }

        m_cum_label.resize(m_lattice.dims());
        m_cum_time.resize(m_lattice.dims());
        m_cum_distance.resize(m_lattice.dims());

        resetFields();
    }

    void
    CubicalRidgeSurfaceFinder::setThresholds(float min, float max)
    {
        m_fm.thresholds(min, max);
    }

    float
    CubicalRidgeSurfaceFinder::getMin()
    {
        return m_fm.minThreshold();
    }

    float
    CubicalRidgeSurfaceFinder::getMax()
    {
        return m_fm.maxThreshold();
    }

    float
    CubicalRidgeSurfaceFinder::getUndefinedValue()
    {
        return m_fm.undefinedValue();
    }

    void
    CubicalRidgeSurfaceFinder::resetFields()
    {
        m_cum_label.fill();
        m_cum_time.fill();
        m_cum_distance.fill();
    }

    const CubicalRidgeSurfaceFinder& CubicalRidgeSurfaceFinder::as_const() const
    {
        return *this;
    }

    const Seed&
    CubicalRidgeSurfaceFinder::getSeed(uint64_t id) const
    {
        return m_seeds.at(id);
    }

    Seed&
    CubicalRidgeSurfaceFinder::getSeed(uint64_t id)
    {
        m_seeds_update.insert(id);
        return m_seeds.at(id);
    }

    const Seed*
    CubicalRidgeSurfaceFinder::getSeedFromVoxel(std::size_t id) const
    {
        if (m_cum_label.data()[id] < 2)
        {
            return nullptr;
        }
        auto seed = label_to_seed(m_cum_label.data()[id]);
        return &m_seeds.at(seed);
    }

    std::optional<uint64_t>
    CubicalRidgeSurfaceFinder::getSeedIndexFromVoxel(std::size_t id) const
    {
        if (m_cum_label.data()[id] < 2)
        {
            return std::optional<std::size_t>();
        }
        return label_to_seed(m_cum_label.data()[id]);
    }

    // Seed changes
    std::size_t
    CubicalRidgeSurfaceFinder::numOfSeeds() const
    {
        return m_seeds.size();
    }

    std::optional<uint64_t>
    CubicalRidgeSurfaceFinder::addSeed(float minDistanceRel, float maxDistanceNewSeed)
    {
        // check if we have a valid seedpoint candidate -> we want the biggest relative distance!
        int64_t max_id = -1;
        float max = 0;
        int index = -1;
        for (std::size_t i = 0; i < m_seedpoints_candidates.size(); ++i)
        {
            for (std::size_t j = 0; j < m_seedpoints_candidates[i].size(); ++j)
            {
                auto face = Face::fromIndex(m_seedpoints_candidates[i][j]);
                auto voxel = face.startId();
                auto distance = m_cum_distance.data()[voxel];
                auto maxDistance = m_seeds[i].getDistance();
                auto f = distance / maxDistance;
                if (f > max)
                {
                    max = f;
                    max_id = voxel;
                    index = i;
                }
            }
        }
        if (max < minDistanceRel || max_id < 0)
        {
            return {};
        }
        auto point = m_lattice.worldPosition(Lattice::gridLocationFromCIndex(static_cast<std::size_t>(max_id), m_lattice.dims()));
        return addSeed(Seed::Seedpoint(point, maxDistanceNewSeed));
    }

    int
    CubicalRidgeSurfaceFinder::removeSeedpointCandidates()
    {
        int cnt = 0;
        for (std::size_t i = 0; i < m_seedpoints_candidates.size(); ++i)
        {
            cnt += m_seedpoints_candidates[i].size();
            m_seedpoints_candidates[i].clear();
        }
        return cnt;
    }

    int
    CubicalRidgeSurfaceFinder::removeSeedpointCandidatesWhenFullyCovered(std::size_t index)
    {
        bool fully_covered = true;
        auto cnt = m_seedpoints_candidates[index].size();
        for (std::size_t j = 0; j < cnt; ++j)
        {
            auto face = Face::fromIndex(m_seedpoints_candidates[index][j]);
            auto voxel = face.startId();
            std::size_t seed = label_to_seed(m_cum_label.data()[voxel]);
            if (seed == index)
            {
                fully_covered = false;
            }
        }
        if (fully_covered)
        {
            m_seedpoints_candidates[index].clear();
            return cnt;
        }
        else
        {
            return 0;
        }
    }

    // this is not used, for good reasons... however the list could get out of hand... deletion may be necessary
    int
    CubicalRidgeSurfaceFinder::removeSeedpointCandidates(float minDistance)
    {
        int cnt = 0;
        for (std::size_t i = 0; i < m_seedpoints_candidates.size(); ++i)
        {
            for (std::size_t j = 0; j < m_seedpoints_candidates[i].size(); ++j)
            {
                auto face = Face::fromIndex(m_seedpoints_candidates[i][j]);
                auto voxel = face.startId();
                auto distance = m_cum_distance.data()[voxel];
                if (distance < minDistance)
                {
                    // delete this entry
                    m_seedpoints_candidates[i].erase(m_seedpoints_candidates[i].begin() + j);
                    --j;
                    ++cnt;
                }
            }
        }
        return cnt;
    }

    uint64_t
    CubicalRidgeSurfaceFinder::addSeed(Seed seed)
    {
        auto id = m_id_counter++;
        m_seeds.insert({id, seed});
        m_seeds_update.insert(id);
        return id;
    }

    void CubicalRidgeSurfaceFinder::updateSeed(uint64_t id, Seed seed)
    {
        m_seeds[id] = seed;
        m_seeds_update.insert(id);
    }

    bool
    CubicalRidgeSurfaceFinder::removeSeed(uint64_t id)
    {
        m_seeds_update.insert(id);
        return m_seeds.erase(id);
    }

    void
    CubicalRidgeSurfaceFinder::clearSeeds()
    {
        for(auto const& pair : m_seeds){
            m_seeds_update.insert(pair.first);
        }
        m_seeds.clear();
        // reset counter ?
        // m_id_counter = 0;
    }

    std::vector<uint64_t>
    CubicalRidgeSurfaceFinder::newSeeds(std::vector<Seed>& seeds)
    {
        // Hash the seeds
        auto map = std::unordered_map<Seed, uint64_t>();
        for(const auto& pair : m_seeds){
            map.insert({pair.second, pair.first});
        }

        // add all new seeds and get the ids of the already existing ones
        auto found = std::unordered_set<uint64_t>();
        auto list = std::vector<uint64_t>();
        for(const auto& seed : seeds){
            if(map.contains(seed)){
                auto id = map[seed];
                list.push_back(id);
                found.insert(id);
            }else{
                list.push_back(addSeed(seed));
            }
        }

        // remove all old seeds
        for(const auto& pair : map){
            if(!found.contains(pair.second))
            {
                removeSeed(pair.second);
            }
        }

        return list;
    }

    void
    CubicalRidgeSurfaceFinder::castToTime(float* time) const
    {
        for(std::size_t i = 0; i < m_lattice.size(); ++i){
            // TODO: we should instead use a get method from m_cum_time
            time[i] = m_cum_time.data()[i];
        }
    }

    void
    CubicalRidgeSurfaceFinder::castToDistance(float* distance) const
    {
        for(std::size_t i = 0; i < m_lattice.size(); ++i){
            // TODO: we should instead use a get method from m_cum_distance
            distance[i] = m_cum_distance.data()[i];
        }
    }

    // void
    // CubicalRidgeSurfaceFinder::castToTwoLabels(unsigned char* labels)
    // {
    //     std::size_t index = 0;
    //     for (int z = 0; z < m_dims[2]; ++z)
    //     {
    //         for (int y = 0; y < m_dims[1]; ++y)
    //         {
    //             for (int x = 0; x < m_dims[0]; ++x)
    //             {
    //                 if (m_cum_label[index] == 0)
    //                 {
    //                     labels[index] = 0;
    //                 }
    //                 else
    //                 {
    //                     bool valid = true;
    //                     std::size_t seed = label_to_seed(m_cum_label[index]);
    //                     if (x > 0)
    //                     {
    //                         std::size_t otherLabel = m_cum_label[index - 1];
    //                         std::size_t otherSeed = label_to_seed(otherLabel);
    //                         valid = valid ? (otherLabel == 0 || seed == otherSeed || m_graph.hasEdge(seed, otherSeed)) : false;
    //                     }
    //                     if (y > 0)
    //                     {
    //                         std::size_t otherLabel = m_cum_label[index - m_dims[0]];
    //                         std::size_t otherSeed = label_to_seed(otherLabel);
    //                         valid = valid ? (otherLabel == 0 || seed == otherSeed || m_graph.hasEdge(seed, otherSeed)) : false;
    //                     }
    //                     if (z > 0)
    //                     {
    //                         std::size_t otherLabel = m_cum_label[index - m_dims[0] * m_dims[1]];
    //                         std::size_t otherSeed = label_to_seed(otherLabel);
    //                         valid = valid ? (otherLabel == 0 || seed == otherSeed || m_graph.hasEdge(seed, otherSeed)) : false;
    //                     }
    //                     labels[index] = valid ? static_cast<unsigned char>((m_cum_label[index] & 1) + 1) : 0;
    //                 }
    //                 // labels[index] = m_cum_label[index] == 0 ? 0 : static_cast<unsigned char>((m_cum_label[index] & 0b1) + 1);
    //                 ++index;
    //             }
    //         }
    //     }
    // }

    // void
    // CubicalRidgeSurfaceFinder::castToSeedLabels(unsigned short* labels)
    // {
    //     for (std::size_t i = 0; i < m_lattice.size(); ++i)
    //     {
    //         labels[i] = m_cum_label[i] == 0 ? 0 : label_to_seed(m_cum_label[i]) + 1;
    //     }
    // }

    void
    CubicalRidgeSurfaceFinder::castToLabels(uint16_t* labels)
    {
        for (std::size_t i = 0; i < m_lattice.size(); ++i)
        {
            labels[i] = static_cast<uint16_t>(m_cum_label.data()[i]);
        }
    }

    // void
    // CubicalRidgeSurfaceFinder::castToSpatialGraph(HxSpatialGraph* graph)
    // {
    //     // graph->startEditing();
    //     graph->clear();

    //     // make the graph sparser
    //     finalize_neighborhood_graph();

    //     for (auto seed : m_seeds)
    //     {
    //         graph->addVertex(seed.firstPoint());
    //     }

    //     for (auto edge : m_graph.edges())
    //     {
    //         graph->addEdge(edge[0], edge[1]);
    //     }

    //     // graph->stopEditing();
    //     // graph->touch();
    // }

    // void
    // CubicalRidgeSurfaceFinder::castCurrentSeedLabels(unsigned char* labels)
    // {
    //     for (std::size_t i = 0; i < m_lattice.size(); ++i)
    //     {
    //         // auto it = m_label.find(i);
    //         // labels[i] = it == m_label.end() ? 0 : it->second;
    //         labels[i] = m_label_view.get(i);
    //     }
    // }

    // void
    // CubicalRidgeSurfaceFinder::computeColoringOfSphere(std::size_t i, std::unordered_map<Face, std::size_t>& map)
    // {
    //     ts_sphere_coloring = &map;
    //     grow(i);
    //     ts_sphere_coloring = nullptr;
    // }

    // void
    // CubicalRidgeSurfaceFinder::computeIntegralCurves(std::size_t i, std::unordered_set<int64_t>& first, std::unordered_set<int64_t>& second)
    // {
    //     ts_integralcurve1 = &first;
    //     ts_integralcurve2 = &second;
    //     grow(i);
    //     ts_integralcurve1 = nullptr;
    //     ts_integralcurve2 = nullptr;
    // }

    // void
    // CubicalRidgeSurfaceFinder::computeStartingGridpoints(std::size_t i, int64_t& first, int64_t& second)
    // {
    //     grow(i);
    //     first = m_first_grid_point;
    //     second = m_second_grid_point;
    // }

    // std::vector<Face>
    // CubicalRidgeSurfaceFinder::copyFacePatch()
    // {
    //     auto vec = std::vector<Face>();
    //     for (auto faceId : m_patch)
    //     {
    //         vec.push_back(Face::fromIndex(faceId));
    //     }
    //     return vec;
    // }

    // float
    // CubicalRidgeSurfaceFinder::distanceOnSphere(Face face)
    // {
    //     return m_fm.distance().get(face.startId());
    // }

    void
    CubicalRidgeSurfaceFinder::updateInsert(uint64_t id)
    {
        grow(id);
        // check m_patch
        // std::cout << "DEBUG: Found " << m_patch.size() << " faces!" << std::endl;
        merge(id);
        addPatch(id);
        // save the voxels we are looking at
        m_seeds_voxels[id].clear();
        for (auto pair : m_label_view.inner())
        {
            m_seeds_voxels[id].push_back(pair.first);
        }
        // other considerations:
        // save possible candidates -> this is done in addPatch!
        // m_graph is changed in update_neighbor_graph
    }

    void
    CubicalRidgeSurfaceFinder::updateClear(uint64_t id)
    {
        // surface update
        m_surface_update.removePatch(id);
        m_seeds_in_surface.erase(id);
        // fields
        for (int64_t j : m_seeds_voxels[id])
        {
            // this assertion might trigger sometimes...
            // assert(m_cum_label[i] != 0);
            auto seed = label_to_seed(m_cum_label.data()[j]);
            if (seed == id)
            {
                m_cum_label.data()[j] = m_cum_label.getUndefinedValue();
                m_cum_time.data()[j] = INFINITY;
                m_cum_distance.data()[j] = INFINITY;
            }
        }
        // graphs
        m_graph.removeNode(id);
        // voxel list
        m_seeds_voxels.erase(id);
        // candidates
        m_seedpoints_candidates.erase(id);
    }

    surface::SurfaceUpdate
    CubicalRidgeSurfaceFinder::calculate()
    {
        m_progressbar.start("Calculate Patches");

        // we assume m_seeds and m_seeds_change to be in sync
        // fields which must be updated accordingly:
        // - m_graph
        // - m_cum_*
        // - m_seeds_update
        // - m_seeds_in_surface
        // - m_surface_update
        // - m_seeds_voxels
        // - m_seedpoints_candidates

        // consume m_seeds_update
        auto counter = 0;
        auto size = m_seeds_update.size();
        auto iter = m_seeds_update.begin();
        while (iter != m_seeds_update.end()){
            const auto id = *iter;
            if (m_progressbar.stop_and_end()) return m_surface_update.build();
            if(m_seeds_in_surface.contains(id)){
                // remove current patch
                updateClear(id);
            }
            if(m_seeds.contains(id)){
                // add/update patch
                updateInsert(id);
            }
            iter = m_seeds_update.erase(iter);
            m_progressbar.update(++counter/static_cast<float>(size));
            m_progressbar.update(fmt::format("Calculate Patches ({:d} / {:d})", counter, size));
        }
        update_patch_orientations();
        m_progressbar.end();
        return m_surface_update.build();
    }

    void
    CubicalRidgeSurfaceFinder::update_patch_orientations()
    {
        // dfs on the graph
        std::unordered_map<uint64_t, bool> orientations;
        orientations.reserve(m_seeds.size());
        std::vector<uint64_t> stack;

        for (auto &pair : m_seeds)
        {
            auto id = pair.first;
            if (!orientations.contains(id))
            {
                // start search
                stack.push_back(id);
                orientations.insert({id, false});
                while (!stack.empty())
                {
                    auto current = stack.back();
                    bool current_orientation = orientations[current];
                    stack.pop_back();
                    // check all neighbors
                    for (auto edge : m_graph.neighbors(current))
                    {
                        if (!orientations.contains(edge.first))
                        {
                            bool new_orientation = current_orientation ^ edge.second;
                            orientations[edge.first] = new_orientation;
                            stack.push_back(edge.first);
                        }
                    }
                }
            }
        }

        // Flip Material if necessary:
        // DIFF orientation and m_patch_orientation -> add the diff to SurfaceUpdate
        for(auto pair : orientations){
            if (!m_patch_orientation.contains(pair.first) || m_patch_orientation[pair.first] != pair.second) {
                m_surface_update.addFlip(pair.first);
            }
        }
        m_patch_orientation = std::move(orientations);
    }

    void
    CubicalRidgeSurfaceFinder::addPatch(uint64_t id)
    {
        // Generate a graph of the connections of the faces of the patch (m_patch)
        auto face_graph = FaceGraph::createWithID(m_patch, m_lattice.dims());

        auto it = face_graph.borderBegin();
        auto end = face_graph.borderEnd();
        while (it != end)
        {
            m_seedpoints_candidates[id].push_back(*it);
            ++it;
        }

        // we may assume that all necessary patches exist
        m_surface_writer.populateSurfacePatch(m_surface_update, id, m_lattice, face_graph.begin().convert(), face_graph.end().convert());
        m_seeds_in_surface.insert(id);
    }

    surface::SurfaceUpdate
    CubicalRidgeSurfaceFinder::recalculate()
    {
        // clear data
        m_graph.clear();
        m_surface_writer.clear();
        m_surface_update.clear();
        m_seeds_update.clear();
        m_seedpoints_candidates.clear();
        m_seeds_voxels.clear();

        resetFields();

        int counter = 0;
        m_progressbar.start("Extracting Ridge Surface Patches");
        for (auto &pair : m_seeds)
        {
            updateInsert(pair.first);
            ++counter;
            if (m_progressbar.stop_and_end()) return m_surface_update.build();
            m_progressbar.update(static_cast<float>(counter) / m_seeds.size());
            m_progressbar.update(fmt::format("Calculating Patches ({:d} / {:d})", counter, m_seeds.size()));
        }
        update_patch_orientations();
        m_progressbar.end();
        return m_surface_update.build();
    }

    bool
    CubicalRidgeSurfaceFinder::faceBetweenVoxels(std::size_t i, std::size_t j) const
    {
        bool validLabels = m_cum_label.data()[i] > 0 && m_cum_label.data()[j] > 0;
        if (!validLabels)
        {
            return false;
        }
        auto seed1 = label_to_seed(m_cum_label.data()[i]);
        auto seed2 = label_to_seed(m_cum_label.data()[j]);
        // Test: Performance improvement?
        // if(seed1 == seed2){
        //     return (m_cum_label[i] % 2) != (m_cum_label[j] % 2);
        // }
        auto polarity = m_graph.getEdge(seed1, seed2);
        if (!polarity)
        {
            // if it is no neighbor, we do not generate any faces
            return false;
        }
        auto border = ((m_cum_label.data()[i] % 2) != (m_cum_label.data()[j] % 2)) ^ polarity.value();
        return border;
    }

    void
    CubicalRidgeSurfaceFinder::finalize(surface::StaticSurface* surface)//, HxSpatialGraph* graph)
    {
        m_progressbar.start("Finalizing Surface");
        surface->clear();
        // graph->clear();

        // not needed for Surface Finder. Only reduces the amount of intersection in the graph.
        // finalize_neighborhood_graph();

        // create a list of all faces, then populate m_surface with it
        auto faces = BigHashSet<Face>();
        const auto dims = m_lattice.dims();

        // x direction
        for (std::size_t z = 0; z < dims[2]; ++z)
        {
            for (std::size_t y = 0; y < dims[1]; ++y)
            {
                for (std::size_t x = 0; x < dims[0] - 1; ++x)
                {
                    std::size_t index = m_lattice.c_index(VecSize(x,y,z));
                    if (faceBetweenVoxels(index, index + 1))
                    {
                        // auto firstSeed = label_to_seed(m_cum_label[index]);
                        // auto secondSeed = label_to_seed(m_cum_label[index + 1]);
                        // if (firstSeed != secondSeed)
                        // {
                        //     graph->addEdge(firstSeed, secondSeed);
                        // }
                        auto orientation = m_patch_orientation[label_to_seed(m_cum_label.data()[index])];
                        if ((m_cum_label.data()[index] % 2 == 1) ^ orientation)
                        {
                            faces.insert(Face(index, Direction::RIGHT));
                        }
                        else
                        {
                            faces.insert(Face(index + 1, Direction::LEFT));
                        }
                    }
                }
            }
        }

        // y direction
        for (std::size_t z = 0; z < dims[2]; ++z)
        {
            for (std::size_t y = 0; y < dims[1] - 1; ++y)
            {
                for (std::size_t x = 0; x < dims[0]; ++x)
                {
                    std::size_t i = m_lattice.c_index(VecSize(x,y,z));
                    std::size_t j = i + dims[0];
                    if (faceBetweenVoxels(i, j))
                    {
                        auto orientation = m_patch_orientation[label_to_seed(m_cum_label.data()[i])];
                        if ((m_cum_label.data()[i] % 2 == 1) ^ orientation)
                        {
                            faces.insert(Face(i, Direction::UP));
                        }
                        else
                        {
                            faces.insert(Face(j, Direction::DOWN));
                        }
                    }
                }
            }
        }

        // z direction
        for (std::size_t z = 0; z < dims[2] - 1; ++z)
        {
            for (std::size_t y = 0; y < dims[1]; ++y)
            {
                for (std::size_t x = 0; x < dims[0]; ++x)
                {
                    std::size_t i = m_lattice.c_index(VecSize(x,y,z));
                    std::size_t j = i + dims[0] * dims[1];
                    if (faceBetweenVoxels(i, j))
                    {
                        auto orientation = m_patch_orientation[label_to_seed(m_cum_label.data()[i])];
                        if ((m_cum_label.data()[i] % 2 == 1) ^ orientation)
                        {
                            faces.insert(Face(i, Direction::FORWARD));
                        }
                        else
                        {
                            faces.insert(Face(j, Direction::BACKWARD));
                        }
                    }
                }
            }
        }

        auto mesh_helper = FacesToCubicalMesh();
        mesh_helper.populateSurface(surface, m_lattice, faces.begin(), faces.end());

        // create all graph vertices by going through the seed points
        // for (auto seed : m_seeds)
        // {
        //     graph->addVertex(seed.firstPoint());
        // }

        // Face neighbors[4];
        // for (auto face : faces)
        // {
        //     auto count = face.neighbors(faces, m_dims, neighbors);
        //     auto firstSeed = label_to_seed(m_cum_label[face.startId()]);

        //     for (std::size_t i = 0; i < count; ++i)
        //     {
        //         auto secondSeed = label_to_seed(m_cum_label[neighbors[i].startId()]);
        //         if (firstSeed != secondSeed)
        //         {
        //             graph->addEdge(firstSeed, secondSeed);
        //         }
        //     }
        // }
        m_progressbar.end();
    }

    void
    CubicalRidgeSurfaceFinder::merge(uint64_t id)
    {
        // For all neighbors, we calculate the correspondence of their labels to the current labels and save that
        auto neighbors = std::unordered_map<uint64_t, int>();
        initialize_possible_neighbors(neighbors);
        track_neighbor_correspondence_by_volume(neighbors);
        merge_cum_fields(id);
        update_neighbor_graph(id, neighbors);
        // update_neighbor_graph_better(i); -> this did not work...find out why?
    }

    void
    CubicalRidgeSurfaceFinder::initialize_possible_neighbors(std::unordered_map<uint64_t, int>& neighbors)
    {
        // check if a label field corresponds to a face of the patch -> if so, we found a neighbor!
        std::unordered_map<uint64_t, std::size_t> seed_counter;
        for (auto faceId : m_patch)
        {
            auto face = Face::fromIndex(faceId);
            // // check if face startId and endId have label from the same seed
            // auto startId = face.startId();
            // auto maybeEndId = face.validEndId(m_lattice.dims());
            // if (!maybeEndId)
            // {
            //     continue;
            // }
            // auto endId = maybeEndId.value();

            auto startLabel = m_cum_label.get_or(face.startVoxel(m_lattice.dims()));
            auto endLabel = m_cum_label.get_or(face.endVoxel(m_lattice.dims()));
            // auto startLabel = m_cum_label[startId];
            // auto endLabel = m_cum_label[endId];
            if ((startLabel != endLabel) && startLabel && endLabel)
            {
                auto seed = label_to_seed(startLabel);
                if (seed == label_to_seed(endLabel))
                {
                    seed_counter[seed] += 1;
                }
            }
        }

        // neighbors must at least intersect 1% of the patch
        for (auto seed : seed_counter)
        {
            if (static_cast<float>(seed.second) >= static_cast<float>(m_patch.size()) * 0.01)
            {
                neighbors.insert({ seed.first, 0 });
            }
        }
    }

    void
    CubicalRidgeSurfaceFinder::update_neighbor_graph_better(std::size_t seed)
    {
        m_graph.addNode(seed);
        m_graph.addEdge(seed, seed, false);

        auto face_cost = [this](Face face) {
            const auto cost = m_cum_time.get_unchecked(face.startVoxel(m_lattice.dims()));
            const auto cost2 = m_cum_time.get_or(face.endVoxel(m_lattice.dims()), cost);
            return cost + cost2;
        };

        // get min face by time
        int64_t min_face = -1;
        float min_cost = INFINITY;
        for (auto face : m_patch)
        {
            // time is average of upper and lower
            auto cost = face_cost(Face::fromIndex(face));
            if (cost < min_cost)
            {
                min_face = face;
            }
        }

        // grow patch out, always checking if the merged label field has our seed point.
        std::unordered_set<int64_t> capped_patch;
        std::vector<int64_t> border_stack;
        capped_patch.insert(min_face);
        border_stack.push_back(min_face);
        while (!border_stack.empty())
        {
            auto face_id = border_stack.back();
            border_stack.pop_back();

            // get neighbors of face
            auto face = Face::fromIndex(face_id);
            std::array<int64_t, 4> neighbors;
            auto neighbors_counter = face.neighbors(m_patch, m_lattice.dims(), neighbors.data());
            for (std::size_t i = 0; i < neighbors_counter; ++i)
            {
                auto neighbor_id = neighbors[i];

                if (capped_patch.find(neighbor_id) != capped_patch.end())
                {
                    continue;
                }

                // check if neighbor is in other seed point territory
                // one voxel on either side is enough, so we just check start_id of face
                auto voxel_id = Face::fromIndex(neighbor_id).startId();
                auto voxel_label = m_cum_label.data()[voxel_id];
                if (voxel_label == 0 || label_to_seed(voxel_label) == i)
                {
                    // grow
                    border_stack.push_back(neighbor_id);
                    capped_patch.insert(neighbor_id);
                }
                else
                {
                    // border
                    auto other_seed = label_to_seed(voxel_label);
                    // check orientation (correct orientation would mean that we have lower label value (0 bit))
                    bool orientation = !(voxel_label & 1);
                    m_graph.addEdge(seed, other_seed, orientation);
                    m_graph.addEdge(other_seed, seed, orientation);
                }
            }
        }
    }

    void
    CubicalRidgeSurfaceFinder::finalize_neighborhood_graph()
    {
        // remove_non_touching_neighors_cube();
    }

    void
    CubicalRidgeSurfaceFinder::track_neighbor_correspondence_by_volume(std::unordered_map<uint64_t, int>& neighbors)
    {
        // calculate correspondence for the orientation of edges
        for (auto pair : m_label_view.inner())
        {
            auto index = pair.first;
            auto cur_label = static_cast<uint64_t>(pair.second & 1);
            if (m_cum_label.data()[index] == 0)
            {
                continue;
            }
            auto cum_label = m_cum_label.data()[index] & 1;
            auto cum_seed = label_to_seed(m_cum_label.data()[index]);
            // check if we have cum_seed in neighbors
            if (neighbors.find(cum_seed) != neighbors.end())
            {
                neighbors[cum_seed] += (cur_label == cum_label) ? 1 : -1;
            }
        }
    }

    // void
    // CubicalRidgeSurfaceFinder::remove_non_touching_neighors()
    // {
    //     auto neighbors = std::vector<std::unordered_set<std::size_t>>();
    //     neighbors.reserve(m_seeds.size());
    //     for (std::size_t i = 0; i < m_seeds.size(); ++i)
    //     {
    //         neighbors.push_back(std::unordered_set<std::size_t>());
    //     }

    //     // find all touching neighbors
    //     int64_t voxel_neighbors[6];
    //     for (std::size_t i = 0; i < m_lattice.size(); ++i)
    //     {
    //         if (m_cum_label[i] == 0)
    //         {
    //             continue;
    //         }
    //         auto cum_seed = label_to_seed(m_cum_label[i]);
    //         // get voxel neighbors, insert all seeds in our set.
    //         auto neighbor_counter = futil::gridNeighbors6(m_dims, i, voxel_neighbors);
    //         for (std::size_t j = 0; j < neighbor_counter; ++j)
    //         {
    //             auto neighbor_label = m_cum_label[voxel_neighbors[j]];
    //             if (neighbor_label == 0)
    //             {
    //                 continue;
    //             }
    //             neighbors[cum_seed].insert(label_to_seed(neighbor_label));
    //         }
    //     }

    //     // filter out all edges which are non-touching
    //     m_graph.filter([&neighbors](std::size_t first, std::size_t second) { return neighbors[first].find(second) != neighbors[first].end(); });
    // }

    // void
    // CubicalRidgeSurfaceFinder::remove_non_touching_neighors_cube()
    // {
    //     auto neighbors = std::vector<std::unordered_map<std::size_t, std::size_t>>();
    //     neighbors.reserve(m_seeds.size());
    //     for (std::size_t i = 0; i < m_seeds.size(); ++i)
    //     {
    //         neighbors.push_back(std::unordered_map<std::size_t, std::size_t>());
    //     }

    //     std::array<std::size_t, 4> sources;
    //     std::size_t source_counter = 0;
    //     std::array<std::size_t, 8> current_labels;
    //     std::size_t current_counter = 0;
    //     // check a 2x2x2 cube -> there must the both labels of one seedpoint and then all other seedpoints inside are counted as neighbors -> if this happens at least a certain amount of time.
    //     for (std::size_t z = 0; z + 1 < static_cast<std::size_t>(m_dims[2]); ++z)
    //     {
    //         for (std::size_t y = 0; y + 1 < static_cast<std::size_t>(m_dims[1]); ++y)
    //         {
    //             for (std::size_t x = 0; x + 1 < static_cast<std::size_t>(m_dims[0]); ++x)
    //             {
    //                 // go through the cube
    //                 source_counter = 0;
    //                 current_counter = 0;
    //                 for (std::size_t zi = 0; zi < 2; ++zi)
    //                 {
    //                     for (std::size_t yi = 0; yi < 2; ++yi)
    //                     {
    //                         for (std::size_t xi = 0; xi < 2; ++xi)
    //                         {
    //                             // for all iterations equal
    //                             std::size_t current_label = m_cum_label[futil::gridPositionToIndex(x + xi, y + yi, z + zi, m_dims)];
    //                             if (current_label == 0 || mutil::contains(current_labels, current_label, current_counter))
    //                             {
    //                                 continue;
    //                             }
    //                             if (mutil::contains(current_labels, label_flip(current_label), current_counter))
    //                             {
    //                                 source_counter = mutil::set_insert(sources, label_to_seed(current_label), source_counter);
    //                             }
    //                             current_counter = mutil::set_insert(current_labels, current_label, current_counter);
    //                         }
    //                     }
    //                 }
    //                 // add neighbors from sources to current_labels
    //                 for (std::size_t i = 0; i < source_counter; ++i)
    //                 {
    //                     std::size_t src_seed = sources[i];
    //                     for (std::size_t j = 0; j < current_counter; ++j)
    //                     {
    //                         std::size_t dst_seed = label_to_seed(current_labels[j]);
    //                         auto found = neighbors[src_seed].find(dst_seed);
    //                         if (found == neighbors[src_seed].end())
    //                         {
    //                             neighbors[src_seed].insert({ dst_seed, 0 });
    //                         }
    //                         else
    //                         {
    //                             neighbors[src_seed][dst_seed] += 1;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }

    //     // filter out all edges which are non-touching
    //     m_graph.filter([&neighbors](std::size_t first, std::size_t second) { return neighbors[first].find(second) != neighbors[first].end() && neighbors[first][second] > 4; });
    // }

    void
    CubicalRidgeSurfaceFinder::update_neighbor_graph(uint64_t id, const std::unordered_map<uint64_t, int>& neighbors)
    {
        m_graph.addNode(id);
        m_graph.addEdge(id, id, false);

        // update graph
        for (auto pair : neighbors)
        {
            m_graph.addEdge(pair.first, id, pair.second <= 0);
            m_graph.addEdge(id, pair.first, pair.second <= 0);
        }
    }

    void
    CubicalRidgeSurfaceFinder::merge_cum_fields(uint64_t id)
    {
        const auto& time_view = m_fm.time();
        const auto& distance_view = m_fm.distance();
        // merge fields depending on time and orientation
        for (auto pair : m_label_view.inner())
        {
            auto index = pair.first;
            // std::size_t cur_label = orientation ? static_cast<std::size_t>((i << 1) + 2 + (pair.second & 1)) : static_cast<std::size_t>((i << 1) + 3 - (pair.second & 1));
            std::size_t cur_label = static_cast<std::size_t>((id << 1) + 2 + (pair.second & 1));
            auto time = time_view.get(index);
            m_cum_label.data()[index] = time < m_cum_time.data()[index] ? cur_label : m_cum_label.data()[index];
            m_cum_distance.data()[index] = time < m_cum_time.data()[index] ? distance_view.get(index) : m_cum_distance.data()[index];
            m_cum_time.data()[index] = time < m_cum_time.data()[index] ? time : m_cum_time.data()[index];
        }
    }

    void
    CubicalRidgeSurfaceFinder::grow(std::size_t id)
    {
        Seed seed = m_seeds[id];
        auto dims = m_lattice.dims();
        auto voxelSize = m_probability->lattice().voxelsize();

        const float gradient_sigma = 3.0f;
        const float tensor_sigma = 2.0f;

        // TODO: collect voxel sources by using Seed Iterators.
        VecFloat point = seed.firstPoint();
        VecSize gridPoint = m_lattice.gridLocation(point);

        // std::cout << "DEBUG: Grow from location " << gridPoint << std::endl;
        //  int64_t pointId = HxLattice3::gridPositionToIndex(dims, gridPoint);

        m_fm.data(m_probability);
        auto iter = seed.getVoxelSources(m_lattice);
        m_fm.setStartVoxels(iter.begin(), iter.end());
        m_fm.maxDistance(seed.getDistance());
        m_fm.march();

        const auto& time_view = m_fm.time();
        const auto& distance_view = m_fm.distance();

        float borderValue = m_fm.maxTimeMarched();

        // after marching, we want to find the OUTER surface and OUTSIDE voxels border
        std::unordered_set<int64_t> outsideBorder;
        std::unordered_set<Face> outsideSphere;

        // ALGO
        {
            Face neighbors[4];
            auto stackFaces = std::vector<Face>();
            // go from borderOfImage to middle and find first border (going lower than borderValue)
            // for now just take the left border image, e.g. 0
            VecSize currPos = gridPoint;
            // McVec3i currPos = McVec3i(gridPoint);
            currPos[0] = 0;
            while (currPos[0] < gridPoint[0])
            {
                if (time_view.get(m_lattice.c_index(currPos)) < borderValue)
                {
                    break;
                }
                ++currPos[0];
            }
            auto face = Face(m_lattice.c_index(currPos), Direction::LEFT);
            stackFaces.push_back(face);
            outsideSphere.insert(face);
            // create border
            while (!stackFaces.empty())
            {
                face = stackFaces.back();
                stackFaces.pop_back();
                borderNeighbor4Classifier([&time_view, dims, borderValue](int64_t voxel) { return time_view.get(voxel) <= borderValue; }, dims, face, neighbors);
                for (int j = 0; j < 4; ++j)
                {
                    auto neighbor = neighbors[j];
                    if (outsideSphere.find(neighbor) == outsideSphere.end())
                    {
                        outsideSphere.insert(neighbor);
                        stackFaces.push_back(neighbor);
                    }
                }
            }

            // fill outsideBorder
            for (auto face : outsideSphere)
            {
                auto option = face.validEndId(dims);
                if (option)
                {
                    outsideBorder.insert(option.value());
                }
            }
        }

        //  std::cout << "outside border collected" << std::endl;

        // After Marching, we want to flood from the start + tensor offset, until we reach the border -> value bigger than border

        // get tensor of pointId
        // TODO: we want to allow more than just float values!
        float* inputData = static_cast<float*>(m_probability->data());
        Eigen::Matrix3f tensor = structure_tensor(inputData, m_lattice.dims(), static_cast<VecInt>(gridPoint), gradient_sigma, tensor_sigma);


        // compute the eigen decomposition
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> solver;
        solver.computeDirect(tensor);
        auto eigenvectors = solver.eigenvectors();

        // Most important Eigenvector:
        VecFloat vector = VecFloat(eigenvectors(0, 2), eigenvectors(1, 2), eigenvectors(2, 2));
        // offset vector with some value (Eigenvalue could be used)
        vector = vector * 3 * voxelSize;
        // Find out which point that would be
        VecSize firstGridPoint = m_lattice.gridLocation(point + vector);
        VecSize secondGridPoint = m_lattice.gridLocation(point - vector);

        // TODO: TS: COPY
        m_first_grid_point = m_lattice.c_index(firstGridPoint);
        m_second_grid_point = m_lattice.c_index(secondGridPoint);

        // std::cout << "structure tensor calculated" << std::endl;

        // We do not know if these gridpoints are still inside the volume!
        auto breaker = [&outsideBorder](int64_t node) {
            return outsideBorder.find(node) != outsideBorder.end();
        };
        auto maybeFirstFace = amanatides(static_cast<VecInt>(gridPoint), static_cast<VecInt>(firstGridPoint), m_lattice, breaker);
        auto maybeSecondFace = amanatides(static_cast<VecInt>(gridPoint), static_cast<VecInt>(secondGridPoint), m_lattice, breaker);

        // Flood gridpoints to find out their crossings with the surface
        // auto stream_connection = [dims](int64_t node, int64_t* neighbors) { return futil::gridNeighbors6(dims, node, neighbors); };
        // auto stream_breaker = [time_view, borderValue](int64_t node) { return time_view.get(node) > borderValue; }; // TODO: || numberNeighbors < 6 (rewrite waterstream or copy it in here)
        auto stream_compare = [&time_view](int64_t left, int64_t right) {
            return time_view.get(left) < time_view.get(right);
        };
        auto waterstream = klenert::GridWaterStream<decltype(breaker), decltype(stream_compare)>(breaker, stream_compare, dims);

        auto firstFace = maybeFirstFace ? maybeFirstFace.value() : waterstream.stream(m_first_grid_point);

        // TODO: TS: COPY
        if (ts_integralcurve1)
        {
            // TODO: Bug of amanatides if gridPoint is exactly the center (coordinates 0,0,0), set IntegralCurves in tmp folder
            // ts_integralcurve1->insert(firstFace.startId());
            // ts_integralcurve1->insert(firstFace.endId(m_dims));
            for (auto& voxel : waterstream.taken())
            {
                ts_integralcurve1->insert(voxel);
            }
        }

        auto secondFace = maybeSecondFace ? maybeSecondFace.value() : waterstream.stream(m_second_grid_point);

        // TODO: TS: COPY
        if (ts_integralcurve2)
        {
            for (auto& voxel : waterstream.taken())
            {
                ts_integralcurve2->insert(voxel);
            }
        }

        waterstream.clear();
        // std::cout << "poles calculated" << std::endl;

        // now we can create the surface of m_time and split them with firstFace and secondFace
        // as we have only two sinks and their places are already known, we can use a flooding algorithm instead of merging of union-find structures.
        // flood surfaces until we found one of two sinks. Color every surface in flood (seen) in the labeling we chose.
        // after that, print the surface for now.

        // TODO: check if first and second Face are really at the border!

        // create a grah for the connectivity? Or use surfaceNeighbor?
        auto flood_connection = [&time_view, dims, borderValue](Face node, Face* neighbors) {
            borderNeighbor4Classifier([&time_view, dims, borderValue](int64_t voxel) { return time_view.get(voxel) <= borderValue; }, dims, node, neighbors);
            return 4;
        };
        auto flood_compare = [&distance_view](Face left, Face right) {
            return distance_view.get(left.startId()) > distance_view.get(right.startId());
        };
        auto flooding = klenert::Flooding<Face, decltype(flood_connection), 4, decltype(flood_compare)>(flood_connection, flood_compare);
        // flooding.flood({ firstFace, secondFace }, { (i << 1) + 1, (i << 1) + 2 }); // If we want to use this, we would need m_label to be bigger than unsigned char
        // std::vector<Face> flooding_startpoints = {firstFace, secondFace};
        // flooding.flood(flooding_startpoints);
        if(flooding.streamflood({ firstFace, secondFace }))
        {
            // TODO: make this error more formal (not just an std::cout output)
            std::cout << "could not generate patch " << id << " , as both poles contain the same locale minima" << std::endl;
        }

        auto& mapping = flooding.regions();

        // TODO: TS: copy mapping to map, if one should save the colored sphere
        if (ts_sphere_coloring)
        {
            for (auto& pair : mapping)
            {
                (*ts_sphere_coloring)[pair.first] = pair.second;
            }
        }

        // we found a contour on the sphere. We now want to make it continuous! -> ridge extraction near the contour through sampling -> read "GPU" paper for that. Sampling rules: point has to lie on an isosurface AND on a ridge line on the isosurface!
        // First find iso-surface (in specified area -> binary search(?)) then find ridge as 0-value of the projected gradient length ONTO iso-surface
        // -> create isosurface around contour, then create ridge on iso-surface by finding 0-value projected gradient on isosurface
        // Isosurface: direction to go to isosurface can be found through gradient.

        // the surface is colored. Now color each voxel inside the volume.
        // Instead of going through all voxels inside-out, we use faces and go outside-in -> does it really make a difference?

        // PERF: we might want to use only the startId to reduce the number of calls to validEndId
        // TODO: instead of time_view we want to make use of the same methods as in RawField -> we want to "get" from VecInt and VecSize....?
        // TODO: this also means we need different HashMaps for different ids? maybe....
        auto face_cost = [&time_view, dims](Face face) {
            const auto cost1 = time_view.get(face.startId());
            // TODO: time_view and other mappings should behave the same as RawField -> fit them together
            // PERF: when using bit_index/BLocation instead of CLocation, we do not have to check "contain" as we have enough "wrapping" on the indices
            // for now, to the steps manually
            const auto loc = face.endVoxel(dims);
            const auto cost2 = dims.contains(loc) ? time_view.get(static_cast<std::size_t>(loc.z()) * dims.x() * dims.y() + static_cast<std::size_t>(loc.y()) * dims.x() + static_cast<std::size_t>(loc.x())) : cost1;
            return cost1 + cost2;
            // auto endId = face.validEndId(dims);
            // return endId ? cost + time_view.get(endId.value()) : cost + cost;
        };
        // other face_cost method:
        // auto face_cost = [this](Face face) {
        //     const auto cost = m_cum_time.get_unchecked(face.startVoxel(m_lattice.dims()));
        //     const auto cost2 = m_cum_time.get_or(face.endVoxel(m_lattice.dims()), cost);
        //     return cost + cost2;
        // };
        auto face_compare = [face_cost](Face left, Face right) {
            return face_cost(left) < face_cost(right);
        };
        auto heap = boost::heap::priority_queue<Face, boost::heap::compare<decltype(face_compare)>>(face_compare);

        m_label_view.clear();
        m_patch.clear();

        // set all faces (could be done by constructor and iterator)
        for (auto pair : mapping)
        {
            heap.push(pair.first);
        }

        // PERF: we might want to "not" check the validity of the face endid -> as we only need to to check labeling here
        // PERF: another option: voxelid not dependent on dims (putting x,y,z together)
        // PERF: using another id would have some strong implications. All hashing would act on that.
        // PERF: only reading RawField would neccesitate constructing the index from the id.
        // PERF: we might also want to save the costs of each face so we do not have to recalculate....
        // PERF: or we wrap the "cost" field with a border of a voxel in each direction so we do not have to care for Out of bounds...
        int64_t neighbors[6];
        while (!heap.empty())
        {
            Face face = heap.top();
            heap.pop();
            if (m_label_view.get(face.startId()))
            {
                continue;
            }
            auto endId = face.validEndId(dims);
            // do we have a face going outside of the sphere?
            if (endId && m_label_view.get(endId.value()))
            {
                m_label_view.set(face.startId(), m_label_view.get(endId.value()));
            }
            else
            {
                m_label_view.set(face.startId(), mapping.at(face));
            }

            // TODO: face.startId() should be face.startLocation(). But then again, face is highly dependend on RawField....
            auto neighbors = m_probability->gridNeighbors6(RawField<float>::FieldLocation(face.startId()));
            // auto count = gridNeighbors6(dims, face.startId(), &neighbors[0]);
            for (std::size_t i = 0; i < 6; ++i)
            {
                if(neighbors[i].empty()){
                    continue;
                }
                auto neighbor_id = neighbors[i].index();
                if (!m_label_view.get(neighbor_id) && time_view.get(neighbor_id) <= borderValue)
                {
                    // go further inside
                    heap.push(Face(neighbor_id, face.startId(), dims));
                }
                else if (m_label_view.get(neighbor_id) < m_label_view.get(face.startId()) && time_view.get(neighbor_id) <= borderValue)
                {
                    // found a face of the surface patch
                    m_patch.insert(Face(neighbor_id, face.startId(), dims).toIndex());
                }
                else if (m_label_view.get(neighbor_id) > m_label_view.get(face.startId()) && time_view.get(neighbor_id) <= borderValue)
                {
                    // found a face of the surface patch
                    m_patch.insert(Face(face.startId(), neighbor_id, dims).toIndex());
                }
            }
        }
    }
} // namespace ridgesurface
