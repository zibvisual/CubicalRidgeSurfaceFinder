#pragma once

#include <vector>
#include <unordered_map>
#include <queue>
#include <optional>

#include <utils/Vec.hpp>
#include <utils/Dims.hpp>
#include <field/RawField.hpp>

#include "Graph.hpp"
#include "ChangeGraph.h"
#include <fastmarching/FastMarching.hpp>
#include <utils/Mapping.hpp>
// #include "MarchingCubes.hpp"
#include <surface/Surface.hpp>
#include <utils/Lattice.hpp>
#include <utils/ProgressbarReportHelper.hpp>

#include "Faces.hpp"
#include "FacesToCubicalMesh.hpp"
#include <field/RawField.hpp>
#include <utils/Range.hpp>

namespace futil
{
    class Face;
}

class HxSpatialGraph;

namespace ridgesurface
{

    enum SeedEnum
    {
        POINT,
        LINE
    };

    class SeedPoint
    {
    public:
        VecFloat point;

        SeedPoint(VecFloat point)
            : point(point)
        {
        }

        bool operator==(const SeedPoint& rhs) const;
    };

    class SeedLine
    {
    public:
        std::vector<VecFloat> points;

        SeedLine(std::vector<VecFloat> points)
            : points(points)
        {}

        bool operator==(const SeedLine& rhs) const;
    };

    class SeedPointIterator
    {
    public:
        SeedPointIterator(VecFloat point, Lattice lattice)
            : point(point)
            , lattice(lattice)
            , inner(true)
        {}

        std::optional<std::size_t>
        next()
        {
            if (inner)
            {
                std::size_t pointId = lattice.c_index(point);
                inner = false;
                return pointId;
            }
            return std::optional<std::size_t>();
        }

        using value_type = std::size_t;
        // ADD_RANGE_TO_RUST_ITERATOR(SeedPointIterator)

    private:
        VecFloat point;
        Lattice lattice;
        bool inner;
    };

    class SeedLineIterator
    {
    public:
        SeedLineIterator(std::vector<VecFloat>& points, Lattice lattice)
            : points(points)
            , index(0)
            , lattice(lattice)
        {}

        std::optional<std::size_t>
        next()
        {
            return index < points.size() ? lattice.c_index(points[index++]) : std::optional<std::size_t>();
        }

        using value_type = std::size_t;
        // ADD_RANGE_TO_RUST_ITERATOR(SeedLineIterator)

        // copy-assignment
        SeedLineIterator&
        operator=(const SeedLineIterator& sli)
        {
            // Guard self assignment
            // if (this == &sli)
            //     return *this;

            // we just copy the reference
            points = sli.points;
            lattice = sli.lattice;

            return *this;
        }

    private:
        std::vector<VecFloat>& points;
        std::size_t index;
        Lattice lattice;
    };

    class SeedIterator;

    class Seed
    {
        friend SeedIterator;

    protected:
        SeedEnum id;

        union
        {
            SeedPoint point;
            SeedLine line;
        };

        float max;

        void
        destroy_value()
        {
            switch (id)
            {
                case LINE:
                    (&line)->SeedLine::~SeedLine();
                    break;
                case POINT: // nothing to do
                default:
                    break;
            }
        }

        void
        copy_value(const Seed& seed)
        {
            switch (id)
            {
                case LINE:
                    line = seed.line;
                    break;
                case POINT:
                    point = seed.point;
                    break;
                default:
                    break;
            }
        }

    public:
        Seed(SeedPoint seedpoint, float distance)
            : id(POINT)
            , point({ seedpoint })
            , max(distance)
        {}

        Seed(SeedLine seedline, float distance)
            : id(LINE)
            , line({ seedline })
            , max(distance)
        {}

        Seed(const Seed& seed)
            : id(seed.id)
            , max(seed.max)
        {
            copy_value(seed);
        }

        ~Seed()
        {
            destroy_value();
        }

        // copy-assignment
        Seed&
        operator=(const Seed& seed)
        {
            // Guard self assignment
            if (this == &seed)
                return *this;
            destroy_value();
            id = seed.id;
            max = seed.max;
            copy_value(seed);
            return *this;
        }

        static Seed
        Seedpoint(VecFloat point, float distance)
        {
            return Seed(SeedPoint(point), distance);
        }

        static Seed
        Seedline(std::vector<VecFloat> points, float distance)
        {
            return Seed(SeedLine(points), distance);
        }

        SeedIterator getVoxelSources(Lattice lattice);

        float
        getDistance() const
        {
            return max;
        }

        void
        setDistance(float distance)
        {
            max = distance;
        }

        bool
        isPoint() const
        {
            return id == POINT;
        }

        bool
        isLine() const
        {
            return id == LINE;
        }

        // return a point representing the seed (for spatial graph, mostly used in seedpoints)
        VecFloat
        firstPoint() const
        {
            switch (id)
            {
                case LINE:
                    return line.points[0];
                case POINT:
                    return point.point;
                default: // not reachable
                    assert(false);
                    throw std::runtime_error("Not Reachable");
            }
        }

        bool
        operator==(const Seed& rhs) const
        {
            bool val = ((id == rhs.id) & (max == rhs.max));
            if (!val)
            {
                return false;
            }

            switch (id)
            {
                case LINE:
                    return line == rhs.line;
                case POINT:
                    return point == rhs.point;
                default: // not reachable
                    assert(false);
                    throw std::runtime_error("Not Reachable");
            }
        }
    };

    class SeedIterator
    {
    protected:
        SeedEnum id;

        union
        {
            SeedPointIterator point;
            SeedLineIterator line;
        };

        void
        destroy_value()
        {
            switch (id)
            {
                case POINT:
                    (&point)->SeedPointIterator::~SeedPointIterator();
                    break;
                case LINE:
                    (&line)->SeedLineIterator::~SeedLineIterator();
                    break;
                default:
                    break;
            }
        }

    public:
        SeedIterator(Seed seed, Lattice lattice)
            : id(seed.id)
        {
            switch (id)
            {
                case POINT:
                    point = SeedPointIterator(seed.point.point, lattice);
                    break;
                case LINE:
                    line = SeedLineIterator(seed.line.points, lattice);
                    break;
                default:
                    break;
            }
        }

        std::optional<std::size_t>
        next()
        {
            switch (id)
            {
                case POINT:
                    return point.next();
                case LINE:
                    return line.next();
                default:
                    return std::optional<std::size_t>();
            }
        }

        using value_type = std::size_t;

        ADD_RANGE_TO_RUST_ITERATOR(SeedIterator)
    };

    /**
     * @brief
     *
     * TODO: Use uuid to track seeds! This allows to remove ChangeGraph and makes updating a lot easier, as we use an indirection.
     * TODO: -> if some performance hit is felt, we can use a sort function to keep track of them.
     * TODO: m_graph is only necessary for the finalization step. Such, use it only there, patches shouldn't have an orienation.
     *
     */
    class CubicalRidgeSurfaceFinder
    {
    public:
        CubicalRidgeSurfaceFinder(progressbar::Progressbar& progress);
        ~CubicalRidgeSurfaceFinder(void);

        void setInput(RawField<float>* probability);
        void setThresholds(float min, float max);
        float getMin();
        float getMax();
        float getUndefinedValue();

        // void saveSeedpoints(FILE* fp);

        // from voxels // THIS IS MORE A HACK! INSTEAD THE SURFACE OBJECT SHOULD HAVE ALL THE DATA
        const Seed* getSeedFromVoxel(std::size_t id) const;
        std::optional<std::size_t> getSeedIndexFromVoxel(std::size_t id) const;
        float getDistance(VecFloat point);

        // indexing
        const Seed& getSeed(std::size_t i) const;

        // Same function again, but easier to use const function in non-const context
        const Seed&
        getConstSeed(std::size_t i) const
        {
            return getSeed(i);
        }

        Seed& getSeed(std::size_t i);

        // Seed changes
        std::size_t numOfSeeds() const;

        void addSeed(Seed seed);
        void addSeed(std::size_t i, Seed seed);
        // add a seedpoint automatically, if possible (if no seedpoint was set, returns false)
        int addSeed(float minDistanceRel, float maxDistanceNewSeed);
        void setSeed(std::size_t i, Seed seed);
        void removeSeed(std::size_t i);
        void clearSeeds();
        // used in load operations, such that we can use a diff algorithm
        void newSeeds(std::vector<Seed>& seeds);

        int removeSeedpointCandidates();
        int removeSeedpointCandidates(float minDistance);
        int removeSeedpointCandidatesWhenFullyCovered(std::size_t index);

        /**
         * @brief synonym for update
         */
        int compute();
        /**
         * @brief calculate anew from all seeds in module.
         *
         */
        int recalculate();
        /**
         * @brief Calculate the surface generated by merching the patches together.
         *
         */
        void finalize(surface::Surface* surface);//, HxSpatialGraph* graph);
        /**
         * @brief update module with the anew seeds.
         */
        int update();

        // Debug
        // void castToTwoLabels(unsigned char* labels);
        // void castCurrentSeedLabels(unsigned char* labels);
        void castToLabels(uint16_t* labels);
        // void castToSeedLabels(unsigned short* labels);
        // void castToSpatialGraph(HxSpatialGraph* graph);
        void castToTime(float* time) const;
        void castToDistance(float* distance) const;

        // void computeColoringOfSphere(std::size_t i, std::unordered_map<futil::Face, std::size_t>& map);
        // void computeIntegralCurves(std::size_t i, std::unordered_set<int64_t>& first, std::unordered_set<int64_t>& second);
        // void computeStartingGridpoints(std::size_t i, int64_t& first, int64_t& second);
        void computePatchSurface(surface::Surface& surface);

        const surface::Surface& patchedSurface() const {
            return m_patched_surface;
        }

        // std::vector<futil::Face> copyFacePatch();

        // float distanceOnSphere(futil::Face face);

    protected:
        void resetFields();
        /**
         * @brief calculated temporary patch of surface for given seed.
         *
         */
        void grow(std::size_t i);
        /**
         * @brief merges two fields together.
         *
         */
        void merge(std::size_t i);

        void merge_cum_fields(std::size_t i);
        void initialize_possible_neighbors(std::unordered_map<std::size_t, int>& neighbors);
        void track_neighbor_correspondence_by_volume(std::unordered_map<std::size_t, int>& neighbors);
        void update_neighbor_graph(std::size_t i, const std::unordered_map<std::size_t, int>& neighbors);
        void update_neighbor_graph_better(std::size_t i);

        void remove_non_touching_neighors();
        // not used right now, but could give better results
        void remove_non_touching_neighors_cube();
        void update_patch_orientations();

        void finalize_neighborhood_graph();

        // calculations done by update/recalculate

        /**
         * @brief Necessary for updateInsert and updateReplace (makes sure all necessary elements exist)
         *
         */
        void updateReserve();
        void updateInsert(std::size_t i);
        void updateInsertTask(std::size_t n, float progress);
        /**
         * @brief After clearing all necessary patches, one must call removeObsoleteTriangles on the surface. (To also clear triangles)
         *
         * @param i
         */
        void updateClear(std::size_t i);
        void updateClearTask(std::size_t n, float progress);
        void updateReplace(std::size_t left, std::size_t right);
        void updateReplaceTask(std::size_t n, float progress);

        void addPatch(std::size_t i);

        /**
         * @brief Given voxel ids, tells you if there exists a face inbetween
         *
         * @param i
         * @param j
         * @return true
         * @return false
         */
        bool faceBetweenVoxels(std::size_t i, std::size_t j) const;

        // progressbar
        progressbar::Progressbar& m_progressbar;

        // Seeds
        std::vector<Seed> m_seeds;
        // Edit status of seeds
        ChangeGraph m_seeds_change;
        // List of voxels changed by a seed
        std::vector<std::vector<int64_t>> m_seeds_voxels;
        // Neighborhood graph of seeds (seeds are only neighbors, if both their colors exist in the overlap)
        klenert::OrderedSparseGraph<bool> m_graph;

        // temporary memory used in grow()
        mutil::HashMapMappingView<SmallHashMap<int64_t, unsigned char>> m_label_view;
        std::unordered_set<int64_t> m_patch;
        int64_t m_first_grid_point;
        int64_t m_second_grid_point;

        // optional memory where data is saved (TODO: refactor, such that we have temporary memory itself as objects which then can get copied if necessary -> problematic as we use closures right now)
        std::unordered_map<Face, std::size_t>* ts_sphere_coloring;
        std::unordered_map<Face, float>* ts_sphere_distance;
        std::unordered_set<int64_t>* ts_integralcurve1;
        std::unordered_set<int64_t>* ts_integralcurve2;

        // Fast Marching object
        using MappingView = mutil::HashMapMappingView<SmallHashMap<int64_t, float>>;
        fastmarching::FastMarching<fastmarching::ObserverDistance<MappingView>> m_fm;

        // Inner labeling which can be transformed to create an unsigned char labeling for generate surface
        RawField<std::size_t> m_cum_label;

        // Input and Outputs
        RawField<float>* m_probability;
        RawField<float> m_cum_time;
        RawField<float> m_cum_distance;

        // patched surface is the surface one works on (no merging of surface itself (delay merging of everything and instead calculate all fields anew?))
        surface::Surface m_patched_surface;
        std::vector<bool> m_patch_orientation;

        FacesToCubicalMesh m_surface_writer;
        // candindates for seeding points (automatic) (patch_id -> candidates)
        std::vector<std::vector<int64_t>> m_seedpoints_candidates;

        // cached lattice
        Lattice m_lattice;

        // void initMemory();
        // void freeMemory();
    };

} // namespace ridgesurface