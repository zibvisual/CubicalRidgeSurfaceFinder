#pragma once

#include <hxfastmarching/private/helpers/FastMarchingObserver.h>
#include <hxfastmarching/private/helpers/GeodesicDistanceNeighbor.h>
#include <mutil/Mapping.hpp>
#include <mutil/HashMap.hpp>
#include <mclib/McDim3l.h>

namespace fastmarching
{
    class FastMarchingGeodesicDistanceObserver
    {
    public:
        using MappingView = mutil::HashMapMappingView<BigHashMap<mcint64, float>>;
        using Distance = DistanceEnabled;

        enum class StoppingCondition
        {
            Running,
            NoMoreData,
            TimeOut,
            NeighborsFull,
            UserInterrupt
        };

    private:
        /// called when a value is set to a non-default value in the initialization phase
        void
        seed(mcint64 id)
        {
            m_time_data.set(id, 0.0f);
            m_distance_data.set(id, 0.0f);
        }

        /// called when a value is updated
        void
        update(mcint64 id, float time, McDim3l dims)
        {
            m_time_data.set(id, time);
            float dis = fastmarching::euclidean(id, dims, m_distance_data);
            m_distance_data.set(id, dis);
        }

        /// called when a voxel is fixed
        void
        fix(mcint64 id)
        {
            if (m_targets.find(id) != m_targets.end())
            {
                m_neighbors.push_back(GeodesicDistanceNeighbor { static_cast<std::size_t>(id), m_distance_data.get(id) });
            }
        }

        /// called when a voxel is a candidate to be fixed next (border of the front)
        void
        candidate(mcint64 id)
        { /* noop */
        }

        /// called to bind data (bindData method)
        void
        bind(MappingView time, MappingView distance)
        {
            m_time_data = time;
            m_distance_data = distance;
        }

        /// called when data fields must be reset
        void
        reset(float default_value)
        {
            m_time_data.clear(default_value);
            m_distance_data.clear(default_value);
            m_neighbors.clear();
        }

        /// called at the start of the marching algorithm (march method)
        void
        start()
        {
            m_start_time = std::chrono::steady_clock::now();
            m_status = StoppingCondition::Running;
        }

        /// are we allowed to start marching?
        bool
        proceed()
        {
            return true;
        }

        /// called at the end of the marching algorithm (march method)
        void
        end()
        { /* noop */
        }

        /// called when no more voxels exists
        void
        halt()
        {
            m_status = StoppingCondition::NoMoreData;
        }

        /// if fast marching should be stopped
        bool
        stop(bool userInterrupt, std::size_t volume, float time)
        {
            if (userInterrupt)
            {
                m_status = StoppingCondition::UserInterrupt;
            }
            // only test timeout every 65536 cycles, as this calculation is expensive
            else if (((volume & 0xFFFF) == 0xFFFF) && (m_until_computation_time < std::chrono::duration<float>(std::chrono::steady_clock::now() - m_start_time).count()))
            {
                m_status = StoppingCondition::TimeOut;
            }
            else if (m_neighbors.size() >= m_until_size_neighbors)
            {
                m_status = StoppingCondition::NeighborsFull;
            }
            return m_status != StoppingCondition::Running;
        }

        /// calculate progress when FM would be finished
        float
        progress(std::size_t volume, float time)
        {
            float fac_time = time / m_until_time;
            float fac_comp_time = std::chrono::duration<float>(std::chrono::steady_clock::now() - m_start_time).count() / m_until_computation_time;
            float fac_targets = static_cast<float>(m_neighbors.size()) / m_until_size_neighbors;
            return std::max(std::max(fac_time, fac_comp_time), fac_targets);
        }

        MappingView&
        time()
        {
            return m_time_data;
        }

        MappingView&
        distance()
        {
            return m_distance_data;
        }

    public:
        const MappingView&
        time() const
        {
            return m_time_data;
        }

        const MappingView&
        distance() const
        {
            return m_distance_data;
        }

        void
        maxComputationTime(float until_computation_time)
        {
            m_until_computation_time = until_computation_time;
        }

        float
        maxComputationTime() const
        {
            return m_until_computation_time;
        }

        void
        maxTime(float until_time)
        {
            m_until_time = until_time;
        }

        float
        maxTime() const
        {
            return m_until_time;
        }

        void
        maxNeighbors(std::size_t neighbors)
        {
            m_until_size_neighbors = neighbors;
        }

        std::size_t
        maxNeighbors() const
        {
            return m_until_size_neighbors;
        }

        StoppingCondition
        status() const
        {
            return m_status;
        }

        const std::vector<GeodesicDistanceNeighbor>&
        neighbors() const
        {
            return m_neighbors;
        }

    private:
        MappingView m_time_data;
        MappingView m_distance_data;
        BigHashSet<mcint64> m_targets;
        std::chrono::steady_clock::time_point m_start_time;
        float m_until_computation_time = FLT_MAX;
        float m_until_time = FLT_MAX;
        std::size_t m_until_size_neighbors = 20;
        StoppingCondition m_status;

        std::vector<GeodesicDistanceNeighbor> m_neighbors;

        template <class O, class R>
        friend class FastMarching;
    };
}
