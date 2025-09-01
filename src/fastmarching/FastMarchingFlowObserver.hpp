#pragma once

#include <hxfastmarching/private/helpers/FastMarchingObserver.h>
#include <mutil/Mapping.hpp>
#include <mutil/HashMap.hpp>
#include <mclib/McDim3l.h>

namespace fastmarching
{
    class FastMarchingFlowObserver
    {
    public:
        using MappingView = mutil::HashMapMappingView<BigHashMap<mcint64, float>>;

        /// no distance field
        using Distance = DistanceDisabled;

        enum class StoppingCondition
        {
            NoStart,
            NoMoreData,
            TimeOut,
            TargetReached,
            UserInterrupt
        };

    private:
        // clang-format off
        /// called when a value is set to a non-default value in the initialization phase
        void seed(mcint64 id){m_time_data.set(id, 0.0f);}
        /// called when a value is updated
        void update(mcint64 id, float time, McDim3l dims){m_time_data.set(id, time);}
        /// called when a voxel is fixed
        void
        fix(mcint64 id)
        {
            m_seen.insert(id);
            if(m_target == id){
                m_status = StoppingCondition::TargetReached;
            }
        }
        /// called when a voxel is a candidate to be fixed next (border of the front)
        void candidate(mcint64 id){ /* noop */}
        // /// called when the MappingView is not owned
        // void bind(MappingView time){
        //     m_time_data = time;
        // }
        /// called when data fields must be reset
        void reset(float default_value){
            m_time_data.clear(default_value);
            m_seen.clear();
        }
        /// called at the start of the marching algorithm (march method)
        void
        start()
        {
            m_start_time = std::chrono::steady_clock::now();
            m_status = StoppingCondition::NoStart;
        }
        /// are we allowed to start marching?
        bool proceed(){ 
            return m_status == StoppingCondition::NoStart;
        }
        /// called at the end of the marching algorithm (march method)
        void end(){ /* noop */}
        /// called when no more voxels exists
        void halt() {
            m_status = StoppingCondition::NoMoreData;
        }
        /// if fast marching should be stopped
        bool
        stop(bool userInterrupt, std::size_t volume, float time)
        {
            if(userInterrupt)
            {
                m_status = StoppingCondition::UserInterrupt;
            }
            // only test timeout every 65536 cycles, as this calculation is expensive
            else if(((volume & 0xFFFF) == 0xFFFF) && (m_until_computation_time < std::chrono::duration<float>(std::chrono::steady_clock::now() - m_start_time).count())){
                m_status = StoppingCondition::TimeOut;
            }
            return m_status != StoppingCondition::NoStart;
        }

        /// calculate progress when FM would be finished
        float progress(std::size_t volume, float time){return 0.5f;}
        MappingView& time() {return m_time_data;}
    public:

        const MappingView& time() const {return m_time_data;}

        // clang-format on

        void
        setTarget(mcint64 id)
        {
            if (!m_seen.contains(id))
            {
                setTargetUnsafe(id);
            }
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

        // not used by Fast Marching but by FastMarchingFlow
        bool
        seen(mcint64 id) const
        {
            return m_seen.contains(id);
        }

        const BigHashSet<mcint64>& seen() const {
            return m_seen;
        }

        StoppingCondition
        status() const
        {
            return m_status;
        }

        void
        setTargetUnsafe(mcint64 id)
        {
            m_target = id;
            if (m_status == StoppingCondition::TargetReached)
            {
                m_status = StoppingCondition::NoStart;
            }
        }

    private:
        MappingView m_time_data;
        BigHashSet<mcint64> m_seen;
        mcint64 m_target;
        std::chrono::steady_clock::time_point m_start_time;
        float m_until_computation_time = FLT_MAX;
        StoppingCondition m_status;

        template <class O, class R>
        friend class FastMarching;
    };
}
