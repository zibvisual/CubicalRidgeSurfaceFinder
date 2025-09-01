#pragma once

#include <hxfastmarching/private/helpers/FastMarchingObserver.h>
#include <mutil/HashMap.hpp>
#include <mclib/McDim3l.h>

namespace fastmarching
{
    template <class MV>
    class FastMarchingAggregateObserver
    {
    public:
        // clang-format off
        using MappingView = MV;
        using Distance = DistanceEnabled;
        
    private:
        void seed(mcint64 id){
            m_time_data.set(id, 0.0f);
            m_distance_data.set(id, 0.0f);
        }
        void update(mcint64 id, float time, McDim3l dims){
            m_time_data.set(id, time);
            float dis = fastmarching::euclidean(id, dims, m_distance_data);
            m_distance_data.set(id, dis);
            m_max_distance = m_max_distance < dis ? dis : m_max_distance;
        }
        void fix(mcint64 id){
            m_fixed.push_back(id);
        }
        void candidate(mcint64 id){/* noop */}
        void bind(MappingView time, MappingView distance){
            m_time_data = time;
            m_distance_data = distance;
        }
        void reset(float default_value){
            m_time_data.clear(default_value);
            m_distance_data.clear(default_value);
            m_max_distance = 0.0f;
        }
        void start(){/* noop */}
        bool proceed(){return true;}
        void end(){/* noop */}
        void halt(){/* noop */}
        bool stop(bool userInterrupt, std::size_t volume, float time){
            return userInterrupt || time >= m_until_time || m_max_distance >= m_until_distance;
        }
        float progress(std::size_t volume, float time){
            float fac_time = time / m_until_time;
            float fac_distance = m_max_distance / m_until_distance;
            return std::max(fac_time, fac_distance);
        }
        MappingView& time() {return m_time_data;}
        MappingView& distance() {return m_distance_data;}
    public:
        const MappingView& time() const {return m_time_data;}
        const MappingView& distance() const {return m_distance_data;}
        void maxTime(float until_time){m_until_time = until_time;}
        void maxDistance(float until_distance){m_until_distance = until_distance;}
        float maxDistanceMarched() const {return m_max_distance;}

        // public methods
        const std::vector<mcint64>& fixed() const{
            return m_fixed;
        }
    private:
        MappingView m_distance_data;
        MappingView m_time_data;
        std::vector<mcint64> m_fixed;

        float m_until_time = FLT_MAX;
        float m_until_distance = FLT_MAX;
        float m_max_distance = 0.0f;

        template <class O, class R>
        friend class FastMarching;
        // clang-format on
    };
}
