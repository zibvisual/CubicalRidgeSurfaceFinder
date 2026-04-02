#pragma once

#include "FastMarchingUtilities.hpp"
#include "FastMarching.hpp"
#include <fastmarching/FastMarchingObserver.hpp>

#include <utils/Dims.hpp>

#include <chrono>
#include <cfloat>

namespace fastmarching
{
    template<class MV>
    class ObserverBitSet {
    public:
        using MappingView = MV;
        using Distance = DistanceEnabled;
        
    private:
        void seed(int64_t id){
            m_time_data.set(id, 0.0f);
            m_distance_data.set(id, 0.0f);
        }
        void update(int64_t id, float time, Dims dims){
            m_time_data.set(id, time);
            // TODO: grid location from index conversion should be easier...
            // float dis = fastmarching::euclidean(static_cast<VecInt>(RawField<float>::FieldLocation(id).location(dims)), dims, m_distance_data);
            float dis = fastmarching::euclidean_view(static_cast<VecInt>(Lattice::gridLocationFromCIndex(id, dims)), dims, m_distance_data);
            m_distance_data.set(id, dis);
            m_max_distance = m_max_distance < dis ? dis : m_max_distance;
        }
        void fix(int64_t id){
            m_bitset[id] = true;
        }
        void candidate(int64_t id){/* noop */}
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
        MappingView& distance() {return m_distance_data;}
        MappingView& time() {return m_time_data;}
    public:
        const MappingView& distance() const {return m_distance_data;}
        const MappingView& time() const {return m_time_data;}
        void maxTime(float until_time){m_until_time = until_time;}
        void maxDistance(float until_distance){m_until_distance = until_distance;}
        float maxDistanceMarched() const {return m_max_distance;}

        void resize(std::size_t size){
            m_bitset.resize(size);
        }
        const std::vector<bool>& bitset() const
        {
            return m_bitset;
        }
        
    private:
        MappingView m_time_data;
        MappingView m_distance_data;
        std::vector<bool> m_bitset;

        float m_until_time = FLT_MAX;
        float m_until_distance = FLT_MAX;
        float m_max_distance = 0.0f;

        template <class O>
        friend class FastMarching;
    };

    // clang-format on
}
