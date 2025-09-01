#pragma once

#include "FastMarchingUtilities.hpp"
#include "FastMarching.hpp"

#include <utils/Dims.hpp>

#include <chrono>
#include <cfloat>

namespace fastmarching
{

    // clang-format off
    // Flag classes
    class DistanceEnabled {};
    class DistanceDisabled {};


        /**
     * @brief Observer with all features.
     * 
     * @tparam MappingView 
     */
    template<class MV, class Set>
    class ObserverAll {
    public:
        using MappingView = MV;
        /// distance field exists
        using Distance = DistanceEnabled;
    private:
        /// called when a starting point should be added
        void seed(int64_t id){
            m_distance_data.set(id, 0.0f);
            m_time_data.set(id, 0.0f);
        }
        /// called when a value is updated
        void update(int64_t id, float time, Dims dims){
            m_time_data.set(id, time);
            // this is so ugly as we currently need to differentiate between mappingviews and fields. We will merge both APIs later on.
            float dis = fastmarching::euclidean_view(static_cast<VecInt>(RawField<float>::FieldLocation(id).location(dims)), dims, m_distance_data);
            m_distance_data.set(id, dis);
            m_max_distance = m_max_distance < dis ? dis : m_max_distance;
        }
        /// called when a voxel is fixed
        void fix(int64_t id){
            m_targets.erase(id);
        }
        /// called when a voxel is a candidate to be fixed next (border of the front)
        void candidate(int64_t id){/* noop */}
        /// called to bind data (bindData method)
        void bind(MappingView time, MappingView distance){
            m_time_data = time;
            m_distance_data = distance;
        }
        /// called in resetDataFields()
        void reset(float default_value){
            m_time_data.clear(default_value);
            m_distance_data.clear(default_value);
            m_max_distance = 0.0f;
        }
        /// called at the start of the marching algorithm (march method)
        void start(){}
        /// can we start the marching algorithm?
        bool proceed(){
            m_start_time = std::chrono::steady_clock::now();
            return true;
        }
        /// called at the end of the marching algorihtm (march method)
        void end(){/* noop */}
        /// called when no more voxels exists
        void halt(){/* noop */}
        /// if fast marching should be stopped
        bool stop(bool userInterrupt, std::size_t volume, float time){
            return userInterrupt
                || (m_targets.empty() && m_number_of_targets > 0) 
                || time >= m_until_time 
                || m_max_distance >= m_until_distance
                // only test timeout every 65536 cycles, as this is expensive
                || (((volume & 0xFFFF) == 0xFFFF) && (std::chrono::duration<float>(std::chrono::steady_clock::now() - m_start_time).count() >= m_until_computation_time));
        }
        /// calculate progress when FM would be finished
        float progress(std::size_t volume, float time){
            float fac_time = time / m_until_time;
            float fac_distance = m_max_distance / m_until_distance;
            float max = std::max(fac_time, fac_distance);
            if(m_number_of_targets > 0){
                float fac_targets = 1.0f - (static_cast<float>(m_targets.size()) / static_cast<float>(m_number_of_targets));
                max = std::max(max, fac_targets);
            }
            return max;
        }
        MappingView& time() {return m_time_data;}
        MappingView& distance() {return m_distance_data;}
    public:
        // special methods for distanceEnabled
        const MappingView& time() const {return m_time_data;}
        const MappingView& distance() const {return m_distance_data;}
        // march_until_* methods
        void maxTime(float until_time){m_until_time = until_time;}
        void maxDistance(float until_distance){m_until_distance = until_distance;}
        float maxDistanceMarched() const {return m_max_distance;}
        template <class InputIter>
        void setTargets(InputIter first, InputIter last){
            m_targets.clear();
            m_targets.insert(first, last);
            m_number_of_targets = m_targets.size();
        }
        void addTarget(int64_t id){
            m_targets.insert(id);
            ++m_number_of_targets;
        }
        void maxComputationTime(float until_computation_time){m_until_computation_time = until_computation_time;}
    private:
        MappingView m_distance_data;
        MappingView m_time_data;
        Set m_targets;
        std::chrono::steady_clock::time_point m_start_time;

        float m_until_time = FLT_MAX;
        float m_until_distance = FLT_MAX;
        float m_until_computation_time = FLT_MAX;

        float m_max_distance = 0.0f;
        std::size_t m_number_of_targets = 0;

        template <class O, class R>
        friend class FastMarching;
    };

    template<class MV, class Set>
    class ObserverTargetSet {
    public:
        using MappingView = MV;
        using Distance = DistanceDisabled;
    private:
        void seed(int64_t id){m_time_data.set(id, 0.0f);}
        void update(int64_t id, float time, Dims dims){m_time_data.set(id, time);}
        void fix(int64_t id){
            // m_number_of_targets -= m_targets.erase(id);
            m_targets.erase(id);
        }
        void candidate(int64_t id){/* noop */}
        void bind(MappingView time){
            m_time_data = time;
        }
        void reset(float default_value){
            m_time_data.clear(default_value);
        }
        void start(){/* noop */}
        bool proceed(){return true;}
        void end(){/* noop */}
        void halt(){/* noop */}
        bool stop(bool userInterrupt, std::size_t volume, float time){
            return userInterrupt || (m_targets.empty() && m_number_of_targets > 0) || time >= m_until_time;
        }
        float progress(std::size_t volume, float time){
            float fac_time = time / m_until_time;
            float max = fac_time;
            if(m_number_of_targets > 0){
                float fac_targets = 1.0f - (static_cast<float>(m_targets.size()) / static_cast<float>(m_number_of_targets));
                max = std::max(max, fac_targets);
            }
            return max;
        }
        MappingView& time() {return m_time_data;}
    public:
        const MappingView& time() const {return m_time_data;}
        void maxTime(float until_time){m_until_time = until_time;}
        template <class InputIter>
        void setTargets(InputIter first, InputIter last){
            m_targets.clear();
            m_targets.insert(first, last);
            m_number_of_targets = m_targets.size();
        }
        void addTarget(int64_t id){
            m_targets.insert(id);
            ++m_number_of_targets;
        }
    private:
        MappingView m_time_data;
        float m_until_time = FLT_MAX;
        Set m_targets;
        std::size_t m_number_of_targets = 0;

        template <class O, class R>
        friend class FastMarching;
    };

    template<class MV>
    class ObserverDistance {
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
        void fix(int64_t id){/* noop */}
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
            // if(userInterrupt){
            //     std::cout << "user Interrupt" << std::endl;
            // }
            // if(time >= m_until_time){
            //     std::cout << "time Interrupt" << std::endl;
            // }
            // if(m_max_distance >= m_until_distance){
            //     std::cout << "distance Interrupt" << std::endl;
            //     std::cout << "max_dist:" << m_max_distance << std::endl;
            //     std::cout << "until_dist:" << m_until_distance << std::endl;
            // }
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
    private:
        MappingView m_time_data;
        MappingView m_distance_data;

        float m_until_time = FLT_MAX;
        float m_until_distance = FLT_MAX;
        float m_max_distance = 0.0f;

        template <class O, class R>
        friend class FastMarching;
    };

    // clang-format on
}
