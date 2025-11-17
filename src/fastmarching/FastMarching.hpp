#pragma once

// #include <mclib/McDim3l.h>
// #include <hxfield/RawField.h>

// #include <unordered_map>
// #include <vector>

// #include <mutil/List.hpp>
// #include <mutil/Mapping.hpp>

#include <fastmarching/FastMarchingUtilities.hpp>
#include <fastmarching/FastMarchingObserver.hpp>
#include <utils/ProgressbarReportHelper.hpp>
#include <utils/Mapping.hpp>
#include <utils/HashMap.hpp>
#include <utils/Utils.hpp>
#include <field/RawField.hpp>


// #include <hxfastmarching/private/helpers/ProgressbarReportHelper.h>
// #include <hxfastmarching/private/helpers/FastMarchingObserver.h>

// // TODO: test which heap data structure is fastest
#include <boost/heap/fibonacci_heap.hpp>
// #include <hxhelpers/internal/HxFloatingPointRepresentation.h>

namespace fastmarching
{
    /**
     * @brief
     *
     * @tparam MappingView is from uint64_t to float (see mutil/mapping.hpp)
     *
     *
     * To use:
     * - Create FastMarching Object
     * - data()
     * - thresholds()
     * - add-/ replace-/ setStartVoxel(s)
     * - set observer conditions (e.g. observer().setTarget())
     * - march
     */
    template <class Observer = ObserverAll<mutil::HashMapMappingView<BigHashMap<uint64_t, float>>, SmallHashSet<uint64_t>>>
    class FastMarching
    {
    public:
        using MappingView = typename Observer::MappingView;

        FastMarching(progressbar::Progressbar& report)
            : m_min_threshold(std::numeric_limits<float>::lowest()), m_max_threshold(std::numeric_limits<float>::max()), m_handles(), m_heap(), m_potential(nullptr), m_max_time(0.0f), m_counter(0), m_report(report)
        {
        }

    private: // inner data shared by all versions
        void
        data_inner(RawField<float> *potential)
        {
            m_potential = potential;

            // we create a new heap with the new timeData comparator (functor can not be changed unless with a const-cast)
            auto comp = [this](const uint64_t &l, const uint64_t &r)
            {
                return m_observer.time().get(l) > m_observer.time().get(r);
            };
            m_heap = heap_type(comp);
        }

    public:
        template <class T = MappingView, std::enable_if_t<std::is_same<typename T::Ownership, mutil::Owned>::value, int> = 0>
        void
        data(RawField<float> *pot)
        {
            data_inner(pot);
        }

        template <class T = Observer, class S = MappingView, std::enable_if_t<std::is_same_v<typename T::Distance, DistanceDisabled> && std::is_same_v<typename S::Ownership, mutil::Borrowed>, int> = 0>
        void
        data(RawField<float> *pot, MappingView timeData)
        {
            data_inner(pot);
            m_observer.bind(timeData);
        }

        template <class T = Observer, class S = MappingView, std::enable_if_t<std::is_same_v<typename T::Distance, DistanceEnabled> && std::is_same_v<typename S::Ownership, mutil::Borrowed>, int> = 0>
        void
        data(RawField<float> *pot, MappingView timeData, MappingView distanceData)
        {
            data_inner(pot);
            m_observer.bind(timeData, distanceData);
        }

        /**
         * @brief Set the threshold of the easing function. If the given parameters are invalid, we use the default threshold.
         *
         * @param min
         * @param max
         *
         * @return true, if params were set. False if params are invalid and default parameters were used instead.
         */
        bool
        thresholds(float min, float max)
        {
            return thresholds(std::pair<float, float>(min, max));
        }

        bool
        thresholds(std::pair<float, float> thresh)
        {
            if (thresh.first >= thresh.second)
            {
                thresholds_unsafe(defaultThresholds());
                return false;
            }
            else
            {
                m_min_threshold = thresh.first;
                m_max_threshold = thresh.second;
                return true;
            }
        }

    private:
        void
        thresholds_unsafe(std::pair<float, float> thresh)
        {
            m_min_threshold = thresh.first;
            m_max_threshold = thresh.second;
        }

    public:
        /**
         * @brief Get the default thresholds of the input image
         *
         * If the input field is not set yet, we return (0,0)
         *
         * @return std::pair<float, float> min and max threshold value
         */
        std::pair<float, float>
        defaultThresholds() const
        {
            return defaultThresholds(m_potential);
        }

        std::pair<float, float>
        defaultThresholds(RawField<float> *potential) const
        {
            if (!potential)
            {
                return std::pair<float, float>(0.0f, 0.0f);
            }
            auto range = potential->getRange();
            return std::pair<float, float>(range.first, range.second * 1.2);
        }

        std::pair<float, float>
        thresholds()
        {
            return std::pair<float, float>(m_min_threshold, m_max_threshold);
        }

        // clang-format off
        float minThreshold() const {return m_min_threshold;}
        float maxThreshold() const {return m_max_threshold;}
        float undefinedValue() const {return INFINITY;}
        // clang-format on

        /**
         * @brief Add a seeding voxel. This has no side effect.
         *
         * @param startVoxel
         */
        void
        addStartVoxel(uint64_t startVoxel)
        {
            m_observer.seed(startVoxel);
            m_handles[startVoxel] = m_heap.push(startVoxel);
        }

        /**
         * @brief Add seeding Voxels. This has no side effect.
         *
         * @tparam Iter
         * @tparam Sentinel
         * @param start
         * @param end
         */
        template <class Iter, class Sentinel>
        void
        addStartVoxels(Iter start, Sentinel end)
        {
            while (start != end)
            {
                addStartVoxel(*start);
                ++start;
            }
        }

        /**
         * @brief Add a seeding voxel. This has no side effect.
         *
         * FastMarching only acts on voxels!
         * This is a convenient method, converting the point to the corresponding voxel!
         *
         * @param startVoxel
         */
        void
        addStartPoint(VecFloat startPoint)
        {
            addStartVoxel(static_cast<uint64_t>(m_potential->createLocation(startPoint).index()));
        }

        /**
         * @brief Add seeding Voxels. This has no side effect.
         *
         * FastMarching only acts on voxels!
         * This is a convenient method, converting the point to the corresponding voxel!
         *
         * @tparam Iter
         * @tparam Sentinel
         * @param start
         * @param end
         */
        template <class Iter, class Sentinel>
        void
        addStartPoints(Iter start, Sentinel end)
        {
            while (start != end)
            {
                addStartPoint(*start);
                ++start;
            }
        }

        /**
         * @brief Replace the seeding voxels. That is, the marching front up to this call gets deleted.
         *
         * @param startVoxel
         */
        void
        replaceStartVoxel(uint64_t startVoxel)
        {
            m_handles.clear();
            m_heap.clear();
            addStartVoxel(startVoxel);
        }

        /**
         * @brief Replace the seeding voxels. That is, the marching front up to this call gets deleted.
         *
         * @tparam Iter
         * @tparam Sentinel
         * @param start
         * @param end
         */
        template <class Iter, class Sentinel>
        void
        replaceStartVoxels(Iter start, Sentinel end)
        {
            m_handles.clear();
            m_heap.clear();
            addStartVoxels(start, end);
        }

        /**
         * @brief Replace the seeding voxels. That is, the marching front up to this call gets deleted.
         *
         * @param startVoxel
         */
        void
        replaceStartPoint(VecFloat startPoint)
        {
            m_handles.clear();
            m_heap.clear();
            addStartPoint(startPoint);
        }

        /**
         * @brief Replace the seeding voxels. That is, the marching front up to this call gets deleted.
         *
         * @tparam Iter
         * @tparam Sentinel
         * @param start
         * @param end
         */
        template <class Iter, class Sentinel>
        void
        replaceStartPoints(Iter start, Sentinel end)
        {
            m_handles.clear();
            m_heap.clear();
            addStartPoints(start, end);
        }

        /**
         * @brief Set a seeding voxel. All other seeding voxels, fronts and information on the fields gets deleted.
         *
         * @param startVoxel
         */
        void
        setStartVoxel(uint64_t startVoxel)
        {
            m_observer.reset(INFINITY);
            replaceStartVoxel(startVoxel);
        }

        /**
         * @brief Set the seeding voxels. All other seeding voxels, fronts and information on the fields gets deleted.
         *
         * @param startVoxel
         */
        template <class Iter, class Sentinel>
        void
        setStartVoxels(Iter start, Sentinel end)
        {
            m_observer.reset(INFINITY);
            replaceStartVoxels(start, end);
        }

        /**
         * @brief Set a seeding voxel. All other seeding voxels, fronts and information on the fields gets deleted.
         *
         * @param startVoxel
         */
        void
        setStartPoint(VecFloat startPoint)
        {
            m_observer.reset(INFINITY);
            replaceStartPoint(startPoint);
        }

        /**
         * @brief Set the seeding voxels. All other seeding voxels, fronts and information on the fields gets deleted.
         *
         * @param startVoxel
         */
        template <class Iter, class Sentinel>
        void
        setStartPoints(Iter start, Sentinel end)
        {
            m_observer.reset(INFINITY);
            replaceStartPoints(start, end);
        }

        /**
         * @brief Set up marching the image. Usually one should only need @c FastMarching::march
         *
         * @c FastMarching::continueMarch should be called afterwards.
         * @c FastMarching::endMarch must be called as well.
         */
        void
        beginMarch()
        {
            m_report.start("Marching...");
            m_observer.start();
            m_counter = 0;
            m_max_time = 0;
        }

        /**
         * @brief If marching stopped, the stopping criteria were changed and one should continue marching, this method can be used.
         *
         * @c FastMarching::continueMarch should be called before calling this method.
         * @c FastMarching::endMarch must be called afterwards.
         */
        void
        continueMarch()
        {
            // std::cout << "DEBUG: start march!" << std::endl;
            if (!m_observer.proceed())
            {
                std::cout << "observer did not want to proceed!" << std::endl;
                return;
            }
            auto dims = m_potential->getDims();
            uint64_t node;

            // START ALGORITHM
            while (!m_heap.empty())
            {
                
                node = m_heap.top();
                m_heap.pop();
                m_handles.erase(node);
                m_observer.fix(node);
                
                // go through all neighbors and look for ones which have to be updated.
                // TODO: field locations of heap are from a specific field (currently m_potential)
                auto neighbors = m_potential->gridNeighbors6(RawField<float>::FieldLocation(node));
                // futil::gridNeighbors6(dims, node, neighbors);
                
                // PERF: we know that node has to be set in data -> use unsafe methods
                m_max_time = m_observer.time().get_unsafe(node);
                // // debug:
                // auto gridPos = Lattice::gridLocationFromCIndex(node, dims);
                // if(gridPos == VecSize(0,6,5) || gridPos == VecSize(1,6,5) || gridPos == VecSize(0,5,5)){
                //     std::cout << "DEBUG:  node: " << Lattice::gridLocationFromCIndex(node, dims) << " and max time: " << m_max_time << ", euclidiean: " << m_observer.distance().get_unsafe(node) << std::endl;
                // }
                
                for (auto neighbor : neighbors)
                {
                    if (neighbor.empty())
                    {
                        continue;
                    }
                    float time_neighbor = m_observer.time().get(neighbor.index());
                    // // debug
                    // if(gridPos == VecSize(0,5,5) && neighbor == m_potential->createLocation(VecSize(0,6,5))){
                    //     std::cout << "index: " << neighbor.index() << std::endl;
                    //     auto myNeighbors = m_potential->gridNeighbors6(neighbor);
                    //     for(auto my_neighbor : myNeighbors){
                    //         std::cout << "DEBUG:  neighbor: " << my_neighbor.location(dims) << " with time: " << m_observer.time().get(my_neighbor.index()) << ", euclidiean: " << m_observer.distance().get_optional(my_neighbor.index()).value_or(INFINITY) << std::endl;
                    //     }
                    // }
                    if (time_neighbor <= m_max_time)
                    {
                        // neighbor is already fixed (or value can't be smaller), so we do not have to calculate the cost
                        //(we assume cost calculation to be more expensive than this check)
                        continue;
                    }

                    float new_cost = fastmarching::eikonal_view(*m_potential, neighbor, m_min_threshold, m_max_threshold, m_observer.time());
                    // std::cout << "DEBUG:  new cost: " << new_cost << std::endl;
                    // TODO: debug assert that new_cost > m_max_time

                    if (new_cost < time_neighbor)
                    {
                        //TODO: only use the voxels we alreay fixed -> where time <= new_cost
                        m_observer.update(neighbor.index(), new_cost, dims);
                        // Check of time_neighbor == INFINITY possible instead, but only if all seedpoints of one front are caluclated at once (instead of added later, which should be possible currently).
                        // That is: we allow the possibility of addSeedpoint() -> march() -> addSeedpoint() -> march() again with the expected values
                        if (m_handles.find(neighbor.index()) == m_handles.end())
                        {
                            // std::cout << "DEBUG:  push new element onto heap" << std::endl;
                            m_handles[neighbor.index()] = m_heap.push(neighbor.index());
                            m_observer.candidate(neighbor.index());
                        }
                        else
                        {
                            // std::cout << "DEBUG:  element updated on heap" << std::endl;
                            // we decrease the value. As we use a min-heap instead of max-heap, we "increase" the value from the perspective of the data structure.
                            m_heap.increase(m_handles.at(neighbor.index()));
                        }
                    }
                }

                ++m_counter;

                // should we stop?
                if (m_observer.stop(m_report.stop(), m_counter, m_max_time)){
                    // std::cout << "DEBUG:  observer says we should stop" << std::endl;
                    return;
                }

                // progress bar calculations every 65536 cycles
                if ((m_counter & 0xFFFF) == 0xFFFF)
                {
                    m_report.update(m_observer.progress(m_counter, m_max_time));
                }
            }
            m_observer.halt();
        }

        /**
         * @brief Must be called after calling @c FastMarching::beginMarch
         *
         * @return float time of last fixed node.
         */
        float
        endMarch()
        {
            m_observer.end();
            m_report.end();
            return m_max_time;
        }

        /**
         * @brief Marches through the image.
         *
         * @return float time of last fixed node.
         */
        float
        march()
        {
            beginMarch();
            continueMarch();
            return endMarch();
        }

        void maxTime(float time)
        {
            m_observer.maxTime(time);
        }

        template <class T = Observer, typename std::enable_if<std::is_same<typename T::Distance, DistanceEnabled>::value, int>::type = 0>
        void maxDistance(float distance)
        {
            m_observer.maxDistance(distance);
        }

        /**
         * @brief Returns the same value as the march function.
         */
        float
        maxTimeMarched() const
        {
            return m_max_time;
        }

        /**
         * @brief Returns the max distance reached.
         */
        template <class T = Observer, typename std::enable_if<std::is_same<typename T::Distance, DistanceEnabled>::value, int>::type = 0>
        float
        maxDistanceMarched() const
        {
            return m_observer.maxDistanceMarched();
        }

        std::size_t
        numberOfIterations() const
        {
            return m_counter;
        }

        /**
         * @brief Returns a list of all voxels seen but not fixed.
         *
         * @return std::vector<uint64_t>
         */
        std::vector<uint64_t>
        borderVoxels() const
        {
            auto vec = std::vector<uint64_t>();
            for (auto pair : m_handles)
            {
                vec.push_back(pair.first);
            }
            return vec;
        }

        /**
         * @brief Returns time field as a mapping view
         *
         */
        const typename Observer::MappingView &
        time() const
        {
            return m_observer.time();
        }

        /**
         * @brief Copies the time field to the data pointer
         *
         * This is only a helper method for fast debugging.
         * Often Fast marching won't touch all voxels of the image.
         * Such, this method is slow, as it looks at all possible voxels.
         *
         * Capping Values is helpful as often visualizations do not care about the set undefined value.
         *
         */
        float
        time(float *dataPtr, std::size_t size, bool capValues = true) const
        {
            float min = std::numeric_limits<float>::max();
            float max = 0.0f;
            float undefined = undefinedValue();

            for (std::size_t i = 0; i < size; ++i)
            {
                float val = m_observer.time().get(i);
                dataPtr[i] = val;
                if (val != undefined)
                {
                    min = std::min(min, val);
                    max = std::max(max, val);
                }
            }

            return capValues ? mutil::capValues(dataPtr, size, max) : undefined;
        }

        /**
         * @brief Returns distance as a mapping view
         */
        template <class T = Observer, typename std::enable_if<std::is_same<typename T::Distance, DistanceEnabled>::value, int>::type = 0>
        const typename Observer::MappingView &
        distance() const
        {
            return m_observer.distance();
        }

        /**
         * @brief Copies the distance field to the data pointer
         *
         * This is only a helper method for fast debugging.
         * Often Fast marching won't touch all voxels of the image.
         * Such, this method is slow, as it looks at all possible voxels.
         *
         */
        template <class T = Observer, typename std::enable_if<std::is_same<typename T::Distance, DistanceEnabled>::value, int>::type = 0>
        float
        distance(float *dataPtr, std::size_t size, bool capValues = true) const
        {
            float min = std::numeric_limits<float>::max();
            float max = 0.0f;
            float undefined = undefinedValue();
            for (std::size_t i = 0; i < size; ++i)
            {
                float val = m_observer.distance().get(i);
                dataPtr[i] = val;
                if (val != undefined)
                {
                    min = std::min(min, val);
                    max = std::max(max, val);
                }
            }

            return capValues ? mutil::capValues(dataPtr, size, max) : undefined;
        }

        /**
         * @brief Return the internally used observer.
         *
         * @return const Observer&
         */
        const Observer &
        observer() const
        {
            return m_observer;
        }

        /**
         * @brief Return the internally used observer.
         *
         * @return Observer&
         */
        Observer &
        observer()
        {
            return m_observer;
        }

        float
        speed(uint64_t voxel)
        {
            float pot = m_potential->get(RawField<float>::Location(voxel)).value();
            return fastmarching::calculate_speed(pot, m_min_threshold, m_max_threshold);
        }

    protected:
        float m_min_threshold;
        float m_max_threshold;

        // We use std::function<>, but this incurs a performance cost (indirection). CompareFunctors are faster.
        // However boost::heap does not allow to change CompareFunctors -> const_cast from const to non-const necessary (see commit b7fad10265ec7928f4ed635c6a7a9fb203ea8cfc)
        // Best Performance Option would be creating own heap (or using heap without being able to change nodes and just mutliply nodes by at most 6)
        typedef boost::heap::fibonacci_heap<uint64_t, boost::heap::compare<std::function<bool(const uint64_t &, const uint64_t &)>>> heap_type;
        std::unordered_map<uint64_t, heap_type::handle_type> m_handles;
        heap_type m_heap;

        RawField<float> *m_potential;
        Observer m_observer;

        // values generated by march
        float m_max_time;
        std::size_t m_counter;

        // Progressbar Report
        progressbar::Progressbar& m_report;
    };

} // namespace fastmarching
