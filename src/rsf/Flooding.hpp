#pragma once

#include <unordered_set>
#include <vector>

#include <utils/Vec.hpp>
#include <utils/Dims.hpp>
#include "Faces.hpp"

// TODO: test which heap data structure is fastest
// #include <boost/heap/fibonacci_heap.hpp>
#include <boost/heap/priority_queue.hpp>

namespace klenert
{
    // /**
    //  * @brief Flooding algorithm which returns the exit, called sink. Inverse as it does not flood low values, but high values.
    //  *
    //  * Given a scalar field and a starting point, the field is flooded (it prefers high values) until a voxel with a value at least as high as the given sink is encountered.
    //  *
    //  * @param value
    //  * @return std::pair<int64_t, int64_t> The voxel on the border and the found sink voxel.
    //  */
    // template <class T>
    // std::pair<int64_t, int64_t> inverseFloodingExit(T* value, McDim3l dims, int64_t start, T border)
    // {
    //     auto compare = [value](int64_t left, int64_t right) { value[left] < value[right] };
    //     typedef std::priority_queue<int64_t, std::vector<int64_t>, boost::heap::compare<decltype(compare)> > heap_type;
    //     heap_type heap = heap_type(compare);

    //     auto seen = std::unordered_set<int64_t>();

    //     heap.push(start);
    //     seen.insert(start);
    //     int64_t node = heap.top();
    //     int64_t neighbors[6];

    //     while (!heap.empty() && value[node] <= border)
    //     {
    //         node = heap.top();
    //         heap.pop();
    //         // get neighbors of node
    //         auto counter = futil::gridNeighbors6(dims, node, &neighbors[0]);
    //         for (int i = 0; i < counter; ++i)
    //         {
    //             int64_t neighbor = neighbors[i];
    //             if (value[neighbor] > border)
    //             {
    //                 return std::pair<int64_t, int64_t>(node, neighbor);
    //             }
    //             if (seen.find(neighbor) == seen.end())
    //             {
    //                 heap.push(neighbor);
    //                 seen.insert(neighbor);
    //             }
    //         }
    //     }
    //     return std::pair<int64_t, int64_t>(node, node);
    // }

    /**
     * @brief WaterStream is flooding of one (possible mutliple) source point until some condition is met.
     *
     * @tparam T
     * @tparam Connectivity
     * @tparam MaxNeighbors
     * @tparam Breaker
     * @tparam Compare
     */
    template <class T, class Connectivity, std::size_t MaxNeighbors, class Breaker, class Compare>
    class WaterStream
    {
    public:
        WaterStream(Connectivity connection, Breaker breaker, Compare compare)
            : m_heap(compare)
            , m_seen()
            , m_connection(connection)
            , m_break(breaker)
        {}

        std::pair<T, T> stream(T start)
        {
            m_heap.clear();
            m_seen.clear();
            m_heap.push(start);
            m_seen.insert(start);
            T node = m_heap.top();
            T neighbors[MaxNeighbors];

            while (!m_heap.empty())
            {
                node = m_heap.top();
                m_heap.pop();
                // get neighbors of node
                auto counter = m_connection(node, &neighbors[0]);
                for (std::size_t i = 0; i < counter; ++i)
                {
                    T neighbor = neighbors[i];
                    if (m_break(neighbor))
                    {
                        return std::pair<T, T>(node, neighbor);
                    }
                    if (m_seen.find(neighbor) == m_seen.end())
                    {
                        m_heap.push(neighbor);
                        m_seen.insert(neighbor);
                    }
                }
            }
            return std::pair<T, T>(node, node);
        }

        const std::unordered_set<T>& touched() const
        {
            return m_seen;
        }

        void clear()
        {
            m_heap.clear();
            m_seen.clear();
        }

    protected:
        typedef boost::heap::priority_queue<T, boost::heap::compare<Compare> > heap_type;
        heap_type m_heap;

        std::unordered_set<T> m_seen;

        Connectivity m_connection;
        Breaker m_break;
    };

    /**
     * @brief Specialized WaterStream is flooding of one (possible mutliple) source point until some condition is met or we are outside the image.
     *
     * @tparam MaxNeighbors
     * @tparam Breaker
     * @tparam Compare
     * 
     * TODO: we should change int64_t to std::size_t (or uint64_t)
     */
    template <class Breaker, class Compare>
    class GridWaterStream
    {
    public:
        GridWaterStream(Breaker breaker, Compare compare, Dims dims)
            : m_heap(compare)
            , m_seen()
            , m_taken()
            , m_break(breaker)
            , m_dims(dims)
        {}

        Face stream(int64_t start)
        {
            m_heap.clear();
            m_seen.clear();
            m_taken.clear();
            m_heap.push(start);
            m_seen.insert(start);

            // heap should never be empty, as the break if we find the image border
            for (;;)
            {
                int64_t node_id = m_heap.top();
                m_taken.insert(node_id);

                VecSize node = Lattice::gridLocationFromCIndex(node_id, m_dims);
                m_heap.pop();

                if (node[0] == 0)
                    return Face(node_id, Direction::LEFT);
                if (node[0] == m_dims[0] - 1)
                    return Face(node_id, Direction::RIGHT);
                if (node[1] == 0)
                    return Face(node_id, Direction::DOWN);
                if (node[1] == m_dims[1] - 1)
                    return Face(node_id, Direction::UP);
                if (node[2] == 0)
                    return Face(node_id, Direction::BACKWARD);
                if (node[2] == m_dims[2] - 1)
                    return Face(node_id, Direction::FORWARD);

                std::vector<int64_t> neighbors = {
                    node_id - 1,
                    node_id + 1,
                    node_id - static_cast<int64_t>(m_dims[0]),
                    node_id + static_cast<int64_t>(m_dims[0]),
                    node_id - static_cast<int64_t>(m_dims[0] * m_dims[1]),
                    node_id + static_cast<int64_t>(m_dims[0] * m_dims[1]),
                };

                for (auto neighbor : neighbors)
                {
                    if (m_break(neighbor))
                    {
                        return Face(node_id, neighbor, m_dims);
                    }
                    if (m_seen.find(neighbor) == m_seen.end())
                    {
                        m_seen.insert(neighbor);
                        m_heap.push(neighbor);
                    }
                }
            }
        }

        const std::unordered_set<int64_t>& touched() const
        {
            return m_seen;
        }

        const std::unordered_set<int64_t>& taken() const
        {
            return m_taken;
        }

        void clear()
        {
            m_heap.clear();
            m_seen.clear();
            m_taken.clear();
        }

    protected:
        typedef boost::heap::priority_queue<int64_t, boost::heap::compare<Compare> > heap_type;
        heap_type m_heap;

        std::unordered_set<int64_t> m_seen;
        std::unordered_set<int64_t> m_taken;

        Breaker m_break;
        Dims m_dims;
    };

    /**
     * @brief Flooding is has multiple source points and "colors" each region, depending on the source point.
     *
     * @tparam T
     * @tparam Connectivity
     * @tparam MaxNeighbors
     * @tparam Compare
     */
    template <class T, class Connectivity, std::size_t MaxNeighbors, class Compare>
    class Flooding
    {
    public:
        Flooding(Connectivity connection, Compare compare)
            : m_connection(connection)
            , m_compare(compare)
            , m_heap(compare)
            , m_colors()
        {}

        /**
         * Flood from local minimas of the startpoints. 
         * If startpoints with different labelings share the same locale minima, we return true (can be seen as an error).
         */
        bool streamflood(std::vector<T> startPoints, std::vector<std::size_t> labeling)
        {
            m_heap.clear();
            m_colors.clear();

            T node;
            T neighbors[MaxNeighbors];
            std::size_t counter = 0;

            // we first need to stream down to local minima
            for (std::size_t i = 0; i < startPoints.size(); ++i)
            {
                node = startPoints[i];
                // m_colors[node] = labeling[i];
                for(;;){
                    // check if we at minima
                    T minima = node;
                    counter = m_connection(node, &neighbors[0]);
                    for (std::size_t i = 0; i < counter; ++i)
                    {
                        T neighbor = neighbors[i];
                        if(m_compare(minima, neighbor)){
                            minima = neighbor;
                        }
                    }
                    if(minima == node){
                        break;
                    }else{
                        // we found a node to stream down
                        node = minima;
                    }
                }

                // add locale minima to heap
                if(m_colors.find(node) != m_colors.end()){
                    return true;
                }else{
                    m_heap.push(node);
                    m_colors[node] = labeling[i];
                }
            }

            inner_flood();
            return false;
        }

        bool streamflood(std::vector<T> startPoints)
        {
            auto numbers = std::vector<std::size_t>();
            numbers.reserve(startPoints.size());
            for (std::size_t i = 1; i <= startPoints.size(); ++i)
            {
                numbers.push_back(i);
            }
            return streamflood(startPoints, numbers);
        }

        void flood(std::vector<T> startPoints, std::vector<std::size_t> labeling)
        {
            m_heap.clear();
            m_colors.clear();

            for (std::size_t i = 0; i < startPoints.size(); ++i)
            {
                auto startPoint = startPoints[i];
                m_heap.push(startPoint);
                m_colors[startPoint] = labeling[i];
            }

            inner_flood();
        }

        void flood(std::vector<T> startPoints)
        {
            m_heap.clear();
            m_colors.clear();

            for (std::size_t i = 0; i < startPoints.size(); ++i)
            {
                auto startPoint = startPoints[i];
                m_heap.push(startPoint);
                m_colors[startPoint] = i;
            }

            inner_flood();
        }

        const std::unordered_map<T, std::size_t>& regions() const
        {
            return m_colors;
        }

        void clear()
        {
            m_heap.clear();
            m_colors.clear();
        }

    protected:
        void inner_flood()
        {
            T node;
            T neighbors[MaxNeighbors];
            std::size_t counter = 0;

            while (!m_heap.empty())
            {
                node = m_heap.top();
                m_heap.pop();
                // get neighbors of node
                counter = m_connection(node, &neighbors[0]);
                for (std::size_t i = 0; i < counter; ++i)
                {
                    T neighbor = neighbors[i];
                    if (m_colors.find(neighbor) == m_colors.end())
                    {
                        m_heap.push(neighbor);
                        m_colors[neighbor] = m_colors[node];
                    }
                }
            }
        }

        Connectivity m_connection;
        Compare m_compare;
        typedef boost::heap::priority_queue<T, boost::heap::compare<Compare> > heap_type;
        heap_type m_heap;
        std::unordered_map<T, std::size_t> m_colors;
    };

} // namespace klenert
