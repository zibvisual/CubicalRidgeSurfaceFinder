#pragma once

#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include <optional>

#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t

namespace klenert
{
    template <typename E>
    class SparseGraph
    {
    public:
        class EdgeIterator
        {
        public:
            using value_type = std::array<std::size_t, 2>;

            EdgeIterator(const std::vector<std::unordered_map<std::size_t, E>>& vec)
                : m_source(vec)
                , m_source_current(0)
                , m_target_current(m_source[m_source_current].cbegin())
                , m_target_end(m_source[m_source_current].cend())
            {}

            std::optional<value_type>
            next()
            {
                if (m_source_current >= m_source.size())
                {
                    return std::optional<value_type>();
                }

                while (m_target_current == m_target_end)
                {
                    ++m_source_current;
                    if (m_source_current >= m_source.size())
                    {
                        return std::optional<value_type>();
                    }
                    m_target_current = m_source[m_source_current].cbegin();
                    m_target_end = m_source[m_source_current].cend();
                }

                std::size_t from = m_source_current;
                std::size_t to = m_target_current->first;
                ++m_target_current;

                return std::optional<value_type>({ from, to });
            }

            // ADD_RANGE_TO_RUST_ITERATOR(EdgeIterator)

        protected:
            const std::vector<std::unordered_map<std::size_t, E>>& m_source;
            std::size_t m_source_current;
            typename std::unordered_map<std::size_t, E>::const_iterator m_target_current;
            typename std::unordered_map<std::size_t, E>::const_iterator m_target_end;
        };

        class NeighborIterator
        {
        public:
            using value_type = std::pair<std::size_t, E>;

            NeighborIterator(const std::unordered_map<std::size_t, E>& map)
                : map(map)
            {}

            typename std::unordered_map<std::size_t, E>::const_iterator
            begin()
            {
                return map.cbegin();
            }

            typename std::unordered_map<std::size_t, E>::const_iterator
            end()
            {
                return map.cend();
            }

        protected:
            const std::unordered_map<std::size_t, E>& map;
        };

        using value_type = std::size_t;

        SparseGraph()
        {}

        SparseGraph(std::size_t nVertices)
            : m_edges()
        {
            m_edges.reserve(nVertices);
        }

        // Returns rust iterator
        EdgeIterator
        edges() const
        {
            return EdgeIterator(m_edges);
        }

        NeighborIterator
        neighbors(uint64_t id) const
        {
            return NeighborIterator(m_edges.at(id));
        }

        void
        addEdge(uint64_t first, uint64_t second, E payload)
        {
            m_edges[first].insert({ second, payload });
        }

        void
        removeEdge(uint64_t first, uint64_t second)
        {
            m_edges[first].erase(second);
        }

        /**
         * Filte all edges with a predicate.
         */
        template <typename BinaryPredicate>
        void
        filter(BinaryPredicate pred)
        {

            for (auto out = m_edges.begin(); out != m_edges.end(); ++out)
            {
                auto& edge = out->second;
                for (auto it = edge.begin(); it != edge.end();)
                {
                    if (pred(out->first, it->first))
                    {
                        ++it;
                    }
                    else
                    {
                        it = edge.erase(it);
                    }
                }
            }
        }

        void
        addNode(uint64_t id)
        {
            m_edges.insert({id, std::unordered_map<uint64_t, E>()});
        }

        void
        removeNode(std::size_t id)
        {
            m_edges.erase(id);
            for(auto& edge : m_edges){
                edge.second.erase(id);
            }
        }

        // void
        // replaceNodes(std::size_t left, std::size_t right)
        // {
        //     m_edges[left] = m_edges[right];
        //     for (auto& map : m_edges)
        //     {
        //         auto it = map.find(right);
        //         if (it != map.end())
        //         {
        //             // erase first to make sure that the iterator is still valid
        //             auto payload = it->second;
        //             map.erase(it);
        //             map.insert({ left, payload });
        //         }
        //     }
        // }

        void
        clear()
        {
            m_edges.clear();
        }

        std::size_t
        numberOfNodes() const
        {
            return m_edges.size();
        }

        bool
        hasNode(std::size_t node) const
        {
            return m_edges.contains(node);
        }

        bool
        hasEdge(std::size_t from, std::size_t to) const
        {
            return m_edges.contains(from) && m_edges[from].contains(to);
        }

        std::optional<E>
        getEdge(std::size_t from, std::size_t to) const
        {
            if (!hasNode(from))
            {
                return std::optional<E>();
            }
            auto it = m_edges.at(from).find(to);
            if (it == m_edges.at(from).end())
            {
                return std::optional<E>();
            }
            return it->second;
        }

        E
        getEdgeUnsafe(std::size_t from, std::size_t to) const
        {
            return m_edges.at(from).at(to);
        }

    protected:
        std::unordered_map<uint64_t,std::unordered_map<uint64_t, E>> m_edges;
    };

} // namespace klenert
