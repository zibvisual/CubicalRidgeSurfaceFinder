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
    class OrderedSparseGraph
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

        OrderedSparseGraph()
        {}

        OrderedSparseGraph(std::size_t nVertices)
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
        neighbors(std::size_t i) const
        {
            return NeighborIterator(m_edges[i]);
        }

        void
        addEdge(std::size_t first, std::size_t second, E payload)
        {
            m_edges[first].insert({ second, payload });
        }

        void
        removeEdge(std::size_t first, std::size_t second)
        {
            m_edges[first].erase(second);
        }

        template <typename BinaryPredicate>
        void
        filter(BinaryPredicate pred)
        {
            for (std::size_t i = 0; i < m_edges.size(); ++i)
            {
                for (auto it = m_edges[i].begin(); it != m_edges[i].end();)
                {
                    if (pred(i, it->first))
                    {
                        ++it;
                    }
                    else
                    {
                        it = m_edges[i].erase(it);
                    }
                }
            }
        }

        void
        addNode(std::size_t val)
        {
            if (val >= m_edges.size())
            {
                m_edges.resize(val);
            }
        }

        void
        removeNode(std::size_t val)
        {
            m_edges.erase(m_edges.begin() + val);
            // every node with a value bigger than val must be reduced by 1. every node with value of val must be deleted.
            // the easiest way is to create new maps
            for (std::size_t i = 0; i < m_edges.size(); ++i)
            {
                auto& map = m_edges[i];
                auto new_map = std::unordered_map<std::size_t, E>();
                new_map.reserve(map.size());
                for (auto iter = map.begin(); iter != map.end(); ++iter)
                {
                    auto key = iter->first;
                    if (key == val)
                    {
                        continue;
                    }
                    else if (key > val)
                    {
                        --key;
                    }
                    new_map.insert({ key, iter->second });
                }
                m_edges[i] = new_map;
            }
        }

        void
        clearNode(std::size_t val)
        {
            m_edges[val].clear();
            for (auto& map : m_edges)
            {
                map.erase(val);
            }
        }

        void
        replaceNodes(std::size_t left, std::size_t right)
        {
            m_edges[left] = m_edges[right];
            for (auto& map : m_edges)
            {
                auto it = map.find(right);
                if (it != map.end())
                {
                    // erase first to make sure that the iterator is still valid
                    auto payload = it->second;
                    map.erase(it);
                    map.insert({ left, payload });
                }
            }
        }

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
            return node < m_edges.size();
        }

        bool
        hasEdge(std::size_t from, std::size_t to) const
        {
            if (!hasNode(from))
            {
                return false;
            }
            return m_edges[from].find(to) != m_edges[from].end();
        }

        std::optional<E>
        getEdge(std::size_t from, std::size_t to) const
        {
            if (!hasNode(from))
            {
                return std::optional<E>();
            }
            auto it = m_edges[from].find(to);
            if (it == m_edges[from].end())
            {
                return std::optional<E>();
            }
            return m_edges[from].at(to);
            // return *it;
        }

        E
        getEdgeUnsafe(std::size_t from, std::size_t to) const
        {
            return m_edges[from].at(to);
        }

    protected:
        std::vector<std::unordered_map<std::size_t, E>> m_edges;
    };

} // namespace klenert
