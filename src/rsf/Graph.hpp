#pragma once

#include <utils/IteratorEndSentinel.hpp>

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
            using value_type = std::array<uint64_t, 2>;
            using iterator_category = std::forward_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using pointer           = value_type;
            using reference         = value_type;

            using inner = std::unordered_map<uint64_t, E>;
            using outer = std::unordered_map<uint64_t, inner>;

            EdgeIterator(const outer& data)
                : m_source_current(data.cbegin())
                , m_source_end(data.cend())
                , m_target_current()
                , m_target_end()
            {
                if(m_source_current != m_source_end){
                    m_target_current = m_source_current->second.cbegin();
                    m_target_end = m_source_current->second.cend();
                    clip();
                }
            }

            EdgeIterator(const outer::const_iterator& outer_begin, const outer::const_iterator& outer_end, const inner::const_iterator& inner_begin, const inner::const_iterator& inner_end)
            : m_source_current(outer_begin)
            , m_source_end(outer_end)
            , m_target_current(inner_begin)
            , m_target_end(inner_end)
            {
                if(m_source_current != m_source_end){
                    clip();
                }
            }

            /**
             * We want to guarantee that current() always returns a valid edge as long as the iterator returns true (e.g. is not finished).
             * To achieve this, we must check for invalid states (m_target_current != m_target_end)
             */
            void clip() {
                while (m_target_current == m_target_end)
                {
                    ++m_source_current;
                    if (m_source_current == m_source_end)
                    {
                        return;
                    }
                    m_target_current = m_source_current->second.cbegin();
                    m_target_end = m_source_current->second.cend();
                }
            }

            void advance()
            {
                // only advance if we are not finished yet
                if (m_source_current == m_source_end)
                {
                    return;
                }

                ++m_target_current;

                // move to next valid m_target_current
                clip();
            }

            bool finished() const
            {
                return m_source_current == m_source_end;
            }

            /** UB if called without checking if iterator is already finished */
            value_type current() const
            {
                uint64_t from = m_source_current->first;
                uint64_t to = m_target_current->first;
                return { from, to }; 
            }

            std::optional<value_type>
            next()
            {
                if (finished())
                {
                    return std::optional<value_type>();
                }

                value_type val = current();
                advance();

                return std::optional<value_type>(val);
            }

            operator bool() const
            {
                return !finished();
            }

            reference operator*() const {
                return current();
            }
            pointer operator->() { 
                return current();
            }

            // Prefix increment
            EdgeIterator& operator++() { advance(); return *this; }  

            // Postfix increment
            EdgeIterator operator++(int) 
            { auto tmp = *this; ++(*this); return tmp; }

            EdgeIterator
            begin() const
            {
                return EdgeIterator(m_source_current, m_source_end, m_target_current, m_target_end);
            }
    
            iter::IteratorEndSentinel
            end() const
            {
                return {};
            }

        protected:
            typename outer::const_iterator m_source_current;
            typename outer::const_iterator m_source_end;
            typename inner::const_iterator m_target_current;
            typename inner::const_iterator m_target_end;
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

        std::size_t
        numberOfEdges() const
        {
            auto sum = 0;
            for(auto list : m_edges){
                sum += list->second.size();
            }
            return sum;
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
