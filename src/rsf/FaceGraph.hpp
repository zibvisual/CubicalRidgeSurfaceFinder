#pragma once

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <queue>

#include <iterator> // For std::forward_iterator_tag
#include <cstddef>  // For std::ptrdiff_t

#include <utils/Dims.hpp>
#include "Faces.hpp"


class FaceGraph {
    public:
        class NodeIterator {
            using iterator_category = std::forward_iterator_tag;
            using difference_type = std::ptrdiff_t;
            using value_type = int64_t;

            public:
                NodeIterator(std::unordered_map<int64_t, std::array<int64_t, 4> >::const_iterator pointer)
                : m_current(pointer)
                {}

                int64_t operator*() const { return m_current->first; }
                int64_t operator->() { return m_current->first; }
                NodeIterator& operator++() { ++m_current; return *this; }
                NodeIterator operator++(int) {NodeIterator tmp = *this; ++(*this); return tmp;}

                friend bool operator== (const NodeIterator& a, const NodeIterator& b) { return a.m_current == b.m_current; }
                friend bool operator!= (const NodeIterator& a, const NodeIterator& b) { return a.m_current != b.m_current; }

                IndexToFaceConverter<NodeIterator> convert() { return IndexToFaceConverter<NodeIterator>(*this); }

            protected:
                typename std::unordered_map<int64_t, std::array<int64_t, 4> >::const_iterator m_current;
        };

        class BorderIterator {
            using iterator_category = std::forward_iterator_tag;
            using difference_type = std::ptrdiff_t;
            using value_type = int64_t;

            public:
                BorderIterator(std::unordered_set<int64_t>::const_iterator pointer)
                : m_current(pointer)
                {}

                int64_t operator*() const { return *m_current; }
                int64_t operator->() { return *m_current; }
                BorderIterator& operator++() { ++m_current; return *this; }
                BorderIterator operator++(int) {BorderIterator tmp = *this; ++(*this); return tmp;}

                friend bool operator== (const BorderIterator& a, const BorderIterator& b) { return a.m_current == b.m_current; }
                friend bool operator!= (const BorderIterator& a, const BorderIterator& b) { return a.m_current != b.m_current; }

            protected:
                typename std::unordered_set<int64_t>::const_iterator m_current;
        };


        FaceGraph() : m_edges(), m_border(), none(std::numeric_limits<int64_t>::max()) {}

        FaceGraph(std::size_t size) : m_edges(size), m_border(), none(std::numeric_limits<int64_t>::max()){}

        template<typename Set>
        static FaceGraph createWithID(Set& set, Dims dims){
            auto graph = FaceGraph(set.size());
            // Go through the set
            for(auto faceId : set){
                graph.m_edges[faceId] = {graph.none, graph.none, graph.none, graph.none};

                // for each face, get all candidates
                auto face = Face::fromIndex(faceId);
                auto neighbors = face.neighbors(set, dims, graph.m_edges[faceId].data());

                if(neighbors < 4){
                    graph.m_border.insert(faceId);
                }
            }
            return graph;
        }

        template<typename Set>
        static FaceGraph createWithFaces(Set& set, Dims dims){
            auto graph = FaceGraph(set.size());
            // Go through the set
            for(auto face : set){
                auto faceId = face.toIndex();
                graph.m_edges[faceId] = {graph.none, graph.none, graph.none, graph.none};

                // for each face, get all candidates
                auto neighbors = face.neighbors(set, dims, graph.m_edges[faceId].data());

                if(neighbors < 4){
                    graph.m_border.insert(faceId);
                }
            }
            return graph;
        }

        NodeIterator begin() { return nodesBegin(); }
        NodeIterator end() { return nodesEnd(); }

        NodeIterator nodesBegin() { return NodeIterator(m_edges.cbegin()); }
        NodeIterator nodesEnd() { return NodeIterator(m_edges.cend()); }

        BorderIterator borderBegin() { return BorderIterator(m_border.cbegin()); }
        BorderIterator borderEnd() { return BorderIterator(m_border.cend()); }


        // METHODS NECESSARY
    protected:
        std::unordered_map<int64_t, std::array<int64_t, 4>> m_edges;
        std::unordered_set<int64_t> m_border;
        int64_t none;
};