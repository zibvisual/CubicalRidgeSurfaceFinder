#pragma once

#include <utils/Vec.hpp>
#include <io/tgf.hpp>
#include <io/sls.hpp>

#include <vector>
#include <unordered_map>

class SpatialGraph {
public:

    SpatialGraph() {}
    // TODO: load io/tgf

    void save_as_tgf(std::string output_path) const {
        write_graph_data(output_path, toRawGraphData());
    }

    void save_as_lineset(std::string output_path) const {
        write_line_data(output_path, toLineSet());
    }
    
    graph_data_t toRawGraphData() const {
        std::unordered_map<std::string, VecFloat> points;
        std::vector<std::array<std::string, 2>> edges;

        points.reserve(m_points.size());
        edges.reserve(m_edges.size());

        auto mapping = names();
        for(std::size_t i = 0; i < mapping.size(); ++i){
            points[mapping[i]] = m_points[i];
        }
        for(auto edge : m_edges){
            edges.push_back({mapping[edge[0]], mapping[edge[1]]});
        }

        return graph_data_t{points, edges};
    }

    LineSet<std::string> toLineSet() const {
        
        auto builder = LineSetBuilder<std::string>();
        builder.reserve_lines(m_edges.size());

        auto mapping = names();

        for(auto edge : m_edges){
            builder.add_point(m_points[edge[0]], mapping[edge[0]]);
            builder.add_point(m_points[edge[1]], mapping[edge[1]]);
            builder.push_line();
        }

        return builder.build();
    }

    void add_point(std::string name, VecFloat point){
        m_ids.insert({name, m_points.size()});
        m_points.push_back(point);
    }

    void add_point(uint64_t id, VecFloat point){
        add_point(std::to_string(id), point);
    }

    void add_edge(std::string from, std::string to){
        m_edges.push_back({m_ids[from], m_ids[to]});
    }

    void add_edge(uint64_t from, uint64_t to){
        add_edge(std::to_string(from), std::to_string(to));
    }

protected:
    std::vector<std::string> names() const {
        std::vector<std::string> names;
        names.resize(m_ids.size());
        for(auto map : m_ids){
            names[map.second] = map.first;
        }
        return names;
    }

    std::vector<VecFloat> m_points;
    std::vector<std::array<std::size_t, 2>> m_edges;
    std::unordered_map<std::string, std::size_t> m_ids;
};