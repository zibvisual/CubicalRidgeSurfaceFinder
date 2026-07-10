#include <queue>

// template <class T, class S, class C>
// void clearpq(std::priority_queue<T, S, C>& q) {
//     struct HackedQueue : private std::priority_queue<T, S, C> {
//         static S& Container(std::priority_queue<T, S, C>& q) {
//             return q.*&HackedQueue::c;
//         }
//     };
//     HackedQueue::Container(q).clear();
// }

template <class Index>
using Distances = std::vector<std::tuple<VecFloat, float>>;

template <class Adjacency, class Index>
class Dijkstra {
public:
    Dijkstra(){
        m_dist.reserve(128);
        m_result.reserve(128);
    }

    const Distances<Index>& run(const std::vector<VecFloat>& points, const std::vector<Adjacency>& adj, const Index start, float radius){
        m_dist.clear();
        // clearpq(m_queue);

        // startnode with distance 0
        m_queue.push(QE(0.0f, start));
        // djkstra
        while (!m_queue.empty())
        {
            QE top = m_queue.top();
            m_queue.pop();

            const float cur_dist = top.first;
            const Index cur_vert = top.second;
            // if a (better) path was found, ignore
            if (m_dist.find(cur_vert) != m_dist.end())
                continue;
            
            // save the distance
            m_dist[cur_vert] = cur_dist;
            // neighbors
            const Adjacency& nbrs = adj[cur_vert];
            for (std::size_t j = 0; j < nbrs.size(); ++j)
            {
                const Index neighbor = nbrs[j];

                const VecFloat cord_vert = points[cur_vert];
                const VecFloat cord_neighbor = points[neighbor];
                const float w = (cord_vert - cord_neighbor).length();
                const float new_dist = cur_dist + w;

                // only if shorter and in radius
                if (new_dist <= radius && m_dist.find(neighbor) == m_dist.end())
                {
                    m_queue.push(QE(new_dist, neighbor));
                }
            }
        }

        // resultlist out
        m_result.clear();
        for (auto it = m_dist.begin(); it != m_dist.end(); ++it)
        {
            const int idx = it->first;
            const float d = it->second;
            if (idx != start)
            {
                m_result.push_back(std::make_tuple(points[idx], d));
            }
        }

        return m_result;
    }

    const Distances<Index>& result(){
        return m_result;
    }
private:
    // min heap
    using QE = std::pair<float, Index>;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> m_queue;
    // distance hashmap
    std::unordered_map<Index, float> m_dist;

    // result in list form
    Distances<Index> m_result;
};