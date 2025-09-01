#pragma once

#include <vector>
#include <optional>
#include <ostream>

class ChangeGraph
{
public:
    ChangeGraph();
    ChangeGraph(std::size_t size);

    void clear();
    // identity with futureSize
    void next();
    // past is clear, future is basically adding everything
    void reset();
    void identity(std::size_t size);
    void set(std::vector<std::size_t> elements);

    void change(std::size_t);
    void remove(std::size_t);
    void add(std::size_t);
    void add();

    /**
     * @brief Replaced the values of left with the values in right. (left.property = right.property)
     *
     * It is assumed that the data is then deleted in right (or counts as its own instance).
     *
     * @param from
     * @param to
     */
    void replace(std::size_t left, std::size_t right);

    std::size_t pastSize() const;
    std::size_t futureSize() const;
    std::size_t noneValue() const;

    bool wasDeleted(std::size_t) const;
    bool willBeAdded(std::size_t) const;
    bool unchanged(std::size_t) const;

    bool wasRenamed(std::size_t) const;
    bool willBeRenamed(std::size_t) const;

    std::optional<std::size_t> renamedFrom(std::size_t) const;
    std::optional<std::size_t> renamedBy(std::size_t) const;

    std::size_t renamedFromUnsafe(std::size_t) const;
    std::size_t renamedByUnsafe(std::size_t) const;

protected:
    std::vector<std::size_t> forward;
    std::vector<std::size_t> backward;
    std::size_t none_val; // std::numeric_limits<size_t>::max()
};

std::ostream&
operator<<(std::ostream& os, ChangeGraph const& g);

// /**
//  * @brief Class which returns replace events such that at some point all nodes are either unchanged, deleted or added.
//  *
//  * IMPORTANT: This class is unsafe, as it will never terminate for cyclic graphs. Your graph may be cyclic if you used add, delete and replace together.
//  * If one needs an replace iterator for cyclic graphs, one has to also check against the start (to detect cycles) and add two replace events with an element which does not exist.
//  * Thus in reality one needs to generate a temporary node.
//  *
//  */
// class ReplaceIteratorNoCycles
// {
// public:
//     using value_type = std::pair<std::size_t, std::size_t>;

//     ri::Option<value_type>
//     next()
//     {
//         if (start_node > graph.futureSize())
//         {
//             return ri::none;
//         }
//     }

// protected:
//     ChangeGraph& graph;
//     std::size_t start_node;
//     std::size_t current_node;
//     bool replacing;
// };
