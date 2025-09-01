#include "ChangeGraph.h"
#include <limits>

std::ostream&
operator<<(std::ostream& os, ChangeGraph const& g)
{
    // return os << g;
    os << "sizes: " << g.pastSize() << ", " << g.futureSize() << std::endl;
    for (std::size_t i = 0; i < g.pastSize(); ++i)
    {
        os << i << " -> " << g.renamedByUnsafe(i) << std::endl;
    }
    os << "///////////" << std::endl;
    for (std::size_t i = 0; i < g.futureSize(); ++i)
    {
        os << i << " <- " << g.renamedFromUnsafe(i) << std::endl;
    }
    return os;
}

ChangeGraph::ChangeGraph()
    : forward()
    , backward()
    , none_val(std::numeric_limits<size_t>::max())
{
}

ChangeGraph::ChangeGraph(std::size_t size)
    : forward()
    , backward()
    , none_val(std::numeric_limits<size_t>::max())
{
    forward.resize(size);
    clear();
}

void
ChangeGraph::clear()
{
    backward.clear();
    forward.clear();
    // for(std::size_t i = 0; i < forward.size(); ++i){
    //     forward[i] = none_val;
    // }
}

void
ChangeGraph::next()
{
    identity(backward.size());
}

void
ChangeGraph::reset()
{
    forward.clear();
    for(std::size_t i = 0; i < futureSize(); ++i){
        backward[i] = none_val;
    }
}

void
ChangeGraph::identity(std::size_t size)
{
    forward.resize(size);
    backward.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        forward[i] = i;
        backward[i] = i;
    }
}

void
ChangeGraph::set(std::vector<std::size_t> elements)
{
    auto result = std::vector<std::size_t>();
    result.reserve(elements.size());
    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        auto index = (elements[i] == none_val) ? none_val : backward[elements[i]];
        result.push_back(index);
    }
    backward = result;
    for (std::size_t i = 0; i < backward.size(); ++i)
    {
        auto id = backward[i];
        if (id != none_val)
        {
            forward[id] = i;
        }
    }
}

void
ChangeGraph::change(std::size_t id)
{
    if (willBeAdded(id))
    {
        return;
    }
    auto past_id = backward[id];
    forward[past_id] = none_val;
    backward[id] = none_val;
}

void
ChangeGraph::remove(std::size_t id)
{
    if (!willBeAdded(id))
    {
        auto past_id = backward[id];
        forward[past_id] = none_val;
    }
    backward.erase(backward.begin() + id);
    for (std::size_t i = id; i < futureSize(); ++i)
    {
        if (!willBeAdded(i))
        {
            forward[backward[i]] = i;
        }
    }
}

void
ChangeGraph::add(std::size_t pos)
{
    backward.insert(backward.begin() + pos, none_val);
    for (std::size_t i = pos + 1; i < futureSize(); ++i)
    {
        if (!willBeAdded(i))
        {
            forward[backward[i]] = i + 1;
        }
    }
}

void
ChangeGraph::add()
{
    add(futureSize());
}

void
ChangeGraph::replace(std::size_t left, std::size_t right)
{
    auto future_right = forward[right];
    auto future_left = forward[left];
    forward[right] = future_left;
    forward[left] = future_right;
    if(future_right != none_val){
        backward[future_right] = left;
    }
    if(future_left != none_val){
        backward[future_left] = right;
    }

    // auto past_right = backward[right];
    // auto past_left = backward[left];
    // // new connections
    // if (past_right != none_val)
    // {
    //     forward[past_right] = left;
    // }
    // backward[left] = past_right;
    // // remove old connections
    // if (past_left != none_val)
    // {
    //     forward[past_left] = none_val;
    // }
    // backward[right] = none_val;
}

std::size_t
ChangeGraph::pastSize() const
{
    return forward.size();
}

std::size_t
ChangeGraph::futureSize() const
{
    return backward.size();
}

std::size_t
ChangeGraph::noneValue() const
{
    return none_val;
}

bool
ChangeGraph::wasDeleted(std::size_t id) const
{
    return forward[id] == none_val;
}

bool
ChangeGraph::willBeAdded(std::size_t id) const
{
    return backward[id] == none_val;
}

bool
ChangeGraph::unchanged(std::size_t id) const
{
    return backward[id] == id;
}

bool
ChangeGraph::wasRenamed(std::size_t id) const
{
    return backward[id] != id && backward[id] != none_val;
}

bool
ChangeGraph::willBeRenamed(std::size_t id) const
{
    return forward[id] != id && forward[id] != none_val;
}

std::optional<std::size_t>
ChangeGraph::renamedFrom(std::size_t id) const
{
    if (backward[id] != id && backward[id] != none_val)
    {
        return backward[id];
    }
    return std::optional<std::size_t>();
}

std::optional<std::size_t>
ChangeGraph::renamedBy(std::size_t id) const
{
    if (forward[id] != id && forward[id] != none_val)
    {
        return forward[id];
    }
    return std::optional<std::size_t>();
}

std::size_t
ChangeGraph::renamedFromUnsafe(std::size_t id) const
{
    return backward[id];
}

std::size_t
ChangeGraph::renamedByUnsafe(std::size_t id) const
{
    return forward[id];
}
