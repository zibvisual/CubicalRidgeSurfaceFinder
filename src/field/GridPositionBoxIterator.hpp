#pragma once

#include <utils/Vec.hpp>
#include <utils/Dims.hpp>
#include <utils/IteratorEndSentinel.hpp>

class GridPositionBoxIterator {
public:
    using value_type = VecSize;
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type;
    using reference         = value_type;

    GridPositionBoxIterator(VecSize min, VecSize max) : m_pos(min), min(min), max(max) {}
    GridPositionBoxIterator(VecSize pos, VecSize min, VecSize max) : m_pos(pos), min(min), max(max) {}

    void advance()
    {
        m_pos[0] += 1;
        if(m_pos[0] >= max[0]){
            m_pos[0] = min[0];
            m_pos[1] += 1;
            if(m_pos[1] >= max[1]){
                m_pos[1] = min[1];
                m_pos[2] += 1;
            }
        }
    }

    bool finished() const
    {
        return m_pos[2] >= max[2];
    }

    /** UB if called without checking if iterator is already finished */
    value_type current() const
    {
        return m_pos;
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

    reference operator*() const
    {
        return current();
    }
    pointer operator->()
    {
        return current();
    }

    // Prefix increment
    GridPositionBoxIterator &operator++()
    {
        advance();
        return *this;
    }

    // Postfix increment
    GridPositionBoxIterator operator++(int)
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    GridPositionBoxIterator
    begin() const
    {
        return GridPositionBoxIterator(m_pos, min, max);
    }

    iter::IteratorEndSentinel
    end() const
    {
        return {};
    }

private:
    VecSize m_pos;
    VecSize min;
    VecSize max;
};