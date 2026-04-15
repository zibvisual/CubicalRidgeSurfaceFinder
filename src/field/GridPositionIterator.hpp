#pragma once

#include <utils/Vec.hpp>
#include <utils/Dims.hpp>
#include <utils/IteratorEndSentinel.hpp>

class GridPositionIterator {
public:
    using value_type = VecSize;
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type;
    using reference         = value_type;

    GridPositionIterator(Dims dims) : m_pos(VecSize(0)), m_dims(dims) {}
    GridPositionIterator(VecSize pos, Dims dims) : m_pos(pos), m_dims(dims) {}

    void advance()
    {
        m_pos[0] += 1;
        if(m_pos[0] >= m_dims[0]){
            m_pos[0] = 0;
            m_pos[1] += 1;
            if(m_pos[1] >= m_dims[1]){
                m_pos[1] = 0;
                m_pos[2] += 1;
            }
        }
    }

    bool finished() const
    {
        return m_pos[2] >= m_dims[2];
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
    GridPositionIterator &operator++()
    {
        advance();
        return *this;
    }

    // Postfix increment
    GridPositionIterator operator++(int)
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    GridPositionIterator
    begin() const
    {
        return GridPositionIterator(m_pos, m_dims);
    }

    iter::IteratorEndSentinel
    end() const
    {
        return {};
    }

private:
    VecSize m_pos;
    Dims m_dims;
};