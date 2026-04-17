#pragma once

#include <utils/Vec.hpp>
#include <utils/Dims.hpp>
#include <utils/IteratorEndSentinel.hpp>

class GridPositionBorderIterator {
public:
    using value_type = std::pair<VecSize, float>;
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type;
    using reference         = value_type;

    GridPositionBorderIterator(Dims dims, std::size_t thickness = 1) : m_pos(VecSize(0)), m_dims(dims), m_t(thickness) {}
    GridPositionBorderIterator(VecSize pos, Dims dims, std::size_t thickness = 1) : m_pos(pos), m_dims(dims), m_t(thickness) {}

    void advance()
    {
        if(m_mode == 0){
            // first x and y free, z < t or n-1-z < t
            m_pos[0] += 1;
            if(m_pos[0] < m_dims[0]) return;
            m_pos[0] = 0;
            m_pos[1] += 1;
            if(m_pos[1] < m_dims[1]) return;
            m_pos[1] = 0;
            m_pos[2] +=1;
            if(m_pos[2] < m_t) return;
            if(m_pos[2] < m_dims[2] - m_t){
                m_pos[2] = m_dims[2] - m_t;
                return;
            }
            if(m_pos[2] < m_dims[2]) return;
            // next mode (starts with 0,0,t)
            m_pos[2] = m_t;
            m_mode = 1;
        }
        else if(m_mode == 1)
        {
            // next step: x and z free, y < t or n-1-y < t AND z NOT
            m_pos[0] += 1;
            if(m_pos[0] < m_dims[0]) return;
            m_pos[0] = 0;
            m_pos[1] += 1;
            if(m_pos[1] < m_t) return;
            if(m_pos[1] < m_dims[1] - m_t){
                m_pos[1] = m_dims[1] - m_t;
                return;
            }
            if(m_pos[1] < m_dims[1]) return;
            m_pos[1] = 0;
            m_pos[2] += 1;
            if(m_pos[2] < m_dims[2] - m_t) return;
            // next mode (starts with 0,t,t)
            m_pos[2] = m_t;
            m_pos[1] = m_t;
            m_mode = 2;
        }
        else if(m_mode == 2)
        {
            // next step: x < t or n-1-x < t AND y,z NOT
            m_pos[0] += 1;
            if(m_pos[0] < m_t) return;
            if(m_pos[0] < m_dims[0] - m_t){
                m_pos[0] = m_dims[0] - m_t;
                return;
            }
            if(m_pos[0] < m_dims[0]) return;
            m_pos[0] = 0;
            m_pos[1] += 1;
            if(m_pos[1] < m_dims[1] - m_t) return;
            m_pos[1] = m_t;
            m_pos[2] += 1;
            if(m_pos[2] < m_dims[2] - m_t) return;
            // finished
            m_mode = 3;
        }
    }

    bool finished() const
    {
        return m_mode == 3;
    }

    /** UB if called without checking if iterator is already finished */
    value_type current() const
    {
        // calculate current thickness (returns [0.0,1.0))
        auto left = static_cast<int>(m_pos.min());
        auto right = (static_cast<VecInt>(m_dims) - static_cast<VecInt>(m_pos)).min() - 1;
        auto dist = std::min(left, right);
        auto thickness = static_cast<float>(dist) / m_t;
        return {m_pos, thickness};
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
    GridPositionBorderIterator &operator++()
    {
        advance();
        return *this;
    }

    // Postfix increment
    GridPositionBorderIterator operator++(int)
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    GridPositionBorderIterator
    begin() const
    {
        return GridPositionBorderIterator(m_pos, m_dims, m_t);
    }

    iter::IteratorEndSentinel
    end() const
    {
        return {};
    }

    // other functions:
    std::size_t thickness() const
    {
        return m_t;
    }

private:
    VecSize m_pos;
    Dims m_dims;
    // thickness
    std::size_t m_t;
    // mode:
    int m_mode = 0;
};