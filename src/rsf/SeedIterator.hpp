#pragma once

#include <unordered_map>
#include <unordered_set>

#include <utils/IteratorEndSentinel.hpp>
#include <rsf/Seed.h>

namespace ridgesurface {

class SeedInSurfaceIterator {
public:
    using value_type = Seed;
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type;
    using reference         = value_type;

    SeedInSurfaceIterator(typename std::unordered_map<uint64_t, Seed>::const_iterator current, 
        typename std::unordered_map<uint64_t, Seed>::const_iterator end,
        const std::unordered_set<uint64_t>& seeds_in_surface)
        : m_current(current)
        , m_end(end)
        , m_seeds_in_surface(seeds_in_surface)
        {
            clip();
        }

    /**
     * We want to guarantee that current() always returns a valid edge as long as the iterator returns true (e.g. is not finished).
     * To achieve this, we must check for invalid states (m_target_current != m_target_end)
     */
    void clip()
    {
        while (m_current != m_end && !m_seeds_in_surface.contains(m_current->first))
        {
            ++m_current;
        }
    }

    void advance()
    {
        // only advance if we are not finished yet
        if (m_current == m_end)
        {
            return;
        }

        ++m_current;

        // move to next valid m_current
        clip();
    }

    bool finished() const
    {
        return m_current == m_end;
    }

    /** UB if called without checking if iterator is already finished */
    value_type current() const
    {
        return m_current->second;
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
    SeedInSurfaceIterator &operator++()
    {
        advance();
        return *this;
    }

    // Postfix increment
    SeedInSurfaceIterator operator++(int)
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    SeedInSurfaceIterator
    begin() const
    {
        return SeedInSurfaceIterator(m_current, m_end, m_seeds_in_surface);
    }

    iter::IteratorEndSentinel
    end() const
    {
        return {};
    }

private:
    typename std::unordered_map<uint64_t, Seed>::const_iterator m_current;
    typename std::unordered_map<uint64_t, Seed>::const_iterator m_end;
    const std::unordered_set<uint64_t>& m_seeds_in_surface;
};

class SeedAdditionIterator {
    public:
        using value_type = Seed;
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using pointer           = value_type;
        using reference         = value_type;
    
        SeedAdditionIterator(typename std::unordered_map<uint64_t, Seed>::const_iterator current, 
            typename std::unordered_map<uint64_t, Seed>::const_iterator end,
            const std::unordered_set<uint64_t>& seeds_in_surface)
            : m_current(current)
            , m_end(end)
            , m_seeds_in_surface(seeds_in_surface)
            {
                clip();
            }
    
        /**
         * We want to guarantee that current() always returns a valid edge as long as the iterator returns true (e.g. is not finished).
         * To achieve this, we must check for invalid states (m_target_current != m_target_end)
         */
        void clip()
        {
            while (m_current != m_end && m_seeds_in_surface.contains(m_current->first))
            {
                ++m_current;
            }
        }
    
        void advance()
        {
            // only advance if we are not finished yet
            if (m_current == m_end)
            {
                return;
            }
    
            ++m_current;
    
            // move to next valid m_current
            clip();
        }
    
        bool finished() const
        {
            return m_current == m_end;
        }
    
        /** UB if called without checking if iterator is already finished */
        value_type current() const
        {
            return m_current->second;
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
        SeedAdditionIterator &operator++()
        {
            advance();
            return *this;
        }
    
        // Postfix increment
        SeedAdditionIterator operator++(int)
        {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }
    
        SeedAdditionIterator
        begin() const
        {
            return SeedAdditionIterator(m_current, m_end, m_seeds_in_surface);
        }
    
        iter::IteratorEndSentinel
        end() const
        {
            return {};
        }
    
    private:
        typename std::unordered_map<uint64_t, Seed>::const_iterator m_current;
        typename std::unordered_map<uint64_t, Seed>::const_iterator m_end;
        const std::unordered_set<uint64_t>& m_seeds_in_surface;
    };

} // end namespace ridgesurface