#pragma once

#include <utils/Dims.hpp>
#include <utils/IteratorEndSentinel.hpp>

namespace field {

class CLocationPairScanner {
public:
    using value_type = std::pair<std::size_t, std::size_t>;
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type;
    using reference         = value_type;

    CLocationPairScanner(Dims dims) : dims(dims) {}
    // GridPositionIterator(VecSize pos, Dims dims) : m_pos(pos), m_dims(dims) {}

    void advance()
    {
        // advance normally
        ++x;
        if(x < dims[0] - dir_x){
            return;
        }
        x = 0;
        ++y;
        if(y < dims[1] - dir_y){
            return;
        }
        y = 0;
        ++z;
        if(z < dims[2] - dir_z){
            return;
        }
        z = 0;

        // go to next direction
        if(dir_x){
            dir_x = false;
            dir_y = true;
            offset = dims[0];
        }else if(dir_y){
            dir_y = false;
            dir_z = true;
            offset = dims[0] * dims[1];
        }else if(dir_z){
            offset = 0;
        }
    }

    bool finished() const
    {
        return offset == 0;
    }

    /** UB if called without checking if iterator is already finished */
    value_type current() const
    {
        std::size_t index = z * dims.x() * dims.y() + y * dims.x() + x;;
        return {index, index + offset};
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
    CLocationPairScanner &operator++()
    {
        advance();
        return *this;
    }

    // Postfix increment
    CLocationPairScanner operator++(int)
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    CLocationPairScanner
    begin() const
    {
        return CLocationPairScanner(dims);
    }

    iter::IteratorEndSentinel
    end() const
    {
        return {};
    }

private:
    Dims dims;
    std::size_t x = 0;
    std::size_t y = 0;
    std::size_t z = 0;
    bool dir_x = true;
    bool dir_y = false;
    bool dir_z = false;
    std::size_t offset = 1;
};

class CLocationPairScannerX {
    public:
        using value_type = std::pair<std::size_t, std::size_t>;
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using pointer           = value_type;
        using reference         = value_type;
    
        CLocationPairScannerX(Dims dims) : dims(dims) {}
        // GridPositionIterator(VecSize pos, Dims dims) : m_pos(pos), m_dims(dims) {}
    
        void advance()
        {
            // advance normally
            ++x;
            if(x < dims[0] - 1){
                return;
            }
            x = 0;
            ++y;
            if(y < dims[1]){
                return;
            }
            y = 0;
            ++z;
            if(z < dims[2]){
                return;
            }
        }
    
        bool finished() const
        {
            return z >= dims[2];
        }
    
        /** UB if called without checking if iterator is already finished */
        value_type current() const
        {
            std::size_t index = z * dims.x() * dims.y() + y * dims.x() + x;;
            return {index, index + 1};
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
        CLocationPairScannerX &operator++()
        {
            advance();
            return *this;
        }
    
        // Postfix increment
        CLocationPairScannerX operator++(int)
        {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }
    
        CLocationPairScannerX
        begin() const
        {
            return CLocationPairScannerX(dims);
        }
    
        iter::IteratorEndSentinel
        end() const
        {
            return {};
        }
    
    private:
        Dims dims;
        std::size_t x = 0;
        std::size_t y = 0;
        std::size_t z = 0;
    };

    class CLocationPairScannerY {
    public:
        using value_type = std::pair<std::size_t, std::size_t>;
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using pointer           = value_type;
        using reference         = value_type;
    
        CLocationPairScannerY(Dims dims) : dims(dims) {}
        // GridPositionIterator(VecSize pos, Dims dims) : m_pos(pos), m_dims(dims) {}
    
        void advance()
        {
            // advance normally
            ++x;
            if(x < dims[0]){
                return;
            }
            x = 0;
            ++y;
            if(y < dims[1] - 1){
                return;
            }
            y = 0;
            ++z;
            if(z < dims[2]){
                return;
            }
        }
    
        bool finished() const
        {
            return z >= dims[2];
        }
    
        /** UB if called without checking if iterator is already finished */
        value_type current() const
        {
            std::size_t index = z * dims.x() * dims.y() + y * dims.x() + x;;
            return {index, index + dims[0]};
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
        CLocationPairScannerY &operator++()
        {
            advance();
            return *this;
        }
    
        // Postfix increment
        CLocationPairScannerY operator++(int)
        {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }
    
        CLocationPairScannerY
        begin() const
        {
            return CLocationPairScannerY(dims);
        }
    
        iter::IteratorEndSentinel
        end() const
        {
            return {};
        }
    
    private:
        Dims dims;
        std::size_t x = 0;
        std::size_t y = 0;
        std::size_t z = 0;
    };

class CLocationPairScannerZ {
public:
    using value_type = std::pair<std::size_t, std::size_t>;
    using iterator_category = std::forward_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type;
    using reference         = value_type;

    CLocationPairScannerZ(Dims dims) : dims(dims) {}
    // GridPositionIterator(VecSize pos, Dims dims) : m_pos(pos), m_dims(dims) {}

    void advance()
    {
        // advance normally
        ++x;
        if(x < dims[0]){
            return;
        }
        x = 0;
        ++y;
        if(y < dims[1]){
            return;
        }
        y = 0;
        ++z;
        if(z < dims[2] - 1){
            return;
        }
    }

    bool finished() const
    {
        return z >= dims[2] - 1;
    }

    /** UB if called without checking if iterator is already finished */
    value_type current() const
    {
        std::size_t index = z * dims.x() * dims.y() + y * dims.x() + x;;
        return {index, index + dims[0] * dims[1]};
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
    CLocationPairScannerZ &operator++()
    {
        advance();
        return *this;
    }

    // Postfix increment
    CLocationPairScannerZ operator++(int)
    {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }

    CLocationPairScannerZ
    begin() const
    {
        return CLocationPairScannerZ(dims);
    }

    iter::IteratorEndSentinel
    end() const
    {
        return {};
    }
    
private:
    Dims dims;
    std::size_t x = 0;
    std::size_t y = 0;
    std::size_t z = 0;
};

} // namespace field
