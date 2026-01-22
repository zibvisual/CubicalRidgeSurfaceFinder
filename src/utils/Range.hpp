#pragma once

#include <optional>

// clang-format off
#define ADD_RANGE_TO_RUST_ITERATOR(class)                     \
    ::iterator::ForRange<class> begin(){return ::iterator::ForRange<class>(this);}           \
    ::iterator::ForRange<class> end(){return ::iterator::ForRange<class>();}
// clang-format on

namespace iterator {
    /**
     * @brief Implementation of a range to be able to use rust iterators inside for loops
     *
     * The use of this adoptor might create a warning that a value might be uninitialized. This is probably created by the usage of boost::optional.
     * See https://stackoverflow.com/questions/21755206/how-to-get-around-gcc-void-b-4-may-be-used-uninitialized-in-this-funct
     *
     * Could be a false positive...
     */
    template <class Iter>
    class ForRange
    {
        using iterator_category = std::input_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = typename Iter::value_type;
        using pointer = value_type*;
        using reference = value_type&;

    public:
        ForRange()
            : current(std::optional<value_type>())
            , iter(nullptr)
        {}

        ForRange(Iter* iter)
            : current(iter->next())
            , iter(iter)
        {}

        ForRange(const ForRange& copy)
            : current(copy.current)
            , iter(copy.iter)
        {}

        value_type operator*() const
        {
            return current.value();
        }

        value_type operator->() const
        {
            return current.value();
        }

        ForRange& operator++()
        {
            current = iter->next();
            return *this;
        }

        ForRange& operator++(int)
        {
            ForRange tmp = *this;
            ++(*this);
            return tmp;
        }

        // This function should only be used for the for loop!
        friend bool operator==(const ForRange& a, const ForRange& b)
        {
            // Only iterators which are empty are the same!
            return !a.current.has_value() && !b.current.has_value();
        };

        friend bool operator!=(const ForRange& a, const ForRange& b)
        {
            // Only iterators which are empty are the same!
            return a.current.has_value() || b.current.has_value();
        };

    protected:
        std::optional<value_type> current;
        Iter* iter;
    };

}
