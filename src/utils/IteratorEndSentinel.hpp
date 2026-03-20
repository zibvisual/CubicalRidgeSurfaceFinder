#pragma once


namespace iter {

/**
 * @brief The IteratorEndSentinel can be compared to every iterator which is convertable to a boolean.
 * 
 * The other iterator should be convertable to true as long as it did not reach the end, false otherwise.
 *
 */
struct IteratorEndSentinel
{
    template <typename Iter>
    friend bool
    operator==(const Iter& a, IteratorEndSentinel)
    {
        return !a;
    }

    template <typename Iter>
    friend bool
    operator!=(const Iter& a, IteratorEndSentinel)
    {
        return a;
    }

    template <typename Iter>
    friend bool
    operator==(IteratorEndSentinel, const Iter& a)
    {
        return !a;
    }

    template <typename Iter>
    friend bool
    operator!=(IteratorEndSentinel, const Iter& a)
    {
        return a;
    }

    template <typename Iter>
    friend bool
    operator==(IteratorEndSentinel, IteratorEndSentinel)
    {
        return true;
    }

    template <typename Iter>
    friend bool
    operator!=(IteratorEndSentinel, IteratorEndSentinel)
    {
        return false;
    }
};

}
