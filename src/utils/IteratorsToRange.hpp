#pragma once

template <typename I>
struct iter_pair : std::pair<I, I>
{ 
    using std::pair<I, I>::pair;

    I begin() { return this->first; }
    I end() { return this->second; }
};
