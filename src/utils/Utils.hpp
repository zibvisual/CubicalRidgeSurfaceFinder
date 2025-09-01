#pragma once

#include <cmath>

namespace mutil {

    // CAPPING
    inline double roundUpWithBiggestPowerOf10(double input)
    {
        double log10val = std::floor(log10(input));
        double modulo = std::pow(10.0f, log10val);
        return input - fmod(input, modulo) + modulo;
    }

    template <typename T>
    T capValues(T* dataPtr, std::size_t size, T max)
    {
        T cap = static_cast<T>(roundUpWithBiggestPowerOf10(static_cast<double>(max)));
        for (std::size_t i = 0; i < size; ++i)
        {
            dataPtr[i] = dataPtr[i] < cap ? dataPtr[i] : cap;
        }
        return cap;
    }

}