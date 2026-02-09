#pragma once

#include <vector>
#include <cstdint>

class UnionFind {
public:
    static constexpr uint32_t leaf_bit = 0b1 << 31;
    UnionFind(uint32_t size){
        uf.resize(size, leaf_bit | 0b1);
    }

    /**
     * Returns the id of the whole union set.
     */
    uint32_t find(uint32_t x){
        // find root
        auto r = x;
        while((uf[r] & leaf_bit) != leaf_bit){
            r = uf[r];
        }

        // compress path to root
        while((uf[x] & leaf_bit) != leaf_bit){
            auto tmp = uf[x];
            uf[x] = r;
            x = tmp;
        }

        return r;
    }

    /**
     * Returns the id of the joined union set.
     */
    uint32_t join(uint32_t x, uint32_t y){
        x = find(x);
        y = find(y);
        if(x != y){
            // join x to y
            if(uf[x] < uf[y])
            {
                uf[y] += uf[x] ^ leaf_bit;
                uf[x] = y;
                return y;
            }
            // join y to x
            else {
                uf[x] += uf[y] ^ leaf_bit;
                uf[y] = x;
                return x;
            }
        }
        return x;
    }

    /**
     * Returns the number of elements inside the union set.
     */
    uint32_t size(uint32_t x) {
        x = find(x);
        return uf[x] ^ leaf_bit;
    }
    
private:
    std::vector<uint32_t> uf;
};