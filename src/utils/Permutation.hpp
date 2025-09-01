#pragma once

#include <vector>
#include <optional>
#include <iostream>

namespace mutil {
    class Permutation {
    public:
        /**
         * Zero Size Permutation. Used to represent fallable calculations.
         */
        Permutation() : m_perm(){}
        Permutation(std::vector<std::size_t> perm) : m_perm(perm) {
            //TODO: check if it is really a permutation (we have values 0..n)
        }

        template <class T>
        static Permutation sortPermutation(std::vector<T> data){
            std::vector<std::size_t> index(data.size(), 0);
            for (std::size_t i = 0 ; i < index.size() ; ++i) {
                index[i] = i;
            }
            sort(index.begin(), index.end(),
                [&](const int& a, const int& b) {
                    return (data[a] < data[b]);
                }
            );
            return Permutation{index};
        }

        template <class T, class Compare>
        static Permutation sortPermutation(std::vector<T> data, Compare comp){
            std::vector<std::size_t> index(data.size(), 0);
            for (std::size_t i = 0 ; i < index.size() ; ++i) {
                index[i] = i;
            }
            sort(index.begin(), index.end(), 
                [&](const int& a, const int& b) {
                    return comp(data[a], data[b]);
                }
            );
            return Permutation{index};
        }

        template <class T>
        static Permutation mapPermutation(std::vector<T> source, std::vector<T> target){
            if(source.size() != target.size()){
                return Permutation();
            }
            auto src_perm = sortPermutation(source);
            auto trg_perm = sortPermutation(target);
            //check if both source and target have the same values!
            for(std::size_t i = 0; i < source.size(); ++i){
                const auto a = source[src_perm[i]];
                const auto b = target[trg_perm[i]];
                if(a != b){
                    return Permutation();
                }
            }
            auto map = src_perm.apply(trg_perm.inverse());
            return map;
        }

        template <class T, class Compare>
        static Permutation mapPermutation(std::vector<T> source, std::vector<T> target, Compare comp){
            if(source.size() != target.size()){
                return Permutation();
            }
            auto src_perm = sortPermutation(source, comp);
            auto trg_perm = sortPermutation(target, comp);
            //check if both source and target have the same values!
            for(std::size_t i = 0; i < source.size(); ++i){
                const auto a = source[src_perm[i]];
                const auto b = target[trg_perm[i]];
                // if a != b then return empty permutation
                if(comp(a, b) || comp(b, a)){
                    std::cout << "---- a:" << a << " and b: " << b << std::endl;
                    return Permutation();
                }
            }
            auto map = trg_perm.apply(src_perm.inverse());
            return map;
        }

        std::size_t size() const {
            return m_perm.size();
        }

        bool empty() const {
            return m_perm.size() == 0;
        }

        Permutation inverse() const {
            std::vector<std::size_t> res(size(), 0);
            for(std::size_t i = 0; i < size(); ++i){
                res[m_perm[i]] = i;
            }
            return Permutation(res);
        }

        /**
         * Apply rhs permutation onto current permutation
         */
        Permutation apply(const Permutation& rhs){
            // apply rhs permutation
            std::vector<std::size_t> res = m_perm;
            for(std::size_t i = 0; i < rhs.size(); ++i){
                res[i] = m_perm[rhs.m_perm[i]];
            }
            return Permutation(res);
        }

        const std::size_t& operator[](std::size_t index) const {
            return m_perm[index];
        }

        const std::size_t& operator[](int index) const {
            return m_perm[index];
        }

        // // same as / (we do not overload as we allocate a new vec)
        // Permutation apply_inverse(const Permutation& rhs){
        //     // TODO apply the inverse of the rhs (without allocating twice)?
        //     std::vector<std::size_t> res = m_perm;
        //     for(std::size_t i = 0; i < rhs.size(); ++i){
        //         res[i] = m_perm[rhs.m_perm[i]];
        //     }
        //     return Permutation(res);
        // }

    private:
        std::vector<std::size_t> m_perm;
    };

    inline std::ostream& operator<<(std::ostream& os, const mutil::Permutation& perm)
    {
        os << "[";
        for(std::size_t i = 0; i < perm.size(); ++i){
            os << perm[i];
            if(i + 1 < perm.size()){
                os << " ";
            }
        }
        os << "]";
        return os;
    }
}
