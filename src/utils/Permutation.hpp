#pragma once

#include <vector>
#include <optional>
#include <iostream>
#include <nanoflann.hpp>
#include <utils/PointCloudAdaptor.hpp>

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
            std::sort(index.begin(), index.end(),
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
            std::sort(index.begin(), index.end(), 
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

        template <class T>
        static Permutation mapPermutationNearestNeighbor(std::vector<T> source, std::vector<T> target, T::value_type threshold = 0.0001){
            if(source.size() != target.size()){
                std::cout << "unequal number of elements" << std::endl;
                return Permutation();
            }
            using pointcloud_t = PointCloudAdaptor<T>;
            auto src_cloud = pointcloud_t(source);
            auto tgt_cloud = pointcloud_t(target);

            // find pairing source<->target <distance, dataset adaptor, dim>
            using kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<typename T::value_type, pointcloud_t>,pointcloud_t,3>;
            kd_tree_t src_tree(3 /* dim (gets ignored because of template dim) */, src_cloud, {50 /* max leaf */});
            kd_tree_t tgt_tree(3 /* dim (gets ignored because of template dim) */, tgt_cloud, {50 /* max leaf */});


            auto perm = std::vector<std::size_t>();
            perm.reserve(source.size());

            uint32_t src_index;
            uint32_t tgt_index;
            typename T::value_type sqr_distance;
            for(std::size_t i = 0; i < source.size(); ++i){
                tgt_tree.knnSearch(&source[i][0], 1, &tgt_index, &sqr_distance);
                if(sqr_distance > threshold){
                    // matching distance threshold to big
                    return Permutation();
                }
                src_tree.knnSearch(&target[tgt_index][0], 1, &src_index, &sqr_distance);
                if(src_index != i){
                    // there is no definite closest pair matching
                    return Permutation();
                }
                perm.push_back(tgt_index);
            }
            return Permutation(perm);
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
