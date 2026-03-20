#pragma once

#include <array>
// #include <vector>
// #include <limits>
// #include <iostream>
// #include <unordered_set>

#include <utils/Vec.hpp>
#include <utils/Permutation.hpp>
// #include <io/wavefront.hpp>

namespace surface {
    /**
     * A list of point indices and nothing else. Used for the finalized Surface.
     */
    class SimpleTriangle {
    public:
        SimpleTriangle(std::array<std::size_t, 3> point_indices) : m_indices(point_indices){}
        SimpleTriangle(std::size_t first, std::size_t second, std::size_t third) : m_indices({first, second, third}){}

        std::array<std::size_t, 3> toRaw() const {
            return m_indices;
        }

        void index_translate(std::size_t translate){
            for(auto i = 0; i < 3; ++i){
                m_indices[i] += translate;
            }
        }

        SimpleTriangle map(const mutil::Permutation& perm) const {
            return SimpleTriangle(perm[m_indices[0]], perm[m_indices[1]], perm[m_indices[2]]);
        }

        bool operator==(const SimpleTriangle& rhs) const {
            // find first matching point of rhs
            std::size_t match = 3;
            for(std::size_t i = 0; i < 3; ++i){
                if(m_indices[0] == rhs.m_indices[i]){
                    match = i;
                    break;
                }
            }
            if(match > 2){
                return false;
            }
            return m_indices[1] == rhs.m_indices[(match + 1) % 3] && m_indices[2] == rhs.m_indices[(match + 2) % 3];
        }

        bool lexicographic_order(const SimpleTriangle& rhs) const {
            // find min of both
            auto lmin = min_index();
            auto rmin = rhs.min_index();

            for(std::size_t i = 0; i < 3; ++i){
                if(m_indices[(lmin + i) % 3] < rhs.m_indices[(rmin + i) % 3]){
                    return true;
                }else if(m_indices[(lmin + i) % 3] > rhs.m_indices[(rmin + i) % 3]){
                    return false;
                }
            }
            return false;
        }

        const std::size_t& operator[](int index) const {
            return m_indices[index];
        }

    private:
        std::size_t min_index() const {
            std::size_t index = 0;
            std::size_t min = m_indices[0];
            for(std::size_t i = 1; i < 3; ++i){
                if(m_indices[i] < min){
                    min = m_indices[i];
                    index = i;
                }
            }
            return index;
        }
        std::array<std::size_t, 3> m_indices;
    };
    
    inline std::ostream& operator<<(std::ostream& os, const SimpleTriangle& triangle)
    {
        os << "<" << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << ">";
        return os;
    }

    /**
     * Patch Id of -1 means we want to remove the triangle.
     */
    class Triangle {
    public:
        Triangle(std::array<std::size_t, 3> point_indices, int64_t patch) : m_points(point_indices), m_patch(patch){}
        std::array<std::size_t, 3> toRaw() const {
            return m_points;
        }

        void index_translate(std::size_t translate){
            for(auto i = 0; i < 3; ++i){
                m_points[i] += translate;
            }
        }

        Triangle map(const mutil::Permutation& perm) const {
            return Triangle({perm[m_points[0]], perm[m_points[1]], perm[m_points[2]]}, m_patch);
        }

        bool operator==(const Triangle& rhs) const {
            if(m_patch != rhs.m_patch){
                return false;
            }
            // find first matching point of rhs
            std::size_t match = 3;
            for(std::size_t i = 0; i < 3; ++i){
                if(m_points[0] == rhs.m_points[i]){
                    match = i;
                    break;
                }
            }
            if(match > 2){
                return false;
            }
            return m_points[1] == rhs.m_points[(match + 1) % 3] && m_points[2] == rhs.m_points[(match + 2) % 3];
        }

        bool lexicographic_order(const Triangle& rhs) const {
            if(m_patch < rhs.m_patch){
                return true;
            }else if(m_patch > rhs.m_patch){
                return false;
            }

            // find min of both
            auto lmin = min_index();
            auto rmin = rhs.min_index();

            for(std::size_t i = 0; i < 3; ++i){
                if(m_points[(lmin + i) % 3] < rhs.m_points[(rmin + i) % 3]){
                    return true;
                }else if(m_points[(lmin + i) % 3] > rhs.m_points[(rmin + i) % 3]){
                    return false;
                }
            }
            return false;
        }

        const std::size_t& operator[](int index) const {
            return m_points[index];
        }

        const int64_t patch_id() const {
            return m_patch;
        }

        uint64_t m_patch;
    private:
        std::size_t min_index() const {
            std::size_t index = 0;
            std::size_t min = m_points[0];
            for(std::size_t i = 1; i < 3; ++i){
                if(m_points[i] < min){
                    min = m_points[i];
                    index = i;
                }
            }
            return index;
        }
        std::array<std::size_t, 3> m_points;
    };

    inline std::ostream& operator<<(std::ostream& os, const Triangle& triangle)
    {
        os << "<" << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << ">";
        return os;
    }
}