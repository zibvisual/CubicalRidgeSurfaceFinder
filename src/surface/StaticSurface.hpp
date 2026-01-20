#pragma once

#include <array>
#include <vector>
#include <limits>
#include <iostream>
#include <unordered_set>

#include <utils/Vec.hpp>
#include <utils/Permutation.hpp>
#include <io/wavefront.hpp>

namespace surface {
    /**
     * A list of point indices and nothing else. Used for the finalized Surface.
     */
    class SimpleTriangle {
    public:
        SimpleTriangle(std::array<std::size_t, 3> point_indices) : m_indices(point_indices){}
        std::array<std::size_t, 3> toRaw() const {
            return m_indices;
        }

        SimpleTriangle map(const mutil::Permutation& perm) const {
            return SimpleTriangle({perm[m_indices[0]], perm[m_indices[1]], perm[m_indices[2]]});
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
     * A simple triangulated surface mesh. Used for the finalized output.
     */
    class StaticSurface {
    public:
        StaticSurface() : m_points(), m_triangles(){}

        static StaticSurface load(std::string input_path){
            auto data = read_wavefront(input_path);
            return StaticSurface::fromRaw(data);
        }

        void save(std::string output_path) const {
            write_wavefront(output_path, toRaw());
        }

        static StaticSurface fromRaw(wavefront_data_t data){
            auto surface = StaticSurface(data.points, data.triangles);
            return surface;
        }
        
        wavefront_data_t toRaw() const {
            // convert triangles and use empty patch (reinterpret_cast should work here)
            std::vector<std::array<std::size_t, 3>> triangles;
            std::vector<std::size_t> patches;
            triangles.reserve(m_triangles.size());
            for(std::size_t i = 0; i < m_triangles.size(); ++i){
                triangles.push_back(m_triangles[i].toRaw());
            }
            return wavefront_data_t{m_points, triangles, patches};
        }

        StaticSurface clone() const {
            return StaticSurface(
                m_points,
                m_triangles
            );
        }

        void removeUnusedPoints(){
            // TODO: go through all triangles and collect indices of the points. 
        }

        bool validSurface() const;
        bool validManifold() const;
        // void addPatch() --> ?
        // void replacePatch() --> ?

        mutil::Permutation mapEqualPoints(const StaticSurface& rhs, float epsilon = 0.0001) const {
            // first check if we have the same amount of points
            if(m_points.size() != rhs.m_points.size()){
                return mutil::Permutation();
            }

            auto perm = mutil::Permutation::mapPermutation(m_points, rhs.m_points, 
                [&](const VecFloat& a, const VecFloat& b) {
                    return a.lexicographic_order_epsilon(b, epsilon);
                }
            );
            return perm;
        }


        mutil::Permutation mapEqualTriangles(const StaticSurface& rhs) const {
            auto perm = mutil::Permutation::mapPermutation(m_triangles, rhs.m_triangles, 
                [&](const SimpleTriangle& a, const SimpleTriangle& b) {
                    return a.lexicographic_order(b);
                }
            );
            return perm;
        }


        bool equals(const StaticSurface& rhs, float epsilon = 0.0001) const {
            // TODO: delete all unused points and all obsolete triangles

            auto pointMap = mapEqualPoints(rhs, epsilon);
            if(pointMap.empty()){
                return false;
            }

            auto triangleMap = permutePoints(pointMap).mapEqualTriangles(rhs);
            if(triangleMap.empty()){
                std::cout << "triangleMap null" << std::endl;
                return false;
            }
            return true;
        }

        const std::vector<SimpleTriangle>& triangles() const {
            return m_triangles;
        }

        const std::vector<VecFloat>& points() const {
            return m_points;
        }

        std::size_t number_of_trianlges() const {
            return m_triangles.size();
        }

        std::size_t number_of_points() const {
            return m_points.size();
        }

        // Methods to modify Surface. Usually not used directly.
        void clear() {
            m_points.clear();
            m_triangles.clear();
        }

        void addPoint(VecFloat point){
            m_points.push_back(point);
        }

        void addTriangle(std::size_t point1, std::size_t point2, std::size_t point3){
            m_triangles.push_back(SimpleTriangle({point1, point2, point3}));
        }

        StaticSurface permutePoints(const mutil::Permutation& perm) const{
            if(perm.size() != m_points.size()){
                return StaticSurface();
            }
            // permute all points
            auto points = m_points;
            for(std::size_t i = 0; i < points.size(); ++i){
                points[i] = m_points[perm[i]];
            }
            // permute indices of triangles
            auto triangles = m_triangles;
            for(std::size_t i = 0; i < triangles.size(); ++i){
                triangles[i] = m_triangles[i].map(perm);
            }
            return StaticSurface(points, triangles);
        }

    private:
        StaticSurface(std::vector<VecFloat> points, std::vector<std::array<std::size_t, 3>> triangles)
        : m_points(points)
        , m_triangles()
        {
            // convert triangles
            m_triangles.reserve(triangles.size());
            for(std::size_t i = 0; i < triangles.size(); ++i){
                m_triangles.push_back(SimpleTriangle(triangles[i]));
            }
            // TODO: check if surface is valid
        }
        StaticSurface(std::vector<VecFloat> points, std::vector<SimpleTriangle> triangles)
        : m_points(points)
        , m_triangles(triangles)
        {
        }
    private:
        std::vector<VecFloat> m_points;
        // contain the indices of the points
        std::vector<SimpleTriangle> m_triangles;
    };
}