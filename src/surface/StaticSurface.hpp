#pragma once

#include <array>
#include <vector>
#include <limits>
#include <iostream>
#include <unordered_set>

#include <utils/Vec.hpp>
#include <utils/Permutation.hpp>
#include <io/wavefront.hpp>
#include <surface/Triangle.hpp>

namespace surface {
    /**
     * A simple triangulated surface mesh. Used for the finalized output.
     */
    class StaticSurface {
    public:
        StaticSurface() : m_points(), m_triangles(){}

        static StaticSurface load(std::filesystem::path input){
            auto data = read_wavefront(input);
            return StaticSurface::fromRaw(data);
        }

        void save(std::filesystem::path output) const {
            write_wavefront(output, toRaw());
        }

        static StaticSurface fromRaw(wavefront_data_t data){
            auto surface = StaticSurface(data.points, data.triangles);
            return surface;
        }
        
        wavefront_data_t toRaw() const {
            // convert triangles and use empty patch (reinterpret_cast should work here)
            std::vector<std::array<std::size_t, 3>> triangles;
            triangles.reserve(m_triangles.size());
            for(std::size_t i = 0; i < m_triangles.size(); ++i){
                triangles.push_back(m_triangles[i].toRaw());
            }
            return wavefront_data_t{m_points, triangles};
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
            return mutil::Permutation::mapPermutationNearestNeighbor(m_points, rhs.m_points, epsilon);
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
            m_triangles.push_back(SimpleTriangle(point1, point2, point3));
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