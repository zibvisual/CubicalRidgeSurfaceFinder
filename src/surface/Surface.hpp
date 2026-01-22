#pragma once

#include <array>
#include <vector>
#include <limits>
#include <iostream>
#include <unordered_set>

#include <utils/Vec.hpp>
#include <utils/Permutation.hpp>
#include <io/wavefront.hpp>

// TODO: write wavefront files with this surface! This is the most important, so we can debug!
// TODO: loading would also be nice for automatic testing!

// wavefront:
// # point:
// v float float float
// # normal:
// vn float float float
// # face
// f index index index
// # patches
// g groupname/groupnumber


namespace surface {
    /**
     * Patch Id of zero means we want to remove the triangle.
     */
    class Triangle {
    public:
        Triangle(std::array<std::size_t, 3> point_indices, int64_t patch) : m_points(point_indices), m_patch(patch){}
        std::array<std::size_t, 3> toRaw() const {
            return m_points;
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

    struct PatchInformation {
        uint64_t m_patch;
        bool orientation;
    };

    /**
     * A triangulated surface mesh with patches (triangle groups). Triangles have a patch id.
     * A list of PatchInformation is also given. They contain the patch_id and if the orientation of the patch is flipped.
     * If a patch id is not contained in the PatchInformation list, one may assume the patch to be deleted and not used anymore.
     */
    class Surface {
    public:
        Surface() : m_points(), m_triangles(){}

        static Surface load(std::string input_path){
            auto data = read_wavefront(input_path);
            return Surface::fromRaw(data);
        }

        void save(std::string output_path) const {
            write_wavefront(output_path, toRaw());
        }

        static Surface fromRaw(wavefront_data_t data){
            auto surface = Surface(data.points, data.triangles, data.patches);
            return surface;
        }
        
        wavefront_data_t toRaw() const {
            // TODO: sort triangles depending on their patch, remove all obsolete triangles and unused points
            // convert triangles
            std::vector<std::array<std::size_t, 3>> triangles;
            std::vector<std::size_t> patches;
            triangles.reserve(m_triangles.size());
            std::size_t lastPatch = 0;
            for(std::size_t i = 0; i < m_triangles.size(); ++i){
                triangles.push_back(m_triangles[i].toRaw());
                if(lastPatch != m_triangles[i].patch_id()){
                    lastPatch = m_triangles[i].patch_id();
                    patches.push_back(i);
                }
            }
            return wavefront_data_t{m_points, triangles, patches};
        }

        Surface clone() const {
            return Surface(
                m_points,
                m_triangles
            );
        }

        void removeUnusedPoints(){
            // TODO: go through all triangles and collect indices of the points. 
        }

        /**
         * If mulitple patches should be removed, prefer removePatches()
         */
        void removePatch(uint64_t id){

            m_triangles.erase(std::remove_if(begin(m_triangles), end(m_triangles), [id](const auto& triangle) { return triangle.patch_id() == id;}));
            removeUnusedPoints();
        }

        void removePatches(const std::unordered_set<uint64_t>& ids)
        {
            m_triangles.erase(std::remove_if(begin(m_triangles), end(m_triangles), [&ids](const auto& triangle) { return ids.contains(triangle.patch_id());}));
            removeUnusedPoints();
        }

        bool validSurface() const;
        bool validManifold() const;
        // void addPatch() --> ?
        // void replacePatch() --> ?

        // // usually you want to use equals with an epsilon instead
        // bool operator==(const Surface& rhs) const {
        //     // first check if we have the same amount of points and triangles and patches
        //     if(m_points.size() != rhs.m_points.size() || m_triangles.size() != rhs.m_triangles.size() || m_patches != rhs.m_patches){
        //         return false;
        //     }

        //     // given map permutation, we can check if we have the same triangles!
        //     auto perm = mutil::Permutation::mapPermutation(m_points, rhs.m_points, 
        //         [&](const VecFloat& a, const VecFloat& b) {
        //             return a.lexicographic_order_exact(b);
        //         }
        //     );

        //     if(perm.empty()){
        //         return false;
        //     }

        //     // we have the same points, so now we can check triangles
        //     for(std::size_t i = 0; i < m_triangles.size(); ++i){
        //         // if(m_triangles[i] != rhs.m_triangles[]) 
        //     }
        // }

        mutil::Permutation mapEqualPoints(const Surface& rhs, float epsilon = 0.0001) const {
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


        mutil::Permutation mapEqualTriangles(const Surface& rhs) const {
            auto perm = mutil::Permutation::mapPermutation(m_triangles, rhs.m_triangles, 
                [&](const Triangle& a, const Triangle& b) {
                    return a.lexicographic_order(b);
                }
            );
            return perm;
        }


        bool equals(const Surface& rhs, float epsilon = 0.0001) const {
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

            //TODO: check if we have the same patches (with regards to triangleMap)
            // we assume that there are not duplicate triangles with different patches
            // // first check if we have the same amount of patches
            // if(m_patches.size() != rhs.m_patches.size()){
            //     return false;
            // }
            return true;
        }

        const std::vector<Triangle>& triangles() const {
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

        /** we count */
        std::size_t number_of_patches() const {
            auto set = std::unordered_set<int64_t>();
            for(const auto triangle : m_triangles){
                set.insert(triangle.patch_id());
            }
            return set.size();
        }

        // Methods to modify Surface. Usually not used directly.
        void clear() {
            m_points.clear();
            m_triangles.clear();
        }

        void addPoint(VecFloat point){
            m_points.push_back(point);
        }

        void addTriangle(std::size_t point1, std::size_t point2, std::size_t point3, std::size_t patch){
            m_triangles.push_back(Triangle({point1, point2, point3}, patch));
        }

        Surface permutePoints(const mutil::Permutation& perm) const{
            if(perm.size() != m_points.size()){
                return Surface();
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
            return Surface(points, triangles);
        }

        /**
         * All triangles are moved from source to target (unneccasary)
         */
        void replacePatch(std::size_t source, std::size_t target){
            for(auto triangle : m_triangles){
                if(triangle.patch_id() == source){
                    triangle.m_patch = target;
                }
            }
        }

        void flipPatchOrientation(std::size_t patch){
            //TODO
        }


    private:
        Surface(std::vector<VecFloat> points, std::vector<std::array<std::size_t, 3>> triangles, std::vector<std::size_t> patches)
        : m_points(points)
        , m_triangles()
        {
            // convert triangles
            std::size_t patch_counter = 0;
            m_triangles.reserve(triangles.size());
            for(std::size_t i = 0; i < triangles.size(); ++i){
                if(patch_counter < patches.size() && patches[patch_counter] == i){
                    ++patch_counter;
                }
                m_triangles.push_back(Triangle(triangles[i], patch_counter + 1));
            }
            // TODO: check if surface is valid
        }
        Surface(std::vector<VecFloat> points, std::vector<Triangle> triangles)
        : m_points(points)
        , m_triangles(triangles)
        {
        }
    private:
        std::vector<VecFloat> m_points;
        // contain the indices of the points
        std::vector<Triangle> m_triangles;
        std::vector<PatchInformation> m_patches;
    };

    // Keeps track of their points and their indices
    class SurfaceBuilder {
    public:
        SurfaceBuilder() {}
    };
}