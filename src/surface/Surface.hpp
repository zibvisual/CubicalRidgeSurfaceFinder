#pragma once

#include <array>
#include <vector>
#include <limits>
#include <iostream>
#include <unordered_set>

#include <surface/Triangle.hpp>
#include <surface/SurfaceUpdate.hpp>
#include <utils/Vec.hpp>
#include <utils/Permutation.hpp>
#include <io/common.hpp>
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
     * Patch Information. All triangles and points of a patch are in a subslice of the vectors.
     * All triangles and points to belong to exactly one patch. If a patch is deleted, the triangles and points can be so too.
     * The order of the patches and of the triangle slices and point slices stays the same!
     */
    struct PatchInformation {
        uint64_t id;
        bool orientation;
        bool deleted;
        std::size_t point_start;
        std::size_t point_end;
        std::size_t triangle_start;
        std::size_t triangle_end;

        PatchInformation(uint64_t id, std::size_t point_start, std::size_t point_end, std::size_t triangle_start, std::size_t triangle_end)
        : id(id)
        , orientation(false)
        , deleted(false)
        , point_start(point_start)
        , point_end(point_end)
        , triangle_start(triangle_start)
        , triangle_end(triangle_end) {}

        bool empty() const {
            return point_start == point_end;
        }
    };

    /**
     * A triangulated surface mesh with patches (triangle groups).
     * 
     * Invariant: patches are in a block within m_trianlges. That is, except for patch -1, all patches are contained in a range/slice.
     * Invariant: Each point and each triangle belong to exactly one patch.
     * Invariant: The indices of the triangles are relative to the patch start!
     * 
     */
    class Surface {
    public:
        Surface() : m_points(), m_triangles(), m_patches(), m_patchlocation(){}

        // single patch surface (id of 0)
        Surface(std::vector<VecFloat> points, std::vector<std::array<std::size_t, 3>> triangles)
        : m_points(points)
        , m_triangles()
        , m_patches()
        , m_patchlocation()
        {
            // convert triangles --> transmute should be possible
            for(auto triangle : triangles){
                m_triangles.push_back(SimpleTriangle(triangle));
            }
            m_patches.push_back(PatchInformation(0,0,m_points.size(), 0, m_triangles.size()));
            m_patchlocation[0] = 0;
        }

        static Surface load(std::filesystem::path input_path){
            auto data = read_wavefront(input_path);
            return Surface(data.points, data.triangles);
        }

        void save(std::filesystem::path output_path) const {
            write_wavefront(output_path, m_points, m_triangles);
        }

        void save_each_patch(std::filesystem::path path, bool save_empty_patches = false) const {
            for(const auto& patch : m_patches){
                if(patch.deleted || (!save_empty_patches && patch.empty())){
                    continue;
                }

                auto output = add_suffix(path, "_patch_" + std::to_string(patch.id));
                // iterarate and create spans
                auto points = std::span(m_points.cbegin() + patch.point_start, m_points.cbegin() + patch.point_end);
                auto triangles = std::span(m_triangles.cbegin() + patch.triangle_start, m_triangles.cbegin() + patch.triangle_end);
                write_wavefront(output, points, triangles);
            }
        }

        Surface clone() const {
            return Surface(
                m_points,
                m_triangles,
                m_patches,
                m_patchlocation
            );
        }

        void clean_data(){
            // TODO: do real deletion of points/triangles/patches
        }

        void removePatch(uint64_t id){
            m_patches[m_patchlocation[id]].deleted = true;
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
            return mutil::Permutation::mapPermutationNearestNeighbor(m_points, rhs.m_points, epsilon);
        }


        mutil::Permutation mapEqualTriangles(const Surface& rhs) const {
            auto perm = mutil::Permutation::mapPermutation(m_triangles, rhs.m_triangles, 
                [&](const SimpleTriangle& a, const SimpleTriangle& b) {
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

        /** we count */
        std::size_t number_of_patches() const {
            std::size_t counter = 0;
            for(const auto patch : m_patches){
                counter += static_cast<std::size_t>(!patch.deleted);
            }
            return counter;
        }

        // Methods to modify Surface. Usually not used directly.
        void clear() {
            m_points.clear();
            m_triangles.clear();
            m_patches.clear();
            m_patchlocation.clear();
        }

        /**
         * Should only be used when also adding another patch. Invariants must hold.
         */
        void addPoint(VecFloat point){
            m_points.push_back(point);
        }

        /**
         * Should only be used when also adding another patch. Invariants must hold.
         */
        void addTriangle(std::size_t point1, std::size_t point2, std::size_t point3){
            m_triangles.push_back(SimpleTriangle(point1, point2, point3));
        }

        void finishNewPatch(uint64_t id, std::size_t point_start, std::size_t triangle_start){
            m_patchlocation[id] = m_patches.size();
            m_patches.push_back(PatchInformation(id, point_start, m_points.size(), triangle_start, m_triangles.size()));
        }

        void flipPatchOrientation(std::size_t patch){
            // flip boolean
            m_patches[m_patchlocation[patch]].orientation ^= true;
        }

        void update(SurfaceUpdate update){
            for(auto deletion : update.m_removals){
                m_patches[m_patchlocation[deletion]].deleted = true;
            }
            
            // add points
            auto point_start = m_points.size();
            m_points.reserve(m_points.size() + update.m_points.size());
            m_points.insert(m_points.end(), update.m_points.begin(), update.m_points.end());

            // add triangles
            auto triangle_start = m_triangles.size();
            m_triangles.reserve(m_triangles.size() + update.m_triangles.size());
            m_triangles.insert(m_triangles.end(), update.m_triangles.begin(), update.m_triangles.end());

            // set patch id
            for(const auto& addition : update.m_additions){
                auto loc = m_patches.size();
                m_patches.push_back(
                    PatchInformation(
                        addition.patch_id, 
                        addition.point_start + point_start,
                        addition.point_end + point_start,
                        addition.triangle_start + triangle_start,
                        addition.triangle_end + triangle_start
                    ));
                m_patchlocation[addition.patch_id] = loc;
            }

            // orientation flips
            for(auto patch : update.m_orientation_flips){
                flipPatchOrientation(patch);
            }
            // TODO: if debug is activated, we keep a patched surface up to date.
        }


    private:
        /**
         * Only use permutation which permute points/triangles inside their own patch. Do not cross over!
         */
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
            return Surface(points, triangles, m_patches, m_patchlocation);
        }

        Surface(std::vector<VecFloat> points, std::vector<SimpleTriangle> triangles, std::vector<PatchInformation> patches)
        : m_points(points)
        , m_triangles(triangles)
        , m_patches(patches)
        , m_patchlocation()
        {
            // generate patch location
            std::size_t counter = 0;
            for(const auto& patch : m_patches){
                m_patchlocation[patch.id] = counter++; 
            }
        }

        Surface(std::vector<VecFloat> points, std::vector<SimpleTriangle> triangles, std::vector<PatchInformation> patches, std::unordered_map<uint64_t, std::size_t> patchlocation)
        : m_points(points)
        , m_triangles(triangles)
        , m_patches(patches)
        , m_patchlocation(patchlocation)
        {}

    private:
        std::vector<VecFloat> m_points;
        // contain the indices of the points
        std::vector<SimpleTriangle> m_triangles;
        std::vector<PatchInformation> m_patches;
        std::unordered_map<uint64_t, std::size_t> m_patchlocation;
    };

    // class SurfaceBuilder {
    // public:
    //     void startPatch(uint64_t patch_id)
    //     {
    //         m_current_patch = patch_id;
    //         m_first_point = m_points.size();
    //         m_first_triangle = m_triangles.size();
    //     }
    //     std::size_t addPoint(VecFloat point)
    //     {
    //         auto size = m_points.size();
    //         m_points.push_back(point);
    //         return size;
    //     }
    //     void addTriangle(std::size_t first, std::size_t second, std::size_t third)
    //     {
    //         m_triangles.push_back(SimpleTriangle(first, second, third));
    //     }
    //     void endPatch()
    //     {
    //         m_additions.push_back(
    //             PatchAddition {
    //                 m_current_patch,
    //                 m_first_point,
    //                 m_points.size(),
    //                 m_first_triangle,
    //                 m_triangles.size()
    //             }
    //         );
    //     }
    
    //     void addFlip(uint64_t patch_id){
    //         m_orientation_flips.push_back(patch_id);
    //     }
    
    //     void removePatch(uint64_t patch_id){
    //         m_removals.push_back(patch_id);
    //     }
    
    //     /**
    //      * Consumes current object and returns SurfaceUpdate. SurfaceUpdateBuilder is empty afterwards!
    //      */
    //     SurfaceUpdate build(){
    //         m_current_patch = 0;
    //         m_first_point = 0;
    //         m_first_triangle = 0;
    
    //         return SurfaceUpdate {
    //             std::move(m_additions),
    //             std::move(m_orientation_flips),
    //             std::move(m_removals),
    //             std::move(m_points),
    //             std::move(m_triangles)
    //         };
    //     }
    
    //     /**
    //      * Clear Builder without creating SurfaceUpdate
    //      */
    //     void clear(){
    //         m_additions.clear();
    //         m_orientation_flips.clear();
    //         m_removals.clear();
    //         m_points.clear();
    //         m_triangles.clear();
    
    //         m_current_patch = 0;
    //         m_first_point = 0;
    //         m_first_triangle = 0;
    //     }
    
    // private:
    //     std::vector<PatchAddition> m_additions;
    //     std::vector<uint64_t> m_orientation_flips;
    //     std::vector<uint64_t> m_removals;
    //     std::vector<VecFloat> m_points;
    //     std::vector<SimpleTriangle> m_triangles;
    
    //     uint64_t m_current_patch;
    //     std::size_t m_first_point;
    //     std::size_t m_first_triangle;
    // };
}