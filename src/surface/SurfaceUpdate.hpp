#pragma once

#include <surface/StaticSurface.hpp>

namespace surface {

struct PatchAddition {
    uint64_t patch_id;
    std::size_t point_start;
    std::size_t point_end;
    std::size_t triangle_start;
    std::size_t triangle_end;

    bool empty() const {
        return point_start == point_end;
    }
};

/**
 * Represents the information necessary to update a patched surfaces.
 * 
 * Invariant: The indices of the triangles are relative to the patch start!
 */
struct SurfaceUpdate {
    std::vector<PatchAddition> m_additions;
    std::vector<uint64_t> m_orientation_flips;
    std::vector<uint64_t> m_removals;
    std::vector<VecFloat> m_points;
    std::vector<SimpleTriangle> m_triangles;

    // SurfaceUpdate() : m_additions(), m_orientation_flips(), m_removals(), m_points(), m_triangles(){}

    void save_each_patch(std::filesystem::path path, bool save_empty_patches = false) const {
        for(const auto& patch : m_additions){
            if(!save_empty_patches && patch.empty()){
                continue;
            }
            auto output = add_suffix(path, "_patch_" + std::to_string(patch.patch_id));

            // iterarate and create spans
            auto points = std::span(m_points.cbegin() + patch.point_start, m_points.cbegin() + patch.point_end);
            auto triangles = std::span(m_triangles.cbegin() + patch.triangle_start, m_triangles.cbegin() + patch.triangle_end);
            write_wavefront(output, points, triangles);
        }
    }
};

class SurfaceUpdateBuilder {
public:
    void startPatch(uint64_t patch_id)
    {
        m_current_patch = patch_id;
        m_first_point = m_points.size();
        m_first_triangle = m_triangles.size();
    }
    std::size_t addPoint(VecFloat point)
    {
        auto relative_offset = m_points.size() - m_first_point;
        m_points.push_back(point);
        return relative_offset;
    }
    void addTriangle(std::size_t first, std::size_t second, std::size_t third)
    {
        m_triangles.push_back(SimpleTriangle(first, second, third));
    }
    void endPatch()
    {
        m_additions.push_back(
            PatchAddition {
                m_current_patch,
                m_first_point,
                m_points.size(),
                m_first_triangle,
                m_triangles.size()
            }
        );
    }

    void addFlip(uint64_t patch_id){
        m_orientation_flips.push_back(patch_id);
    }

    void removePatch(uint64_t patch_id){
        m_removals.push_back(patch_id);
    }

    /**
     * Consumes current object and returns SurfaceUpdate. SurfaceUpdateBuilder is empty afterwards!
     */
    SurfaceUpdate build(){
        m_current_patch = 0;
        m_first_point = 0;
        m_first_triangle = 0;

        return SurfaceUpdate {
            std::move(m_additions),
            std::move(m_orientation_flips),
            std::move(m_removals),
            std::move(m_points),
            std::move(m_triangles)
        };
    }

    /**
     * Clear Builder without creating SurfaceUpdate
     */
    void clear(){
        m_additions.clear();
        m_orientation_flips.clear();
        m_removals.clear();
        m_points.clear();
        m_triangles.clear();

        m_current_patch = 0;
        m_first_point = 0;
        m_first_triangle = 0;
    }

private:
    std::vector<PatchAddition> m_additions;
    std::vector<uint64_t> m_orientation_flips;
    std::vector<uint64_t> m_removals;
    std::vector<VecFloat> m_points;
    std::vector<SimpleTriangle> m_triangles;

    uint64_t m_current_patch;
    std::size_t m_first_point;
    std::size_t m_first_triangle;
};

} // namespace surface