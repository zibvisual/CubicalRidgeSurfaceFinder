#pragma once

#include <surface/StaticSurface.hpp>

namespace surface {

struct PatchAddition {
    uint64_t patch_id;
    std::size_t point_start;
    std::size_t point_end;
    std::size_t triangle_start;
    std::size_t triangle_end;
};

/**
 * Represents the information necessary to update a patched surfaces.
 */
struct SurfaceUpdate {
    std::vector<PatchAddition> m_additions;
    std::vector<uint64_t> m_orientation_flips;
    std::vector<uint64_t> m_removals;
    std::vector<VecFloat> m_points;
    std::vector<SimpleTriangle> m_triangles;
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
        auto size = m_points.size();
        m_points.push_back(point);
        return size;
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