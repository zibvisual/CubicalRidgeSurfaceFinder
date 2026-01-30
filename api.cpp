#include "src/rsf/CubicalRidgeSurfaceFinder.h"

namespace ridgesurface {
    class CubicalRidgeSurfaceFinder;
}

// Class for points.
struct VecFloat_C {
    float x;
    float y;
    float z;
};

// Class for triangles.
struct Triangle_C {
    size_t v0;
    size_t v1;
    size_t v2;
};

// Class for surfaces.
struct Surface_C {
    VecFloat_C* points;
    size_t num_points;
    Triangle_C* triangles;
    size_t num_triangles;
};

extern "C" {

    /**
     * @brief Transforms a `StaticSurface` object into a `Surface_C` type. 
     * @param s The `StaticSurface` object containing the points and triangles to transform.
     * @return A `Surface_C` object representing the same data from `StaticSurface`, 
     *         but in a C-compatible format.
     */
    Surface_C surface_to_c(const surface::StaticSurface& s) {
        Surface_C out;

        out.num_points = s.points().size();
        out.points = new VecFloat_C[out.num_points];

        const auto& points = s.points();
        // out.points = reinterpret_cast<VecFloat_C*>(const_cast<VecFloat*>(points.data()));
        for (size_t i = 0; i < out.num_points; ++i) {
            out.points[i].x = points[i].x();
            out.points[i].y = points[i].y();
            out.points[i].z = points[i].z();
        }

        out.num_triangles = s.number_of_trianlges();
        out.triangles = new Triangle_C[out.num_triangles];

        const auto& triangles = s.triangles();
        // out.triangles = reinterpret_cast<Triangle_C*>(const_cast<surface::Triangle*>(triangles.data()));
        for (size_t i = 0; i < out.num_triangles; ++i) {
            out.triangles[i].v0 = triangles[i][0];
            out.triangles[i].v1 = triangles[i][1];
            out.triangles[i].v2 = triangles[i][2];
        }

        return out;
    }


    /**
     * @brief Frees the memory allocated for a `Surface_C` object.
     * @param s A pointer to a `Surface_C` object whose memory will be freed.
     */
    void free_surface(Surface_C* s) {
        if (!s) return;

        delete[] s->points;
        delete[] s->triangles;

        s->points = nullptr;
        s->triangles = nullptr;
        s->num_points = 0;
        s->num_triangles = 0;
    }

    /**
     * @brief Constructor for `CubicalRidgeSurfaceFinder` class .
     */ 
    ridgesurface::CubicalRidgeSurfaceFinder* CRSF_new() { 
        auto* progress = new progressbar::Progressbar(); // This is needed by the constructor
        auto* rsf = new ridgesurface::CubicalRidgeSurfaceFinder(*progress);
        return rsf;
    }

    /**
     * @brief Destructor for `CubicalRidgeSurfaceFinder` class.
     * @param rsf A pointer to the `CubicalRidgeSurfaceFinder` object.
     */ 
    void CRSF_delete(ridgesurface::CubicalRidgeSurfaceFinder* rsf) {
        delete rsf;
    }

    /**
     * @brief Set the input image for a `CubicalRidgeSurfaceFinder` object.
     * @param rsf A pointer to the `CubicalRidgeSurfaceFinder` object.
     * @param data A pointer to a 1D array containing the image data.
     * @param width, height, depth Dimensions of the input image (in voxels).
     */ 
    void CRSF_setImage(ridgesurface::CubicalRidgeSurfaceFinder* rsf, float* data, int width, int height, int depth) {

        // Create a 1D vector that copies the 'data' array.
        std::vector<float> vec(data, data + width * height * depth);

        // Create a Rawfield object with the image data and specified dimensions.
        RawField<float>* raw = new RawField<float>(vec, Dims(width, height, depth));

        // Set the RawField object as the input data.
        rsf->setInput(&raw->get_view());
    }

    /**
     * @brief Sets the minimum and maximum threshold for a `CubicalRidgeSurfaceFinder`.
     * @param rsf Pointer to the `CubicalRidgeSurfaceFinder` object.
     * @param min Minimum threshold value.
     * @param max Maximum threshold value.
     */
    void CRSF_setThreshold(ridgesurface::CubicalRidgeSurfaceFinder* rsf, float min, float max) {
        rsf->setThresholds(min, max);
    }

    /**
     * @brief Gets the minimum threshold for a `CubicalRidgeSurfaceFinder.
     * @param rsf Pointer to the `CubicalRidgeSurfaceFinder` object.
     * @return The minimum threshold value.
     */
    float CRSF_getMinThreshold(ridgesurface::CubicalRidgeSurfaceFinder* rsf) {
        return rsf->getMin();
    }
    /**
     * @brief Gets the maximum threshold for a `CubicalRidgeSurfaceFinder` object.
     * @param rsf Pointer to the `CubicalRidgeSurfaceFinder` object.
     * @return The maximum threshold value.
     */
    float CRSF_getMaxThreshold(ridgesurface::CubicalRidgeSurfaceFinder* rsf) {
        return rsf->getMax();
    }

    /**
     * @brief Adds a seed point for the `CubicalRidgeSurfaceFinder`.
     * @param rsf Pointer to the `CubicalRidgeSurfaceFinder` object.
     * @param x, y, z The coordinates of the seed point in the 3D lattice.
     * @param distance The maximum distance to propagate during the Fast Marching algorithm.
     */
    void CRSF_addSeed(ridgesurface::CubicalRidgeSurfaceFinder* rsf, float x, float y, float z, float distance) {

        ridgesurface::SeedPoint sp(VecFloat(x, y, z));
        ridgesurface::Seed seed(sp, distance);

        rsf->addSeed(seed);
    }

    /**
     * @brief Adds a new seed point based on the largest relative distance from candidate points.
     * 
     * This function iterates over all candidate seed points and selects the one with the largest 
     * relative distance to the previously existing seeds. The selected seed point is then added 
     * to the list of seed points for the surface finding algorithm.
     * 
     * @param minDistanceRel The minimum relative distance required to consider a candidate as a valid seed point.
     * @param maxDistanceNewSeed The maximum distance to propagate during the Fast Marching algorithm by the new seed.
     * 
     * @return The index of the seed point in the candidate list, or -1 if no valid seed was found. 
     *         A valid seed is one with a relative distance greater than `minDistanceRel`.
     */
    int CRSF_addAutomaticSeeds(ridgesurface::CubicalRidgeSurfaceFinder* rsf, float minDistanceRel, float maxDistanceNewSeed) {
        return rsf->addSeed(minDistanceRel, maxDistanceNewSeed).value_or(-1);
    }

    /**
     * @brief Clears all seed points from a `CubicalRidgeSurfaceFinder` object.
     *
     * @param rsf Pointer to the `CubicalRidgeSurfaceFinder` object.
     */
    void CRSF_clearSeedPoints(ridgesurface::CubicalRidgeSurfaceFinder* rsf) {
        rsf->clearSeeds();
    }

    /**
     * @brief Recalculates the ridge surface patches based on the current seed points.
     * 
     * This function clears any previous surface data, including seed points and surface patches,
     * and then recalculates the ridge surface patches based on the current seeds.
     * 
     * @return The number of processed seeds (either the last index if interrupted or the total number of seeds).
     */
    void CRSF_recalculate(ridgesurface::CubicalRidgeSurfaceFinder* rsf) {
        rsf->recalculate();
    }

    /**
     * @brief Computes the surface using Fast Marching algorithm and merges all the patches.
     * @param rsf Pointer to the `CubicalRidgeSurfaceFinder` object.
     * @return A pointer to a `Surface_C` struct containing the final computed surface data.
     */
    Surface_C* CRSF_compute_finalize(ridgesurface::CubicalRidgeSurfaceFinder* rsf) {
        rsf->calculate();
        surface::StaticSurface surface;
        rsf->finalize(&surface);

        // Convert to a `Surface_C` struct
        Surface_C* out = new Surface_C;
        *out = surface_to_c(surface);

        return out;
    }

}
