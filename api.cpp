#include "src/rsf/CubicalRidgeSurfaceFinder.h"
#include "src/utils/SeedSampler.hpp"
#include <io/vertexset.hpp>

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
    surface::StaticSurface* surf_ptr;
};

// Class for Surface Updates
struct SurfaceUpdate_C {
    VecFloat_C* points;
    size_t num_points;
    Triangle_C* triangles;
    size_t num_triangles;
    // TODO: further fields
};

// Class for handling memory
struct CRSF_Handle {
    progressbar::Reporter progress;
    ridgesurface::CubicalRidgeSurfaceFinder rsf;
    RawFieldView<float> raw;
    surface::SurfaceUpdate update;

    CRSF_Handle()
    : progress(progressbar::Reporter())
    , rsf(ridgesurface::CubicalRidgeSurfaceFinder(progress))
    , raw((RawFieldView<float>(nullptr, Dims(0))))
    , update()
    {}

    // ~CRSF_Handle(){
    //     if(raw){
    //         delete raw;
    //     }
    // }
};

extern "C" {

    /**
     * @brief Transforms a `StaticSurface` object into a `Surface_C` type. 
     * @param s The `StaticSurface` object containing the points and triangles to transform.
     * @return A `Surface_C` object representing the same data from `StaticSurface`, 
     *         but in a C-compatible format.
     * 
     * Both objects point to the same memory location. Ownership belongs to `Surface_C` 
     * and `free_surface` must be called.
     */
    Surface_C surface_to_c(surface::StaticSurface* s) {
        Surface_C out; 
        out.num_points = s->points().size(); 
        out.points = reinterpret_cast<VecFloat_C*>(const_cast<VecFloat*>(s->points().data()));

        out.num_triangles = s->number_of_trianlges(); 
        out.triangles = reinterpret_cast<Triangle_C*>(const_cast<surface::SimpleTriangle*>(s->triangles().data())); 
        out.surf_ptr = s; 
        return out; 
    }

    /**
     * @brief Frees the memory allocated for a `Surface_C` object.
     * @param s A pointer to a `Surface_C` object whose memory will be freed.
     */
    void free_surface(Surface_C* s) {
        if (!s) return;

        if (s->surf_ptr) {
            delete s->surf_ptr;
            s->surf_ptr = nullptr;
        } 
        delete s;
    }

    /**
     * @brief Transforms a `SurfaceUpdate` object into a `SurfaceUpdate_C` type. 
     * @param s The `SurfaceUpdate` object.
     * @return A `SurfaceUpdate_C` object representing the same data from `SurfaceUpdate`, 
     *         but in a C-compatible format.
     * 
     * Both objects point to the same memory location. Ownership has the original `SurfaceUpdate`.
     */
    SurfaceUpdate_C surface_update_to_c(surface::SurfaceUpdate* update) {
        SurfaceUpdate_C out; 
        out.num_points = update->m_points.size();
        const auto& points = reinterpret_cast<VecFloat_C*>(update->m_points.data()); 

        out.num_triangles = update->m_triangles.size(); 
        out.triangles = reinterpret_cast<Triangle_C*>(update->m_triangles.data());  
        return out;
    }

    /**
     * @brief Constructor for `CubicalRidgeSurfaceFinder` and `CRSF_Handle` classes.
     * @return Pointer to the `CRSF_Handle` object.
     */ 
    CRSF_Handle* CRSF_new() { 
        return new CRSF_Handle();
    }

    /**
     * @brief Destructor for `CRSF_Handle` object.
     * @param h A pointer to the `CRSF_Handle` object.
     */ 
    void CRSF_delete(CRSF_Handle* h) {
        if (!h) return;
        delete h;
    }

    /**
     * @brief Set the input image for a `CubicalRidgeSurfaceFinder` object.
     * @param h A pointer to the `CRSF_Handle` object.
     * @param data A pointer to a 1D array containing the image data.
     * @param width, height, depth Dimensions of the input image (in voxels).
     */ 
    void CRSF_setImage(CRSF_Handle* h, float* data, int width, int height, int depth) {
        h->raw = RawFieldView<float>(data, Dims(width, height, depth));
        h->rsf.setInput(&h->raw);
    }

    /**
     * @brief Set the input image for a `CubicalRidgeSurfaceFinder` object.
     * @param h A pointer to the `CRSF_Handle` object.
     * @param data A pointer to a 1D array containing the image data.
     * @param width, height, depth Dimensions of the input image (in voxels).
     * @param origin Image origin point.
     * @param vsize Image voxel size.
     */ 
    void CRSF_setImage_bbox(CRSF_Handle* h, float* data, int width, int height, int depth, VecFloat_C origin_c, VecFloat_C vsize_c) {
        VecFloat origin(origin_c.x, origin_c.y, origin_c.z);
        VecFloat vsize(vsize_c.x, vsize_c.y, vsize_c.z);

        h->raw = RawFieldView<float>(data, Lattice(Dims(width, height, depth), origin, vsize));

        h->rsf.setInput(&h->raw);
    }

    /**
     * @brief Sets the minimum and maximum threshold for a `CubicalRidgeSurfaceFinder`.
     * @param h Pointer to the `CRSF_Handle` object.
     * @param min Minimum threshold value.
     * @param max Maximum threshold value.
     */
    void CRSF_setThreshold(CRSF_Handle* h, float min, float max) {
        h->rsf.setThresholds(min, max);
    }

    /**
     * @brief Gets the minimum threshold for a `CubicalRidgeSurfaceFinder.
     * @param h Pointer to the `CRSF_Handle` object.
     * @return The minimum threshold value.
     */
    float CRSF_getMinThreshold(CRSF_Handle* h) {
        return h->rsf.getMin();
    }

    /**
     * @brief Gets the maximum threshold for a `CubicalRidgeSurfaceFinder` object.
     * @param h Pointer to the `CRSF_Handle` object.
     * @return The maximum threshold value.
     */
    float CRSF_getMaxThreshold(CRSF_Handle* h) {
        return h->rsf.getMax();
    }

    /**
     * @brief Adds a seed point for the `CubicalRidgeSurfaceFinder`.
     * @param h Pointer to the `CRSF_Handle` object.
     * @param x, y, z The coordinates of the seed point in the 3D lattice.
     * @param distance The maximum distance to propagate during the Fast Marching algorithm.
     */
    void CRSF_addSeed(CRSF_Handle* h, float x, float y, float z, float distance) {

        ridgesurface::SeedPoint sp(VecFloat(x, y, z));
        ridgesurface::Seed seed(sp, distance);

        h->rsf.addSeed(seed);
    }

    /**
     * @brief Adds a new seed point based on the largest relative distance from candidate points.
     * 
     * This function iterates over all candidate seed points and selects the one with the largest 
     * relative distance to the previously existing seeds. The selected seed point is then added 
     * to the list of seed points for the surface finding algorithm.
     * 
     * @param h Pointer to the `CRSF_Handle` object.
     * @param minDistanceRel The minimum relative distance required to consider a candidate as a valid seed point.
     * @param maxDistanceNewSeed The maximum distance to propagate during the Fast Marching algorithm by the new seed.
     * 
     * @return The index of the seed point in the candidate list, or -1 if no valid seed was found. 
     *         A valid seed is one with a relative distance greater than `minDistanceRel`.
     */
    int CRSF_addAutomaticSeeds(CRSF_Handle* h, float minDistanceRel, float maxDistanceNewSeed) {
        auto candidate = h->rsf.getSeedpointCandidate(minDistanceRel, maxDistanceNewSeed*0.05);
        if(!candidate.has_value()){
            return -1;
        }
        // seedpoint creation
        ridgesurface::SeedPoint sp(candidate.value());
        ridgesurface::Seed seed(sp, maxDistanceNewSeed);
        return h->rsf.addSeed(seed);
    }

    /**
     * @brief Clears all seed points from a `CubicalRidgeSurfaceFinder` object.
     *
     * @param h Pointer to the `CRSF_Handle` object.
     */
    void CRSF_clearSeedPoints(CRSF_Handle* h) {
        h->rsf.clearSeeds();
    }

    /**
     * @brief Generates sample points to use as seed points in CRSF and adds them to the CRSF instance.
     * It chooses the voxel with the highest value for each conected graph.
     * @param h Pointer to the `CRSF_Handle` object.
     * @param data A pointer to a 1D array containing the image or distance map data.
     * @param threshold Threshold for a voxel to be considered as a candidate seed.
     * @param distance The maximum distance to propagate during the Fast Marching algorithm.
     */
    uint64_t CRSF_seedSampler(CRSF_Handle* h, float threshold, float distance)
    {
        if (h->rsf.data() == nullptr){
            throw std::runtime_error("CRSF_seedSampler error: no image set in CRSF.");
        }

        Dims dims = h->raw.getDims();
        auto lattice = h->raw.lattice();

        std::vector<uint32_t> result = sample(h->raw.data(), dims, threshold);

        for (std::size_t i = 0; i < result.size(); ++i) {
            VecSize gridLoc = Lattice::gridLocationFromCIndex(result[i], dims);
            VecFloat worldPos = lattice.worldPosition(gridLoc);
            float x = worldPos.x();
            float y = worldPos.y();
            float z = worldPos.z();
            CRSF_addSeed(h, x, y, z, distance);
        }
        
        return result.size();
    }


    /**
     * @brief Recalculates the ridge surface patches based on the current seed points.
     * 
     * This function clears any previous surface data, including seed points and surface patches,
     * and then recalculates the ridge surface patches based on the current seeds.
     * @param h Pointer to the `CRSF_Handle` object.
     * @return View to the UpdateSurface object. View is invalid if `calculate` or `recalculate` is called.
     */
    SurfaceUpdate_C CRSF_recalculate(CRSF_Handle* h) {
        h->update = h->rsf.recalculate();
        return surface_update_to_c(&h->update);
    }

    /**
     * @brief Computes the surface using Fast Marching algorithm;
     * @param h Pointer to the `CRSF_Handle` object.
     * @return View to the UpdateSurface object. View is invalid if `calculate` or `recalculate` is called.
     */
    SurfaceUpdate_C CRSF_calculate(CRSF_Handle* h) {
        h->update = h->rsf.calculate();
        return surface_update_to_c(&h->update);
    }

    /**
     * @brief Computes the surface using Fast Marching algorithm and merges all the patches.
     * @param h Pointer to the `CRSF_Handle` object.
     * @return A pointer to a `Surface_C` struct containing the final computed surface data.
     * 
     * Caller gets the ownership of the data and must call free_surface when the data is not needed anymore.
     */
    Surface_C* CRSF_finalize(CRSF_Handle* h) {
        h->rsf.calculate();
        // surface::StaticSurface surface;
        // rsf->finalize(&surface);

        // // Convert to a `Surface_C` struct
        // Surface_C* out = new Surface_C;
        // *out = surface_to_c(surface);
        
        surface::StaticSurface* surface = new surface::StaticSurface;
        h->rsf.finalize(surface);

        // Convert to a `Surface_C` struct
        Surface_C* out = new Surface_C;
        *out = surface_to_c(surface);

        return out;
    }

    /**
     * @brief Currently used for debugging purposes. Might be deleted at any time.
     * @param h Pointer to the `CRSF_Handle` object.
     * 
     * Saves all patches given by the current surface update. 
     * If called every time after recalculate/calculate, all generated patches are saved.
     */
    void CRSF_write_surface_update(CRSF_Handle* h, char* output_path) {
        auto filestem = std::string(output_path);
        h->update.save_each_patch(filestem);
    }

    /**
     * @brief Saves debug data.
     * @param h Pointer to the `CRSF_Handle` object.
     */
    void CRSF_save_debug_information(CRSF_Handle* h, char* output_path) 
    {
        auto op = std::string(output_path);
        h->rsf.save_debug_information(op);
    }
}
