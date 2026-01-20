#include "src/rsf/CubicalRidgeSurfaceFinder.h"

extern "C" {

    struct VecFloat_C {
        float x;
        float y;
        float z;
    };

    struct Triangle_C {
    size_t v0;
    size_t v1;
    size_t v2;
    uint64_t patch;
};


    // struct Surface_C {
    //     VecFloat_C* points_begin;
    //     VecFloat_C* points_head;
    //     VecFloat_C* points_end;
    //     Triangle_C* triangles_begin;
    //     Triangle_C* triangles_head;
    //     Triangle_C* triangles_end;
    // };


    struct Surface_C {
        VecFloat_C* points;
        size_t num_points;
        size_t points_capacity;
        Triangle_C* triangles;
        size_t num_triangles;
        size_t triangles_capacity;
    };

    // function to transform the type of Surface to the struct
    Surface_C surface_to_c(const surface::Surface& s) {
        Surface_C out;

        // points
        out.points_capacity = s.points().capacity();
        out.num_points = s.points().size();
        out.points = new VecFloat_C[out.points_capacity];

        const auto& points = s.points();
        // out.points = reinterpret_cast<VecFloat_C*>(const_cast<VecFloat*>(points.data()));
        for (size_t i = 0; i < out.num_points; ++i) {
            out.points[i].x = points[i].x();
            out.points[i].y = points[i].y();
            out.points[i].z = points[i].z();
        }

        // triangles
        out.triangles_capacity = s.triangles().capacity();
        out.num_triangles = s.number_of_trianlges();
        out.triangles = new Triangle_C[out.triangles_capacity];

        const auto& triangles = s.triangles();
        // out.triangles = reinterpret_cast<Triangle_C*>(const_cast<surface::Triangle*>(triangles.data()));
        for (size_t i = 0; i < out.num_triangles; ++i) {
            out.triangles[i].v0 = triangles[i][0];
            out.triangles[i].v1 = triangles[i][1];
            out.triangles[i].v2 = triangles[i][2];
            out.triangles[i].patch = triangles[i].patch_id();
        }

        return out;
    }


    // function to free surface memory
    void free_surface(Surface_C* s) {
        if (!s) return;

        delete[] s->points;
        delete[] s->triangles;

        s->points = nullptr;
        s->triangles = nullptr;
        s->num_points = 0;
        s->num_triangles = 0;
    }

    // constructor
    void* CRSF_new() { // this must be created as a void pointer because python does not understand CubicalRidgeSurfaceFinder*
        auto* progress = new progressbar::ProgressbarReportDynamic(); // this is needed by the constructor
        auto* rsf = new ridgesurface::CubicalRidgeSurfaceFinder(*progress); // construct the object
        return reinterpret_cast<void*>(rsf); // the CubicalRidgeSurfaceFinder* is changed to a void pointer
    }

    // destructor
    void CRSF_delete(void* ptr) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr);
        delete rsf;
    }

    // set the input image
    void CRSF_setImage(void* ptr, float* data, int width, int height, int depth) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr);

        // Create a 1D vector that copies the Python array data
        std::vector<float> vec(data, data + width * height * depth);

        // Create a Rawfield
        RawField<float>* raw = new RawField<float>(vec, Dims(width, height, depth));

        // Set the image as an CRSF input
        rsf->setInput(raw);
    }

    // setters/getters
    void CRSF_setThreshold(void* ptr, float min, float max) {
        reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr)->setThresholds(min, max);
    }

    float CRSF_getMinThreshold(void* ptr) {
        return reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr)->getMin();
    }

    float CRSF_getMaxThreshold(void* ptr) {
        return reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr)->getMax();
    }

    // add seed points
    void CRSF_addSeed(void* ptr, float x, float y, float z, float distance) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr);

        ridgesurface::SeedPoint sp(VecFloat(x, y, z));
        ridgesurface::Seed seed(sp, distance);

        rsf->addSeed(seed);
    }

    // add automatic seed points
    int CRSF_addAutomaticSeeds(void* ptr, float minDistanceRel, float maxDistanceNewSeed) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr);
        return rsf->addSeed(minDistanceRel, maxDistanceNewSeed);
    }

    // clear seed points
    void CRSF_clearSeedPoints(void* ptr) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr);
        rsf->clearSeeds();
    }

    // compute
    void CRSF_compute(void* ptr) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr);
        rsf->compute();
    }

    // recalculate
    void CRSF_recalculate(void* ptr) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr);
        rsf->recalculate();
    }

// get patch surface
    Surface_C* CRSF_getPatchedSurface(void* ptr_rsf) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr_rsf);
        const surface::Surface& s = rsf->patchedSurface();
        // s.save("surfaceC.obj");

        Surface_C* out = new Surface_C;
        *out = surface_to_c(s);

        return out;
    }

    // finalize
    Surface_C* CRSF_finalize(void* ptr_rsf) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr_rsf);
        
        surface::Surface surface;
        rsf->finalize(&surface);

        // Convert to C struct
        Surface_C* out = new Surface_C;
        *out = surface_to_c(surface);

        return out;
    }

}

