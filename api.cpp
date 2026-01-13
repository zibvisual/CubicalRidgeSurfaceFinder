#include "src/rsf/CubicalRidgeSurfaceFinder.h"

extern "C" {

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

    // finalize
    void CRSF_finalize(void* ptr_rsf, void* ptr_surf) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr_rsf);
        auto* surface = reinterpret_cast<surface::Surface*>(ptr_surf);

        rsf->finalize(surface);
    }

    // get patch surface
    void* CRSF_getPatchedSurface(void* ptr_rsf) {
        auto* rsf = reinterpret_cast<ridgesurface::CubicalRidgeSurfaceFinder*>(ptr_rsf);
        const surface::Surface& s = rsf->patchedSurface();
        return reinterpret_cast<void*>(new surface::Surface(s));
    }
}

