#pragma once

#include <utils/Lattice.hpp>

namespace ridgesurface {

struct CubicalRidgeSurfaceFinderSettings
{
    float image_border_threshold = 0.1;

    float border_distance_threshold(const Lattice& lattice, float distance_to_march) const
    {
        return lattice.voxelsize().length() * std::min(distance_to_march, lattice.centerbox().min()) * image_border_threshold;
    }
};

} // namespace ridgesurface