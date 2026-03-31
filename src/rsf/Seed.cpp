#include "Seed.h"

namespace ridgesurface {
    SeedVoxelSourcesIterator
    Seed::getVoxelSources(Lattice lattice)
    {
        return SeedVoxelSourcesIterator(*this, lattice);
    }
}