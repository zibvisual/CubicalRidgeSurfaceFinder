#include "Seed.h"

namespace ridgesurface {
    SeedIterator
    Seed::getVoxelSources(Lattice lattice)
    {
        return SeedIterator(*this, lattice);
    }
}