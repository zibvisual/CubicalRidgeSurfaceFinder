#pragma once

#include <optional>

#include <utils/Lattice.hpp>
#include "Faces.hpp"

/**
 * @brief We assume to start from the middle of the voxel and go into a direction
 *
 * @param pStart
 * @return Face
 */
template <class Breaker>
std::optional<Face>
amanatides(VecInt pStart, VecInt pEnd, Lattice lattice, Breaker breaker)
{
    VecFloat tMax, tDelta;
    VecInt step;
    Dims dims = lattice.dims();
    VecFloat voxelSize = lattice.voxelsize();
    VecInt offset = VecInt(1, dims[0], dims[0] * dims[1]);
    VecInt stepId;

    for (int i = 0; i < 3; ++i)
    {
        int dir = pEnd[i] - pStart[i];
        if (dir > 0)
        {
            step[i] = 1;
            tDelta[i] = fabs(voxelSize[i] / static_cast<float>(dir));
            tMax[i] = tDelta[i] * 0.5;
        }
        else if (dir < 0)
        {
            step[i] = -1;
            tDelta[i] = fabs(voxelSize[i] / static_cast<float>(dir));
            tMax[i] = tDelta[i] * 0.5;
        }
        else
        {
            step[i] = 0;
            tDelta[i] = std::numeric_limits<float>::infinity();
            tMax[i] = std::numeric_limits<float>::infinity();
        }
        stepId[i] = step[i] * offset[i];
    }

    VecInt pNew = pStart;
    //TODO: this might be unsafe as pNew might not be inside lattice?
    std::size_t pNewId = lattice.c_index(pNew.tryVecSize().value());

    while (pNew != pEnd)
    {
        // float currDist = 0.f;
        if (tMax[0] < tMax[1])
        {
            if (tMax[0] < tMax[2])
            {
                pNew[0] = pNew[0] + step[0];
                pNewId += stepId[0];
                if (pNew[0] < 0 || pNew[0] >= dims[0] || breaker(pNewId))
                {
                    return step[0] < 0 ? Face(pNewId - stepId[0], Direction::LEFT) : Face(pNewId - stepId[0], Direction::RIGHT);
                }
                // currDist = tMax[0];
                tMax[0] = tMax[0] + tDelta[0];
            }
            else
            {
                pNew[2] = pNew[2] + step[2];
                pNewId += stepId[2];
                if (pNew[2] < 0 || pNew[2] >= dims[2] || breaker(pNewId))
                {
                    return step[2] < 0 ? Face(pNewId - stepId[2], Direction::BACKWARD) : Face(pNewId - stepId[2], Direction::FORWARD);
                }
                // currDist = tMax[2];
                tMax[2] = tMax[2] + tDelta[2];
            }
        }
        else
        {
            if (tMax[1] < tMax[2])
            {
                pNew[1] = pNew[1] + step[1];
                pNewId += stepId[1];
                if (pNew[1] < 0 || pNew[1] >= dims[1] || breaker(pNewId))
                {
                    return step[1] < 0 ? Face(pNewId - stepId[1], Direction::DOWN) : Face(pNewId - stepId[1], Direction::UP);
                }
                // currDist = tMax[1];
                tMax[1] = tMax[1] + tDelta[1];
            }
            else
            {
                pNew[2] = pNew[2] + step[2];
                pNewId += stepId[2];
                if (pNew[2] < 0 || pNew[2] >= dims[2] || breaker(pNewId))
                {
                    return step[2] < 0 ? Face(pNewId - stepId[2], Direction::BACKWARD) : Face(pNewId - stepId[2], Direction::FORWARD);
                }
                // currDist = tMax[2];
                tMax[2] = tMax[2] + tDelta[2];
            }
        }
    }

    return std::optional<Face>();
}