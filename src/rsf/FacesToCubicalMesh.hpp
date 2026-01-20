#pragma once

#include <unordered_map>
#include <utils/Dims.hpp>
#include <utils/BBox.hpp>
#include <utils/Lattice.hpp>
#include <surface/Surface.hpp>
#include <surface/StaticSurface.hpp>
#include "Faces.hpp"

class FacesToCubicalMesh
{
public:
    FacesToCubicalMesh()
        : m_voxelToMeshPoint()
    {}

    template <class Iter>
    void
    populateSurfacePatch(surface::Surface* surface, std::size_t patch, Lattice lattice, Iter begin, Iter end)
    {
        if (!surface)
            return;

        // add all faces and their corresponding points (use a hashmap to keep track of the points)
        uint64_t pointCounter = surface->points().size();

        Iter it = begin;
        while (it != end)
        {
            Face face = *it;
            auto corners = face.cornerIndices(lattice.dims());
            auto corner_positions = face.cornerPosition(lattice);
            for (std::size_t i = 0; i < 4; ++i)
            {
                const auto corner = corners[i];
                if (m_voxelToMeshPoint.find(corner) == m_voxelToMeshPoint.end())
                {
                    // m_voxelToMeshPoint does not contain point -> add it with a new index
                    m_voxelToMeshPoint[corner] = pointCounter++;
                    surface->addPoint(corner_positions[i]);
                }
            }
            // add the two triangles with the correct orientation and patch number
            surface->addTriangle(
                m_voxelToMeshPoint[corners[0]],
                m_voxelToMeshPoint[corners[1]],
                m_voxelToMeshPoint[corners[3]],
                patch
            );
            surface->addTriangle(
                m_voxelToMeshPoint[corners[1]],
                m_voxelToMeshPoint[corners[2]],
                m_voxelToMeshPoint[corners[3]],
                patch
            );
            ++it;
        }
    }

    template <class Iter>
    void
    populateSurface(surface::Surface* surface, Lattice lattice, Iter begin, Iter end)
    {
        if (!surface)
            return;

        surface->clear();
        populateSurfacePatch(surface, 1, lattice, begin, end);
    }

    template <class Iter>
    void
    populateSurface(surface::StaticSurface* surface, Lattice lattice, Iter begin, Iter end)
    {
        if (!surface)
            return;

        surface->clear();
        // add all faces and their corresponding points (use a hashmap to keep track of the points)
        uint64_t pointCounter = surface->points().size();

        Iter it = begin;
        while (it != end)
        {
            Face face = *it;
            auto corners = face.cornerIndices(lattice.dims());
            auto corner_positions = face.cornerPosition(lattice);
            for (std::size_t i = 0; i < 4; ++i)
            {
                const auto corner = corners[i];
                if (m_voxelToMeshPoint.find(corner) == m_voxelToMeshPoint.end())
                {
                    // m_voxelToMeshPoint does not contain point -> add it with a new index
                    m_voxelToMeshPoint[corner] = pointCounter++;
                    surface->addPoint(corner_positions[i]);
                }
            }
            // add the two triangles with the correct orientation
            surface->addTriangle(
                m_voxelToMeshPoint[corners[0]],
                m_voxelToMeshPoint[corners[1]],
                m_voxelToMeshPoint[corners[3]]
            );
            surface->addTriangle(
                m_voxelToMeshPoint[corners[1]],
                m_voxelToMeshPoint[corners[2]],
                m_voxelToMeshPoint[corners[3]]
            );
            ++it;
        }
    }

    // /**
    //  * @brief
    //  *
    //  * @tparam Iter
    //  * @param field
    //  * @param begin must be the same as populateSurface
    //  * @param end must be the same as populateSurface
    //  * @param data
    //  */
    // template <class Iter>
    // void populateSurfaceScalarField(HxSurfaceScalarField* field, Iter begin, Iter end, std::function<float(Face)> data)
    // {
    //     field->setEncoding(HxSurfaceField::Encoding::ON_TRIANGLES);
    //     Iter it = begin;
    //     std::size_t counter = 0;
    //     auto dataPtr = field->dataPtr();

    //     while (it != end)
    //     {
    //         Face face = *it;
    //         auto distance = data(face);
    //         dataPtr[counter++] = distance;
    //         dataPtr[counter++] = distance;
    //         ++it;
    //     }
    // }

    inline void clear(){
        m_voxelToMeshPoint.clear();
    }

protected:
    std::unordered_map<uint64_t, uint64_t> m_voxelToMeshPoint;
};