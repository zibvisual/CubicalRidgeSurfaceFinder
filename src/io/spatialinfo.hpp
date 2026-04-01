#pragma once

#include <optional>
#include <utils/Vec.hpp>
#include <utils/BBox.hpp>
#include <utils/Lattice.hpp>
#include <io/metadata.hpp>
#include <fmt/format.h>

namespace {
    // more strict parsing
    inline float s2f(const std::string& str){
        std::size_t pos;
        float val = stof(str, &pos);
        if(pos != str.length()){
            throw std::invalid_argument(fmt::format("Error: String {} could not be converted to a float!", str));
        }
        return val;
    }
}

struct SpatialArguments;

/**
 * Class to read and write spatial information.
 */
class SpatialInformation {
public:
    SpatialInformation():m_min(VecFloat(0.0f)),m_other(VecFloat(1.0f)){}
    SpatialInformation(VecFloat origin, VecFloat voxelsize, bool signal_origin = true):m_min(origin),m_other(voxelsize),m_signal(signal_origin){}
    SpatialInformation(CornerBBox bbox):m_min(bbox.min_corner()),m_other(bbox.max_corner()),m_signal(false),m_bbox(true){}
    SpatialInformation(CenterBBox bbox):m_min(bbox.min_center()),m_other(bbox.max_center()),m_signal(true),m_bbox(true){}

    static SpatialInformation fromLattice(const Lattice& lattice)
    {
        return SpatialInformation(lattice.origin(), lattice.voxelsize());
    }
    
    /**
     * Lattice always contains the centered origin and voxelsize.
     */
    Lattice toLattice(Dims dims){
        // first calculate the voxelsize
        VecFloat voxelsize;
        if(m_bbox && m_signal){
            voxelsize = VecFloat(
                dims[0] == 1 ? 0 : (m_other[0] - m_min[0]) / static_cast<float>(dims[0] - 1),
                dims[1] == 1 ? 0 : (m_other[1] - m_min[1]) / static_cast<float>(dims[1] - 1),
                dims[2] == 1 ? 0 : (m_other[2] - m_min[2]) / static_cast<float>(dims[2] - 1)
            );
        }else if(m_bbox && !m_signal){
            voxelsize = (m_other - m_min) / static_cast<VecFloat>(dims);
        }else{
            voxelsize = m_other;
        }
        // now calculate the origin
        VecFloat origin = m_signal ? m_min : m_min + (voxelsize * 0.5f);
        return Lattice(dims, origin, voxelsize);
    }

private:
    friend SpatialArguments;

    VecFloat m_min;
    // either max or voxelsize (if m_bbox true -> max)
    VecFloat m_other;
    // if the m_min and m_other represent the center or the corner of the voxel
    bool m_signal = true;
    bool m_bbox = false;
};

struct SpatialArguments
{
    SpatialArguments(){}
    SpatialArguments(std::optional<GenericBBox> bbox, std::optional<VecFloat> origin, std::optional<VecFloat> voxelsize, std::optional<bool> bbox_signal, std::optional<bool> origin_signal)
    : bbox(bbox)
    , origin(origin)
    , voxelsize(voxelsize)
    , bbox_signal(bbox_signal)
    , origin_signal(origin_signal)
    {}

    SpatialArguments& update(const SpatialArguments& args)
    {
        if(args.bbox)
            bbox = args.bbox;
        if(args.origin)
            origin = args.origin;
        if(args.voxelsize)
            voxelsize = args.voxelsize;
        if(args.bbox_signal)
            bbox_signal = args.bbox_signal;
        if(args.origin_signal)
            origin_signal = args.origin_signal;
        return *this;
    }

    bool empty() const
    {
        return !(bbox || origin || voxelsize || bbox_signal || origin_signal); 
    }

    /**
     * Return a spatial argument which would return the given spatial informatin.
     * 
     * As toInformation() is a lossy process, the returned SpatialArgument might be different to the one which created the input.
     */
    static SpatialArguments fromInformation(const SpatialInformation& info)
    {
        if(info.m_bbox){
            return SpatialArguments(GenericBBox(info.m_min, info.m_other), std::nullopt, std::nullopt, info.m_signal, std::nullopt);
        }else{
            return SpatialArguments(std::nullopt, info.m_min, info.m_other, std::nullopt, info.m_signal);
        }
    }

    /**
     * Generate SpatialInformatin from optional arguments given by bbox, origin, voxelsize.
     * 
     * bbox_signal should be true if the user says that the bounding box start and ends at the center of voxels (which is the case in amira)
     * origin_signal should be true if the user says that the given origin is the center of the first voxel
     */
    SpatialInformation toInformation()
    {
        // use bbox if possible (default being the voxel corners (no signal))
        if(bbox.has_value()){
            if (bbox_signal.value_or(false)){
                return SpatialInformation((CenterBBox) bbox.value());
            }else{
                return SpatialInformation((CornerBBox) bbox.value());   
            }
        }

        // otherwise we use origin and voxelsize
        return SpatialInformation(origin.value_or(VecFloat(0.f)), voxelsize.value_or(VecFloat(1.0f)),origin_signal.value_or(true));
    }

    static SpatialArguments fromMetadata(metadata_t meta){
        std::optional<GenericBBox> bbox;
        std::optional<VecFloat> origin;
        std::optional<VecFloat> voxelsize;
        std::optional<bool> bbox_signal;
        std::optional<bool> origin_signal;

        const auto size = std::min(meta.keywords.size(), meta.data.size());
        for(std::size_t i = 0; i < size; ++i){
            const auto keyword = meta.keywords[i];
            const auto data = meta.data[i];

            // bbox
            if((keyword == "bbox" || keyword == "center-bbox" || keyword == "corner-bbox") && !bbox.has_value()){
                if(data.size() == 2){
                    bbox = GenericBBox(VecFloat(s2f(data[0])), VecFloat(s2f(data[1])));
                }else if(data.size() == 6){
                    bbox = GenericBBox(VecFloat(s2f(data[0]),s2f(data[1]),s2f(data[2])), VecFloat(s2f(data[3]),s2f(data[4]),s2f(data[5])));
                }
            }
            // corner or center?
            if(keyword == "corner-bbox" && !bbox_signal.has_value()){
                bbox_signal = false;
            }else if(keyword == "center-bbox" && !bbox_signal.has_value()){
                bbox_signal = true;
            }

            // origin
            if((keyword == "origin" || keyword == "center" || keyword == "corner")  && !origin.has_value()){
                if(data.size() == 1){
                    origin = VecFloat(s2f(data[0]));
                }else if(data.size() == 3){
                    origin = VecFloat(s2f(data[0]),s2f(data[1]),s2f(data[2]));
                }
            }
            // corner or center?
            if(keyword == "corner" && !origin_signal.has_value()){
                origin_signal = false;
            }else if(keyword == "center" && !origin_signal.has_value()){
                origin_signal = true;
            }

            // voxelsize
            if(keyword == "voxelsize" && !voxelsize.has_value()){
                if(data.size() == 1){
                    voxelsize = VecFloat(s2f(data[0]));
                }else if(data.size() == 3){
                    voxelsize = VecFloat(s2f(data[0]),s2f(data[1]),s2f(data[2]));
                }
            }
        }
        return SpatialArguments(bbox, origin, voxelsize, bbox_signal, origin_signal);
    }

    std::optional<GenericBBox> bbox;
    std::optional<VecFloat> origin;
    std::optional<VecFloat> voxelsize;
    std::optional<bool> bbox_signal;
    std::optional<bool> origin_signal;
};