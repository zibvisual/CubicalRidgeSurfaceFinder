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

/**
 * Class to read and write spatial information.
 */
class SpatialInformation {
public:
    SpatialInformation():m_min(VecFloat(0.0f)),m_other(VecFloat(1.0f)){}
    SpatialInformation(VecFloat origin, VecFloat voxelsize, bool signal_origin = true):m_min(origin),m_other(voxelsize),m_signal(signal_origin){}
    SpatialInformation(CornerBBox bbox, bool signal_bbox = false):m_min(bbox.min_corner()),m_other(bbox.max_corner()),m_signal(signal_bbox),m_bbox(true){}

    static SpatialInformation fromArguments(std::optional<CornerBBox> bbox, std::optional<VecFloat> origin, std::optional<VecFloat> voxelsize, std::optional<bool> bbox_signal = false, std::optional<bool> origin_signal = true){
        // if we have origin and voxelsize, we use these
        if(origin.has_value() && voxelsize.has_value()){
            return SpatialInformation(origin.value(), voxelsize.value(), origin_signal.value_or(true));
        }

        // otherwise use the bbox
        if(bbox.has_value()){
            return SpatialInformation(bbox.value(), bbox_signal.value_or(false));
        }

        // otherwise use default voxelsize
        if(origin.has_value()){
            return SpatialInformation(origin.value(), VecFloat(1.0f), origin_signal.value_or(true));
        }

        // otherwise use default origin
        if(voxelsize.has_value()){
            return SpatialInformation(VecFloat(0.f), voxelsize.value());
        }

        return SpatialInformation();
    }

    static SpatialInformation fromMetadata(metadata_t meta){
        std::optional<CornerBBox> bbox;
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
                    bbox = CornerBBox(VecFloat(s2f(data[0])), VecFloat(s2f(data[1])));
                }else if(data.size() == 6){
                    bbox = CornerBBox(VecFloat(s2f(data[0]),s2f(data[1]),s2f(data[2])), VecFloat(s2f(data[3]),s2f(data[4]),s2f(data[5])));
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
        return SpatialInformation::fromArguments(bbox, origin, voxelsize, bbox_signal, origin_signal);
    }
    
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
    VecFloat m_min;
    VecFloat m_other;
    bool m_signal = true;
    bool m_bbox = false;
};