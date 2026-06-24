#include <vector>
#include <variant>
#include <optional>
#include <io/argparse.hpp>
#include <io/spatialinfo.hpp>

namespace parsing {
    template <class T>
    std::optional<T> flatten_arg(const argparse::ArgumentParser& program, std::string name, T implicit_value){
        if(!program.is_used(name)){
            return {};
        }
        auto arg = program.get<std::vector<T>>(name);
        if(arg.empty()){
            return implicit_value;
        }
        return arg.front();
    }

    template <class T, unsigned long N>
    std::optional<std::array<T,N>> flatten_args(const argparse::ArgumentParser& program, std::string name, std::array<T, N> implicit_values)
    {
        if(!program.is_used(name)){
            return {};
        }
        auto result = implicit_values;
        auto arg = program.get<std::vector<T>>(name);
        std::size_t index = 0;
        while(index < arg.size()){
            result.push_back(arg[index++]);
        }
        return result;
    }

    template <class T>
    T get_arg(const argparse::ArgumentParser& program, std::string name, T implicit_value){
        auto arg = program.get<std::vector<T>>(name);
        if(arg.empty()){
            return implicit_value;
        }
        return arg.front();
    }

    template <class T, unsigned int N>
    std::vector<T> get_args(const argparse::ArgumentParser& program, std::string name, std::array<T, N> implicit_values)
    {
        auto result = std::vector<T>();
        auto arg = program.get<std::vector<T>>(name);
        std::size_t index = 0;
        while(index < arg.size()){
            result.push_back(arg[index++]);
        }
        while(index < N){
            result.push_back(implicit_values[index++]);
        }
        return result;
    }

    std::optional<SpatialArguments> merge_spatial_arguments(const argparse::ArgumentParser& program)
    {
        auto args = SpatialArguments();
        if(program.is_used("--bbox")){
            std::vector<float> vals = program.get<std::vector<float>>("--bbox");
            if(vals.size() == 2){
                args.bbox = GenericBBox(VecFloat(vals[0]), VecFloat(vals[1]));
            }else if(vals.size() == 6){
                args.bbox = GenericBBox(VecFloat(vals[0], vals[1], vals[2]), VecFloat(vals[3],vals[4],vals[5]));
            }else{
                //ERROR: we need either 2 or 6 float values
                std::cerr << "Either 2 or 6 float values are neccesary for a bounding box" << std::endl;
                return {};
            }
        }
        if(program.is_used("--origin")){
            std::vector<float> orig_vals = program.get<std::vector<float>>("--origin");
            if(orig_vals.size() == 1){
                args.origin = VecFloat(orig_vals[0]);
            }else if(orig_vals.size() == 3){
                args.origin = VecFloat(orig_vals[0], orig_vals[1], orig_vals[2]);
            }else{
                //ERROR: we need either 1 or 3 float values
                std::cerr << "Either 1 or 3 float values are neccesary for origin" << std::endl;
                return {};
            }
        }
        if(program.is_used("--voxelsize")){
            std::vector<float> voxel_vals = program.get<std::vector<float>>("--voxelsize");
            if(voxel_vals.size() == 1){
                args.voxelsize = VecFloat(voxel_vals[0]);
            }else if(voxel_vals.size() == 3){
                args.voxelsize = VecFloat(voxel_vals[0], voxel_vals[1], voxel_vals[2]);
            }else{
                //ERROR: we need either 1 or 3 float values
                std::cerr << "Either 1 or 3 float values are neccesary for voxelsize" << std::endl;
                return {};
            }
        }
        if(program.get<bool>("--corner")){
            args.bbox_signal = false;
            args.origin_signal = false;
        }
        if(program.get<bool>("--center")){
            args.bbox_signal = true;
            args.origin_signal = true;
        }
        return args;
    }
}