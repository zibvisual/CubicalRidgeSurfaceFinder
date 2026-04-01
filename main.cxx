/** Cubical Ridge Surface Finder (c) 2025 by Zuse Institute Berlin
 * 
 * Cubical Ridge Surface Finder is licensed under a
 * Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
 * 
 * You should have received a copy of the license along with this
 * work.  If not, see <https://creativecommons.org/licenses/by-nc-sa/4.0/>.
 */

#include <iostream>
#include "Config.h"
#include "src/utils/ProgressbarReportHelper.hpp"
#include "src/io/npy.hpp"
#include "src/io/argparse.hpp"
#include "src/field/RawField.hpp"

#include "src/rsf/CubicalRidgeSurfaceFinder.h"
#include "src/io/vertexset.hpp"

#include <utils/SeedSampler.hpp>

#include <filesystem>
#include <vector>

// #include <chrono>
// #include <thread>

// TODO: make seed input optional, and use automatic seeds and seed sampler if no seed input was set!
// TODO: then we need some seed threshold and maybe other parameters

namespace {
    std::optional<std::string> flatten_arg(const argparse::ArgumentParser& program, std::string name, std::string implicit){
        if(!program.is_used(name)){
            return {};
        }
        auto arg = program.get<std::vector<std::string>>(name);
        if(arg.empty()){
            return implicit;
        }
        return arg.front();
    }
}

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program("Ridge Surface Finder", std::to_string(RidgeSurfaceFinder_VERSION_MAJOR) + "." + std::to_string(RidgeSurfaceFinder_VERSION_MINOR));
    program.add_argument("input").help("Numpy input file to load");
    program.add_argument("-s", "--seed").default_value("0.01").help("Numpy label, wavefront object or float");
    program.add_argument("--min").default_value(0.f).scan<'g', float>().help("Min threshold. Everything under is not touched.");
    program.add_argument("--max").default_value(1.2f).scan<'g', float>().help("Max threshold. Everything above or equal has no cost.");
    program.add_argument("--range").default_value(50.0f).scan<'g', float>().help("How far each seedpoint should be marched for");
    program.add_argument("--corner").default_value(false).implicit_value(true).help("If the bounding box/origin is at the corner of the voxel instead of the center.");
    program.add_argument("--center").default_value(false).implicit_value(true).help("If the bounding box/origin is at the center of the voxel instead of the corner.");
    program.add_argument("--bbox").nargs(2,6).scan<'g', float>().help("The bounding box of the image");
    program.add_argument("--voxelsize").nargs(1,3).scan<'g', float>().help("The uniform voxel size");
    program.add_argument("--origin").nargs(1,3).scan<'g', float>().help("position of the image");
    program.add_argument("-o", "--output").help("Name of output file");
    program.add_argument("-d", "--debug").nargs(0,1).help("Save debug information in given path");
    program.add_argument("-a", "--automatic").nargs(1,2).help("Adds more seed points: (min-distance-between-points, [max seedpoints])");

    try
    {
        program.parse_args(argc, argv);
    }
    catch (const std::exception &err)
    {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    // Arguments
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
            return 1;
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
            return 1;
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
            return 1;
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

    // load input
    RawField<float> img;
    try
    {
        img = RawField<float>::load(program.get("input"));
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    std::cout << "loading seeds..." << std::endl;

    // load seeds (either labels or wavefront-like vertexset)
    std::vector<ridgesurface::Seed> seeds;
    try {
        // vertexset
        auto vertexset = read_vertexset(program.get("--seed"), program.get<float>("--range"));
        for(std::size_t i = 0; i < vertexset.points.size(); ++i){
            auto point = vertexset.points[i];
            seeds.push_back(ridgesurface::Seed(ridgesurface::SeedPoint(VecFloat(point.x(), point.y(), point.z())), vertexset.data[i]));
        }
        std::cout << "Loaded " << seeds.size() << " seeds from vertexset " << program.get("--seed") << std::endl;
    }
    catch (const std::exception &err) {
        // try labels
        try {
            auto labels = RawField<uint8_t>::load(program.get("--seed"));
            // generate seeds from voxels
            for(std::size_t i = 0; i < labels.lattice().size(); ++i){
                if(labels.get(RawField<uint8_t>::FieldLocation(i)) == 0){
                    seeds.push_back(ridgesurface::Seed(ridgesurface::SeedPoint(static_cast<VecFloat>(Lattice::gridLocationFromCIndex(i, labels.dims()))), program.get<float>("--range")));
                }
            }
            std::cout << "Loaded " << seeds.size() << " seeds from label field " << program.get("--seed") << std::endl;
        }
        catch (const std::exception &err) {
            // try integer conversion
            try {
                float threshold = std::stof(program.get("--seed"));
                // generate seeds through seed sampler
                std::vector<uint32_t> result = sample(img.data(), img.dims(), threshold);
    
                for (std::size_t i = 0; i < result.size(); ++i) {
                    VecSize gridLoc = Lattice::gridLocationFromCIndex(result[i], img.dims());
                    VecFloat worldPos = img.lattice().worldPosition(gridLoc);
                    seeds.push_back(ridgesurface::Seed(ridgesurface::SeedPoint(VecFloat(worldPos.x(), worldPos.y(), worldPos.z())), program.get<float>("--range")));   
                }
                std::cout << "Generated " << seeds.size() << " seeds" << std::endl;
            }catch (const std::exception &err) {
                std::cerr << "Neither paths to wavefront object or numpy labels found nor conversion to float possible." << std::endl;
                return 1;
            }
        }
    }

    // automatic values
    std::optional<float> automatic_threshold;
    std::optional<uint64_t> automatic_max_seedpoints;
    if(auto raw = program.present<std::vector<std::string>>("--automatic")){
        try {
            automatic_threshold = std::stof((*raw)[0]);
            if(raw->size() > 1){
                automatic_max_seedpoints = std::stoi((*raw)[1]);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return 1;
        }
    }

    // run
    try
    {
        std::filesystem::path input_name = program.get("input");
        auto debug_setting = flatten_arg(program, "-d", "debug");
        bool debug = debug_setting.has_value();
        std::filesystem::path debug_path = debug_setting.value_or("debug");
        auto patched_surface = surface::Surface();

        // warn if voxelsize are wideley different
        const auto voxelsize = img.getVoxelSize();
        const float max_voxelsize = std::max(std::max(voxelsize.x(), voxelsize.y()), voxelsize.z());
        const float min_voxelsize = std::min(std::min(voxelsize.x(), voxelsize.y()), voxelsize.z());
        if((max_voxelsize - min_voxelsize) / min_voxelsize > 0.5){
            std::cout << "Warning: Voxelsizes are assumed to be equal in all axis for the algorithms. Results might not be satisfactory." << std::endl;
        }

        const auto min = program.get<float>("--min");
        const auto max = program.get<float>("--max");

        std::cout << "Cost intensity range for marching set to [" << min << ", " << max << "]" << std::endl;

        auto progressbar = progressbar::Progressbar(progressbar::ProgressbarReportLevel::FullReport);
        ridgesurface::CubicalRidgeSurfaceFinder m_finder(progressbar);
        m_finder.setInput(&img.get_view());
        m_finder.newSeeds(seeds);
        m_finder.setThresholds(min,max);
        if(debug){
            m_finder.save_debug_information(debug_path / "seed_sample" / input_name.stem().concat("_rsf_"));
            patched_surface.update(m_finder.recalculate());
        }else{
            m_finder.recalculate();
        }

        if(automatic_threshold){
            // TODO support automatic progressbar (we do not know when we finish)
            progressbar.level(progressbar::ProgressbarReportLevel::FullReport);
            auto max = automatic_max_seedpoints.value_or(std::numeric_limits<uint64_t>::max());
            while(m_finder.numOfSeeds() < max){
                auto shift_distance = program.get<float>("--range") * 0.05;
                auto candidate = m_finder.getSeedpointCandidate(automatic_threshold.value(), shift_distance);
                if(!candidate.has_value()){
                    std::cout << "No more candidates found!" << std::endl;
                    break;
                }
                // seedpoint creation
                ridgesurface::SeedPoint sp(candidate.value());
                ridgesurface::Seed seed(sp, program.get<float>("--range"));
                m_finder.addSeed(seed);
                if(debug){
                    patched_surface.update(m_finder.calculate());
                }else{
                    m_finder.calculate();
                }
            }
        }
        
        auto surf = surface::StaticSurface();
        m_finder.finalize(&surf);
        std::filesystem::path default_output_name = input_name;
        default_output_name = default_output_name.replace_filename(input_name.stem().concat("_rsf")).replace_extension("obj");
        const auto output_name = program.present("-o").value_or(default_output_name);
        surf.save(output_name);
        std::cout << "Saved result in " << output_name << std::endl;

        // debug
        if(debug){
            m_finder.save_debug_information(debug_path / input_name.stem().concat("_rsf_"));
            patched_surface.save_each_patch(debug_path / input_name.stem());
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
}