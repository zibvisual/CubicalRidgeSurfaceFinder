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
#include <utils/StructureTensor.hpp>

#include <field/GridPositionBorderIterator.hpp>

#include <filesystem>
#include <vector>
#include <variant>

// #include <chrono>
// #include <thread>


namespace {
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

    void saveSeeds(std::filesystem::path op, const std::vector<ridgesurface::Seed>& seeds){
        auto data = vertex_data<float>{std::vector<VecFloat>(), std::vector<float>()};
        for(auto seed : seeds){
            data.points.push_back(seed.firstPoint());
            data.data.push_back(seed.getDistance());
        }
        write_vertexset(op, data);
    }
}

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program("Ridge Surface Finder", std::to_string(RidgeSurfaceFinder_VERSION_MAJOR) + "." + std::to_string(RidgeSurfaceFinder_VERSION_MINOR));
    program.add_argument("input").help("Input file to load");

    // seeds
    auto& input = program.add_mutually_exclusive_group(true);
    input.add_argument("-s", "--seed").help("Label field or vertex set");
    input.add_argument("--sample-component").default_value(0.01f).scan<'g', float>().help("Generate one seed point per component (threshold creates mask). Useful with --automatic");
    input.add_argument("--sample-greedy").default_value(40.0f).scan<'g', float>().help("Generate many seedpoint via FM, given by radius.");

    // fm settings
    program.add_argument("--min").default_value(0.f).scan<'g', float>().help("Min threshold. Everything under is not touched.");
    program.add_argument("--max").default_value(1.2f).scan<'g', float>().help("Max threshold. Everything above or equal has no cost.");
    program.add_argument("--range").default_value(50.0f).scan<'g', float>().help("How far each seedpoint should be marched for");
    
    // processing
    program.add_argument("-a", "--automatic").nargs(1,2).help("Adds more seed points: (min-distance-between-points, [max seedpoints])");
    program.add_argument("--border-margin").default_value(0.1f).scan<'g', float>().help("Threshold when to skip given seed points. Value between 0 and 1.");
    program.add_argument("--shift").nargs(0,1).scan<'g', float>().help("Distance how far a seed point is allowed to shift to go to the ridge.");
    program.add_argument("--border-padding").nargs(0,1).scan<'i', int>().help("How much the image should be extended. Given by voxel number.");

    // spatial information
    auto &centering = program.add_mutually_exclusive_group();
    centering.add_argument("--corner").default_value(false).implicit_value(true).help("If the bounding box/origin is at the corner of the voxel instead of the center.");
    centering.add_argument("--center").default_value(false).implicit_value(true).help("If the bounding box/origin is at the center of the voxel instead of the corner.");
    program.add_argument("--bbox").nargs(2,6).scan<'g', float>().help("The bounding box of the image");
    program.add_argument("--voxelsize").nargs(1,3).scan<'g', float>().help("The uniform voxel size");
    program.add_argument("--origin").nargs(1,3).scan<'g', float>().help("position of the image");
    
    // filter
    // program.add_argument("-f", "--filter");

    // output
    program.add_argument("-o", "--output").help("Name of output file");
    program.add_argument("-d", "--debug").nargs(0,1).help("Save debug information in given path");

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

    auto time_begin = std::chrono::steady_clock::now();
    auto reporter = progressbar::Reporter();
    // reporter.start("CRSF");

    std::filesystem::path input_name = program.get("input");
    auto debug_setting = flatten_arg(program, "-d", "debug");
    bool debug = debug_setting.has_value();
    std::filesystem::path debug_path = debug_setting.value_or("debug");

    const auto min = program.get<float>("--min");
    const auto max = program.get<float>("--max");

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
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Load " << program.get("input") << std::endl;
    RawField<float> img;
    try
    {
        img = RawField<float>::load(program.get("input"), args);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    std::cout << "Dims: " << img.dims() << std::endl;

    // load seeds (either labels or wavefront-like vertexset)
    std::vector<ridgesurface::Seed> seeds;
    if(program.is_used("--seed")){
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
                std::cerr << "Neither paths to wavefront object nor numpy labels found." << std::endl;
                return 1;
            }
        }
    }

    // generate seeds
    if(program.is_used("--sample-greedy") || program.is_used("--sample-component"))
    {
        std::cout << "Generate seeds" << std::endl;
        std::vector<uint32_t> generated_seed_indices;
        if(program.is_used("--sample-greedy"))
        {
            generated_seed_indices = greedySampling(reporter, img.get_view(), program.get<float>("--sample-greedy"), min, max);
        }else if(program.is_used("--sample-component"))
        {
            float threshold = program.get<float>("--sample-component");
            generated_seed_indices = sample(img.data(), img.dims(), threshold);
        }
    
        for (std::size_t i = 0; i < generated_seed_indices.size(); ++i) {
            VecSize gridLoc = Lattice::gridLocationFromCIndex(generated_seed_indices[i], img.dims());
            VecFloat worldPos = img.lattice().worldPosition(gridLoc);
            seeds.push_back(ridgesurface::Seed(ridgesurface::SeedPoint(VecFloat(worldPos.x(), worldPos.y(), worldPos.z())), program.get<float>("--range")));   
        }
        generated_seed_indices.clear();
        std::cout << "Generated " << seeds.size() << " seeds" << std::endl;

        if(debug){
            saveSeeds(debug_path / input_name.stem() / input_name.stem().concat("_rsf_generatedSeeds.obj"), seeds);
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

    // extend image
    Dims originalDims = img.lattice().dims();
    if(program.is_used("--border-padding")){
        int padding = program.get<int>("--border-padding");
        std::cout << "Extend image by " << padding << " voxels" << std::endl;
        img = img.extend(padding);
        // fade out image
        for(auto pair : GridPositionBorderIterator(img.dims(), padding))
        {
            auto val = img.get_unchecked(pair.first);
            img.set_unchecked(pair.first, val * pair.second);
        }

        if(debug){
            img.save(debug_path / input_name.stem() / input_name.stem().concat("_rsf_extendedField.nrrd"));
        }

    }

    // run
    try
    {
        auto patched_surface = surface::Surface();

        // warn if voxelsize are wideley different
        const auto voxelsize = img.getVoxelSize();
        const float max_voxelsize = std::max(std::max(voxelsize.x(), voxelsize.y()), voxelsize.z());
        const float min_voxelsize = std::min(std::min(voxelsize.x(), voxelsize.y()), voxelsize.z());
        if((max_voxelsize - min_voxelsize) / min_voxelsize > 0.5){
            std::cout << "Warning: Voxelsizes are assumed to be equal in all axis for the algorithms. Results might not be satisfactory." << std::endl;
        }

        std::cout << "Cost intensity range for marching set to [" << min << ", " << max << "]" << std::endl;

        ridgesurface::CubicalRidgeSurfaceFinder m_finder(reporter);
        m_finder.setInput(&img.get_view());
        m_finder.setThresholds(min,max);
        if(program.is_used("--border-padding")){
            VecSize min_clip = VecSize(program.get<int>("--border-padding"));
            VecSize max_clip = static_cast<VecSize>(originalDims) + min_clip;
            m_finder.getTransformer().setClipping(min_clip, max_clip);
        }        

        const float default_shift_distance = program.get<float>("--range") * 0.1;
        const std::optional<float> shift_distance = flatten_arg<float>(program, "shift", default_shift_distance);

        // shift seeds
        if(shift_distance){
            const float shift_distance_val = shift_distance.value() * img.getVoxelSize().length();

            std::cout << "Shift seeds by a maximum of " << shift_distance_val  << " distance" << std::endl;
            if(debug){
                // record structure tensor vector via lineset
                auto builder = LineSetBuilder<std::monostate>();
                const float gradient_sigma = 3.0f;
                const float tensor_sigma = 2.0f;
                for(const auto& seed: seeds)
                {
                    auto gridPoint = img.lattice().gridLocation(seed.firstPoint());
                    VecFloat vector = structure_tensor_direction(img.data(), img.dims(), static_cast<VecInt>(gridPoint), gradient_sigma, tensor_sigma);
                    // project a vector onto 3*voxelsize
                    VecFloat projected = vector.normalize().value_or(VecFloat(0.0f)) * img.getVoxelSize() * 10;
                    builder.add_point(seed.firstPoint() - projected, {});
                    builder.add_point(seed.firstPoint(), {});
                    builder.add_point(seed.firstPoint() + projected, {});
                    builder.push_line();
                }
                builder.build().save(debug_path / input_name.stem() / input_name.stem().concat("_rsf_structure_tensor_direction.sls"));

                // record shift via lineset
                auto builder2 = LineSetBuilder<int8_t>();
                builder2.reserve_lines(seeds.size());

                for(auto& seed: seeds){
                    builder2.add_point(seed.firstPoint(), 0);
                    m_finder.moveSeedToRidge(seed, shift_distance_val);
                    builder2.add_point(seed.firstPoint(), 1);
                    builder2.push_line();
                }
                const auto lineset = builder2.build();
                lineset.save(debug_path / input_name.stem() / input_name.stem().concat("_rsf_shiftlines.sls"));
                saveSeeds(debug_path / input_name.stem() / input_name.stem().concat("_rsf_shiftedSeeds.obj"), seeds);
            }else{
                // just shift
                reporter.start("Shifting all seeds");
                for(auto& seed: seeds){
                    m_finder.moveSeedToRidge(seed, shift_distance_val);
                }
                reporter.end();
            }
        }
        
        // skip invalid seeds (only if not sampled by components?)
        if(program.is_used("--border-margin"))
        {
            m_finder.settings.image_border_threshold = program.get<float>("--border-margin");
            auto skippedSeeds = std::vector<ridgesurface::Seed>();
            for(auto seed: seeds){
                if(m_finder.validSeed(seed)){
                    m_finder.addSeed(seed);
                }else{
                    skippedSeeds.push_back(seed);
                }
            }
            if(skippedSeeds.size() > 0){
                std::cout << "Skipped " << skippedSeeds.size() << " seeds" << std::endl;
                if(debug)
                {
                    saveSeeds(debug_path / input_name.stem() / input_name.stem().concat("_rsf_skippedSeeds.obj"), skippedSeeds);
                }
            }
        }else{
            for(auto seed: seeds){
                m_finder.addSeed(seed);
            }
        }
        
        if(debug){
            patched_surface.update(m_finder.recalculate());
        }else{
            m_finder.recalculate();
        }

        if(automatic_threshold){
            auto max = automatic_max_seedpoints.value_or(std::numeric_limits<uint64_t>::max());
            while(m_finder.numOfSeeds() < max){
                auto candidate = m_finder.getSeedpointCandidate(automatic_threshold.value(), shift_distance.value_or(default_shift_distance));
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

        std::cout << "Number of seeds actively used: " << m_finder.numOfSeeds() << std::endl;
        
        auto surf = surface::StaticSurface();
        m_finder.finalize(&surf);
        std::filesystem::path default_output_name = input_name.parent_path() / input_name.stem().concat("_rsf_finalized.obj");
        const auto output_name = program.present("-o").value_or(default_output_name);
        surf.save(output_name);
        std::cout << "Saved result in " << output_name << std::endl;

        // debug
        if(debug){
            m_finder.save_debug_information(debug_path / input_name.stem() / input_name.stem().concat("_rsf_"));
            patched_surface.save_each_patch(debug_path / input_name.stem() / "patches" / input_name.stem().concat("_rsf_"));
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    // reporter.end();
    auto time_end = std::chrono::steady_clock::now();
    std::cout << "CRSF finished! [" << (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count() / 1000000.0) << " sec]" << std::endl;
    return 0;
}