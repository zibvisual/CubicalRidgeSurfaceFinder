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
 #include <utils/parsing.hpp>
 
 #include <field/GridPositionBorderIterator.hpp>
 
 #include <filesystem>
 #include <vector>
 #include <variant>

void saveSeeds(std::filesystem::path op, const std::vector<ridgesurface::Seed>& seeds){
    auto data = vertex_data<float>{std::vector<VecFloat>(), std::vector<float>()};
    for(auto seed : seeds){
        data.points.push_back(seed.firstPoint());
        data.data.push_back(seed.getDistance());
    }
    write_vertexset(op, data);
}
 
 // #include <chrono>
 // #include <thread>
 
int main(int argc, char *argv[])
{
    argparse::ArgumentParser program("Membrane MidSurface Extractor", std::to_string(RidgeSurfaceFinder_VERSION_MAJOR) + "." + std::to_string(RidgeSurfaceFinder_VERSION_MINOR));
    program.add_argument("input").help("Input file to load");

    // seed generation
    program.add_argument("algorithm").default_value(std::string{"shell"}).choices("fm", "shell");
    program.add_argument("--seed-range").default_value(0.8f).scan<'g', float>().help("Exclusion zone dependent on fast marching range and given multiplicative.");
    program.add_argument("--seed-threshold").scan<'g', float>().help("Min threshold. Everything under is not touched. [default: dynamic]");
    program.add_argument("--seed-shift").nargs(0,1).scan<'g', float>().help("Distance how far a seed point is allowed to shift to go to the ridge. [implicit: --range * 0.1]");

    // fm settings
    program.add_argument("--min").scan<'g', float>().help("Min threshold of FM. Everything under is not touched. [default: min of image]");
    program.add_argument("--max").scan<'g', float>().help("Max threshold of FM. Everything above or equal has no cost. [default: max of image * 1.2 (max - min)]");
    program.add_argument("--range").default_value(50.0f).scan<'g', float>().help("How far each seedpoint should be marched for");
     
    // processing
    program.add_argument("--border-padding").default_value(10).scan<'i', int>().help("How much the image should be extended. Given by voxel number.");
    program.add_argument("--border-margin").default_value(10).scan<'i', int>().help("How far away a seed point must be from the image border.");

    // spatial information
    auto &centering = program.add_mutually_exclusive_group();
    centering.add_argument("--corner").default_value(false).implicit_value(true).help("If the bounding box/origin is at the corner of the voxel instead of the center.");
    centering.add_argument("--center").default_value(false).implicit_value(true).help("If the bounding box/origin is at the center of the voxel instead of the corner.");
    program.add_argument("--bbox").nargs(2,6).scan<'g', float>().help("The bounding box of the image");
    program.add_argument("--voxelsize").nargs(1,3).scan<'g', float>().help("The uniform voxel size");
    program.add_argument("--origin").nargs(1,3).scan<'g', float>().help("position of the image");

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
    std::filesystem::path input_name = program.get("input");
    // reporter.start("CRSF");

    // Arguments
    auto args = parsing::merge_spatial_arguments(program);
    if(!args)
    {
       return 1;
    }
 
    // load input
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Load " << input_name << std::endl;
    RawField<float> img;
    try
    {
        // TODO: check if nrrd values were float/double. If it could be labels, we should warn the user!
       img = RawField<float>::load(input_name, args.value());
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
 
    // Get min and max of image
    const auto img_min_max = img.min_max();
    const auto img_min = img_min_max.first;
    const auto img_max = img_min_max.second;
    const auto img_diff = img_max - img_min;
 
     // Settings:
    auto debug_setting = parsing::flatten_arg<std::string>(program, "--debug", "debug");
    bool debug = debug_setting.has_value();
    std::filesystem::path debug_path = debug_setting.value_or("debug");
     
    const auto fm_range = program.get<float>("--range");
    const auto min = program.present<float>("--min").value_or(img_min);
    const auto max = program.present<float>("--max").value_or(img_max + img_diff * 1.5f);
     
    auto seed_range = program.get<float>("--seed-range");
    auto default_seed_threshold = img_min;
    const auto seed_algo = program.get("algorithm");
    if(seed_algo == "fm"){
       default_seed_threshold = img_min + img.getVoxelSize().length() * 0.1; 
    }
    auto seed_threshold = program.present<float>("--seed-threshold").value_or(default_seed_threshold);
     
    const float default_shift_distance = fm_range * 0.1 * img.getVoxelSize().length();
    const float shift_distance_voxels = program.present<float>("--seed-shift").value_or(fm_range * 0.1);
    // const std::optional<float> shift_distance_voxels = parsing::flatten_arg<float>(program, "--seed-shift", fm_range * 0.1);
    const float shift_distance = shift_distance_voxels * img.getVoxelSize().length();
    
    int border_padding = program.get<int>("--border-padding");
    int border_margin = program.get<int>("--border-margin");
 
    // Print Settings
    std::cout << "Dims: " << img.dims() << std::endl;
    std::cout << "Voxelsize: " << img.getVoxelSize() << std::endl;
    if(debug){
        std::cout << "Debug activated with files saved at " << debug_path << std::endl;
    }
    std::cout << "FM range: " << fm_range << std::endl;
    std::cout << "Cost intensity range for marching set to [" << min << ", " << max << "]" << std::endl;

    std::cout << "seed algorithm of " << seed_algo << std::endl;
    std::cout << "seed generation with voxel range of " << (seed_range * fm_range) << std::endl;
    std::cout << "Seed threshold: " << seed_threshold << std::endl;
    if(shift_distance_voxels > 0){
        std::cout << "seed shifting distance of up to " << shift_distance_voxels << " voxels" << std::endl;
    }
    if(border_padding > 0){
        std::cout << "border padding of " << border_padding << std::endl;
    }
    if(border_margin > 0){
        std::cout << "border margin of " << border_margin << std::endl;
    }
 
     // generate seeds
     std::vector<ridgesurface::Seed> seeds;
        std::cout << "Generate seeds" << std::endl;
        std::vector<uint32_t> generated_seed_indices;

        if(seed_algo == "shell"){
            generated_seed_indices = greedyShellSampling(reporter, img.get_view(), seed_range * fm_range, min, max, seed_threshold, border_margin);
        }else{
            generated_seed_indices = greedyFmSampling(reporter, img.get_view(), seed_range * fm_range, min, max, seed_threshold, border_margin);
        }

        for (std::size_t i = 0; i < generated_seed_indices.size(); ++i) {
            VecSize gridLoc = Lattice::gridLocationFromCIndex(generated_seed_indices[i], img.dims());
            VecFloat worldPos = img.lattice().worldPosition(gridLoc);
            seeds.push_back(ridgesurface::Seed(ridgesurface::SeedPoint(VecFloat(worldPos.x(), worldPos.y(), worldPos.z())), fm_range));   
        }
        generated_seed_indices.clear();
        std::cout << "Generated " << seeds.size() << " seeds" << std::endl;

        if(debug){
            saveSeeds(debug_path / input_name.stem() / input_name.stem().concat("_rsf_generatedSeeds.obj"), seeds);
        }
 
     // extend image
     Dims originalDims = img.lattice().dims();
     if(border_padding > 0){
         std::cout << "Extend image by " << border_padding << " voxels" << std::endl;
         img = img.extend(border_padding);
         // fade out image
         for(auto pair : GridPositionBorderIterator(img.dims(), border_padding))
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
 
         ridgesurface::CubicalRidgeSurfaceFinder m_finder(reporter);
         m_finder.setInput(&img.get_view());
         m_finder.setThresholds(min,max);
         if(border_padding){
             VecSize min_clip = VecSize(border_padding);
             VecSize max_clip = static_cast<VecSize>(originalDims) + min_clip;
             m_finder.getTransformer().setClipping(min_clip, max_clip);
         }
 
         // shift seeds
         if(shift_distance > 0.0f){
             std::cout << "Shift seeds by a maximum of " << shift_distance  << " distance" << std::endl;
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
                     m_finder.moveSeedToRidge(seed, shift_distance);
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
                     m_finder.moveSeedToRidge(seed, shift_distance);
                 }
                 reporter.end();
             }
         }
 
        // if(border_margin > 0.0f)
        // {
        //      m_finder.settings.image_border_threshold = border_margin;
        //      auto skippedSeeds = std::vector<ridgesurface::Seed>();
        //      for(auto seed: seeds){
        //          if(m_finder.validSeed(seed)){
        //              m_finder.addSeed(seed);
        //          }else{
        //              skippedSeeds.push_back(seed);
        //          }
        //      }
        //      if(skippedSeeds.size() > 0){
        //          std::cout << "Skipped " << skippedSeeds.size() << " seeds" << std::endl;
        //          if(debug)
        //          {
        //              saveSeeds(debug_path / input_name.stem() / input_name.stem().concat("_rsf_skippedSeeds.obj"), skippedSeeds);
        //          }
        //      }
        // }else{
             for(auto seed: seeds){
                 m_finder.addSeed(seed);
             }
        // }
        
        if(debug){
            m_finder.recalculate().save_each_patch(debug_path / input_name.stem() / "patches" / input_name.stem().concat("_rsf_"));
        }else{
            m_finder.recalculate();
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
     std::cout << "MidMembrane finished! [" << (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count() / 1000000.0) << " sec]" << std::endl;
     return 0;
 }