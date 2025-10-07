/** Cubical Ridge Surface Finder (c) 2025 by Zuse Institute Berlin
 * 
 * Cubical Ridge Surface Finder is licensed under a
 * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * 
 * You should have received a copy of the license along with this
 * work.  If not, see <http://creativecommons.org/licenses/by-nc-sa/3.0/>.
 */

#include <iostream>
#include "Config.h"
#include "src/utils/ProgressbarReportHelper.hpp"
#include "src/io/npy.hpp"
#include "src/io/argparse.hpp"
#include "src/field/RawField.hpp"

#include "src/rsf/CubicalRidgeSurfaceFinder.h"
#include "src/io/vertexset.hpp"

#include <filesystem>
#include <vector>

// #include <chrono>
// #include <thread>

int main(int argc, char *argv[])
{
    // lets do some progressbar testing
    // auto progress = progressbar::ProgressbarReportDynamic(progressbar::FullReport);
    // progress.start("first task");
    // std::this_thread::sleep_for(std::chrono::seconds(3));
    // progress.update(0.5f);
    // std::this_thread::sleep_for(std::chrono::seconds(3));
    // progress.start("second task");
    // std::this_thread::sleep_for(std::chrono::seconds(3));
    // progress.update(0.3f);
    // std::this_thread::sleep_for(std::chrono::seconds(3));
    // progress.end();
    // std::this_thread::sleep_for(std::chrono::seconds(3));
    // progress.end();


    argparse::ArgumentParser program("Ridge Surface Finder", RidgeSurfaceFinder_VERSION_MAJOR + "." + RidgeSurfaceFinder_VERSION_MINOR);
    program.add_argument("input").help("Numpy input file to load");
    program.add_argument("seed").help("Numpy label for voxel seed or wavefront object");
    program.add_argument("--min").default_value(0.f).scan<'g', float>().help("Min threshold. Everything under is not touched.");
    program.add_argument("--max").default_value(1.2f).scan<'g', float>().help("Max threshold. Everything above or equal has no cost.");
    program.add_argument("--range").default_value(50.0f).scan<'g', float>().help("How far each seedpoint should be marched for");
    program.add_argument("--corner").default_value(false).implicit_value(true).help("If the bounding box/origin is at the corner of the voxel instead of the center.");
    program.add_argument("--bbox").nargs(2,6).scan<'g', float>().help("The bounding box of the image");
    program.add_argument("--voxelsize").nargs(1,3).scan<'g', float>().default_value(std::vector<float>{1.f}).help("The uniform voxel size");
    program.add_argument("--origin").nargs(1,3).scan<'g', float>().default_value(std::vector<float>{0.f}).help("position of the image");
    program.add_argument("-o", "--output").help("Name of output file");

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

    // load seeds (either labels or wavefront-like vertexset)
    std::vector<ridgesurface::Seed> seeds;
    try {
        // vertexset
        auto vertexset = read_vertexset(program.get("seed"), program.get<float>("range"));
        for(std::size_t i = 0; i < vertexset.points.size(); ++i){
            auto point = vertexset.points[i];
            seeds.push_back(ridgesurface::Seed(ridgesurface::SeedPoint(VecFloat(point.x(), point.y(), point.z())), vertexset.data[i]));
        }
    }
    catch (const std::exception &err) {
        // try labels
        try {
            auto labels = RawField<uint8_t>::load(program.get("seed"));
            // generate seeds from voxels
            for(std::size_t i = 0; i < labels.lattice().size(); ++i){
                if(labels.get(RawField<uint8_t>::FieldLocation(i)) == 0){
                    seeds.push_back(ridgesurface::Seed(ridgesurface::SeedPoint(static_cast<VecFloat>(Lattice::gridLocationFromCIndex(i, labels.dims()))), program.get<float>("--range")));
                }
            }
        }
        catch (const std::exception &err) {
            std::cerr << "Neither wavefront object nor numpy labels found for seeds" << std::endl;
            return 1;
        }
    }

    // run
    try
    {
        auto img = RawField<float>::load(program.get("input"));
        if(program.is_used("--bbox") && program.is_used("--voxelsize") 
            || program.is_used("--bbox") && program.is_used("--origin"))
        {
            //ERROR: Either bbox or voxelsize (and origin)
            std::cerr << "When a bounding box is given, no voxelsize or origin is allowed!" << std::endl;
            return 1;
        }

        if(program.is_used("--corner") && !program.is_used("--bbox") 
            || program.is_used("--corner") && !program.is_used("--origin"))
        {
            //INFO
            std::cout << "--corner argument ignored, as no origin or bounding box was given" << std::endl;
        }

        // load bbox or origin/voxelsize
        if(program.is_used("--bbox")){
            std::vector<float> vals = program.get<std::vector<float>>("--bbox");
            CornerBBox bbox = CornerBBox::fromCornerBbox(VecFloat(0.f), VecFloat(1.f));
            if(vals.size() == 2){
                if(program.get<bool>("--corner")){
                    bbox = CornerBBox::fromCornerBbox(VecFloat(vals[0]), VecFloat(vals[1]));
                }else{
                    bbox = CornerBBox::fromCenterBBox(VecFloat(vals[0]), VecFloat(vals[1]), img.dims());
                }
            }else if(vals.size() == 6){
                if(program.get<bool>("--corner")){
                    bbox = CornerBBox::fromCornerBbox(VecFloat(vals[0], vals[1], vals[2]), VecFloat(vals[3],vals[4],vals[5]));
                }else{
                    bbox = CornerBBox::fromCenterBBox(VecFloat(vals[0], vals[1], vals[2]), VecFloat(vals[3],vals[4],vals[5]), img.dims());
                }
            }else{
                //ERROR: we need either 2 or 6 float values
                std::cerr << "Either 2 or 6 float values are neccesary for a bounding box" << std::endl;
                return 1;
            }
            img.setBBox(bbox);
        }else{
            if(program.is_used("--origin")){
                std::vector<float> orig_vals = program.get<std::vector<float>>("--origin");
                VecFloat orig;
                if(orig_vals.size() == 1){
                    orig = VecFloat(orig_vals[0]);
                }else if(orig_vals.size() == 3){
                    orig = VecFloat(orig_vals[0], orig_vals[1], orig_vals[2]);
                }else{
                    //ERROR: we need either 1 or 3 float values
                    std::cerr << "Either 1 or 3 float values are neccesary for origin" << std::endl;
                    return 1;
                }
                img.setOrigin(orig);
            }
            if(program.is_used("--voxelsize")){
                std::vector<float> voxel_vals = program.get<std::vector<float>>("--voxelsize");
                VecFloat voxelsize;
                if(voxel_vals.size() == 1){
                    voxelsize = VecFloat(voxel_vals[0]);
                }else if(voxel_vals.size() == 3){
                    voxelsize = VecFloat(voxel_vals[0], voxel_vals[1], voxel_vals[2]);
                }else{
                    //ERROR: we need either 1 or 3 float values
                    std::cerr << "Either 1 or 3 float values are neccesary for voxelsize" << std::endl;
                    return 1;
                }
                img.setVoxelSize(voxelsize);
            }
        }

        // warn if voxelsize are wideley different
        const auto voxelsize = img.getVoxelSize();
        const float max_voxelsize = std::max(std::max(voxelsize.x(), voxelsize.y()), voxelsize.z());
        const float min_voxelsize = std::min(std::min(voxelsize.x(), voxelsize.y()), voxelsize.z());
        if((max_voxelsize - min_voxelsize) / min_voxelsize > 0.5){
            std::cout << "Warning: Voxelsizes are assumed to be equal in all axis for the algorithms. Results might not be satisfactory." << std::endl;
        }

        const auto min = program.get<float>("--min");
        const auto max = program.get<float>("--max");

        std::cout << "Loaded " << seeds.size() << " seeds" << std::endl;
        std::cout << "Cost intensity range for marching set to [" << min << ", " << max << "]" << std::endl;

        auto progressbar = progressbar::ProgressbarReportDynamic(progressbar::ProgressbarReportLevel::FullReport);
        ridgesurface::CubicalRidgeSurfaceFinder m_finder(progressbar);
        m_finder.setInput(&img);
        m_finder.newSeeds(seeds);
        m_finder.setThresholds(min,max);
        m_finder.recalculate();
        
        auto surf = surface::Surface();
        m_finder.finalize(&surf);
        const auto default_output_name = program.get("input").append("_rsf.obj");
        const auto output_name = program.present("-o").value_or(default_output_name);
        surf.save(output_name);
        std::cout << "Saved result in " << output_name << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
}