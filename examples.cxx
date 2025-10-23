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

#include <filesystem>
#include <vector>

#ifndef __DATAPATH_RAW__
#define __DATAPATH_RAW__ "placeholder"
#endif
#define __DATAPATH__ std::string(__DATAPATH_RAW__)

float monkey_saddle(float x, float y, float z){
    return x*x*x-3.0*x*y*y-z;
}

float smelt_petal(float x, float y, float z){
    return x+y+z+x*y*z;
}

enum ExampleID {
    MONKEY_SADDLE,
    SMELT_PETAL,
    SIZE
};

const char* ExamplesString[] = {
    "monkey-saddle",
    "smelt-petal"
};

int main(int argc, char *argv[])
{
    argparse::ArgumentParser program("Ridge Surface Finder Examples", RidgeSurfaceFinder_VERSION_MAJOR + "." + RidgeSurfaceFinder_VERSION_MINOR);
    program.add_argument("example").help("Example to generate.").choices(ExamplesString[0], ExamplesString[1]);
    // program.add_argument("--min").default_value(0.f).scan<'g', float>().help("Min threshold. Everything under is not touched.");
    // program.add_argument("--max").default_value(1.2f).scan<'g', float>().help("Max threshold. Everything above or equal has no cost.");
    // program.add_argument("--range").default_value(50.0f).scan<'g', float>().help("How far each seedpoint should be marched for");
    //  program.add_argument("--corner").default_value(false).implicit_value(true).help("If the bounding box/origin is at the corner of the voxel instead of the center.");
    //  program.add_argument("--bbox").nargs(2,6).scan<'g', float>().help("The bounding box of the image");
    //  program.add_argument("--voxelsize").nargs(1,3).scan<'g', float>().default_value(std::vector<float>{1.f}).help("The uniform voxel size");
    //  program.add_argument("--origin").nargs(1,3).scan<'g', float>().default_value(std::vector<float>{0.f}).help("position of the image");
    //  program.add_argument("-o", "--output").help("Name of output file");
 
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

    auto example_str = program.get("example");
    auto example = ExampleID::SIZE;
    for(auto i = 0; i < ExampleID::SIZE; ++i){
        if(example_str == ExamplesString[i]){
            example = static_cast<ExampleID>(i);
            break;
        }
    }

    // generate scalarfield
    auto exponent = 8.0f;

    auto dims = Dims(200);
    auto lattice = Lattice(dims,CornerBBox(VecFloat(-1.0f), VecFloat(1.0f)));
    auto img = RawField<float>(lattice);
    for(std::size_t z = 0; z < dims.z(); ++z){
        for(std::size_t y = 0; y < dims.y(); ++y){
            for(std::size_t x = 0; x < dims.x(); ++x){
                auto loc = VecSize(x,y,z);
                auto pos = lattice.worldPosition(loc);
                auto val = 0.0f;
                switch(example){
                    case ExampleID::MONKEY_SADDLE:
                        val = monkey_saddle(pos.x(), pos.y(), pos.z());
                        break;
                    case ExampleID::SMELT_PETAL:
                        val = smelt_petal(pos.x(), pos.y(), pos.z());
                        break;
                    default:
                        break;
                }
                val = powf(1.0f - std::min(abs(val), 1.0f), exponent);
                img.set(loc, val);
            }
        }
    }

    // save scalarfield
    img.save(__DATAPATH__ + fmt::format("/examples/{}_field", example_str));

    // load seeds
    std::vector<ridgesurface::Seed> seeds;
    try {
        // vertexset
        auto vertexset = read_vertexset(__DATAPATH__ + fmt::format("/examples/{}_seeds.obj", example_str), 50.0f);
        for(std::size_t i = 0; i < vertexset.points.size(); ++i){
            auto point = vertexset.points[i];
            seeds.push_back(ridgesurface::Seed(ridgesurface::SeedPoint(VecFloat(point.x(), point.y(), point.z())), vertexset.data[i]));
        }
    }
    catch (const std::exception &err) {
        std::cerr << fmt::format("No seeds found at {}/examples/{}_seeds.obj", __DATAPATH__, example_str) << std::endl;
        return 1;
    }
 
    // run
    try
    {
        const auto min = 0.0f;
        const auto max = 1.5f;

        auto progressbar = progressbar::ProgressbarReportDynamic(progressbar::ProgressbarReportLevel::FullReport);
        ridgesurface::CubicalRidgeSurfaceFinder m_finder(progressbar);
        m_finder.setInput(&img);
        m_finder.newSeeds(seeds);
        m_finder.setThresholds(min,max);
        m_finder.recalculate();
        
        auto surf = surface::Surface();
        m_finder.finalize(&surf);
        surf.save(__DATAPATH__ + fmt::format("/examples/{}_surf.obj", example_str));
        std::cout << fmt::format("Saved result in {}/examples/{}_surf.obj", __DATAPATH__, example_str) << std::endl;
     }
     catch (const std::exception &e)
     {
        std::cerr << e.what() << std::endl;
        return 1;
     }
 
    return 0;
}