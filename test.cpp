/** Cubical Ridge Surface Finder (c) 2025 by Zuse Institute Berlin
 * 
 * Cubical Ridge Surface Finder is licensed under a
 * Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
 * 
 * You should have received a copy of the license along with this
 * work.  If not, see <http://creativecommons.org/licenses/by-nc-sa/3.0/>.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>

// #include "src/field/RawField.hpp"
#include <fastmarching/FastMarchingUtilities.hpp>
#include <fastmarching/FastMarching.hpp>

#include <surface/Surface.hpp>

#include <rsf/Faces.hpp>
#include <rsf/FacesToCubicalMesh.hpp>

#include <utils/Permutation.hpp>

#include <rsf/CubicalRidgeSurfaceFinder.h>
#include <io/vertexset.hpp>

#include <iostream>

using Catch::Matchers::WithinULP;
using Catch::Matchers::UnorderedEquals;

TEST_CASE("Load Files", "[npy][txt][io]")
{
    // check if type is correct
    // auto stype = npy::probe_npy_scalartype("test_data/volumes/grayimage_234");
    // std::cout << stype << std::endl;

    // load image
    auto img = RawField<float>::load("test_data/volumes/grayimage_234");
    REQUIRE(img.dims() == Dims(4, 3, 2));

    // check if fortran order works correctly
    auto fortran = RawField<float>::load("test_data/volumes/grayimage_234_fortran");
    REQUIRE(fortran.dims() == Dims(4, 3, 2));

    std::size_t counter = 0;
    for (std::size_t z = 0; z < img.dims().z(); ++z)
    {
        for (std::size_t y = 0; y < img.dims().y(); ++y)
        {
            for (std::size_t x = 0; x < img.dims().x(); ++x)
            {
                REQUIRE(img.get(VecSize(x, y, z)).value() == fortran.get(VecSize(x, y, z)).value());
                ++counter;
            }
        }
    }

    // metadata
    REQUIRE(img.getCornerBoundingBox().equals(CornerBBox(VecFloat(0), VecFloat(4,3,2))));
    REQUIRE(img.getVoxelSize() == VecFloat(1.0));
    REQUIRE(img.getCenterBoundingBox().min_center() == VecFloat(0.5));

    REQUIRE(fortran.getCenterBoundingBox().equals(CenterBBox(VecFloat(12.3), VecFloat(12.6,12.5,12.4))));
    REQUIRE(fortran.getVoxelSize() == VecFloat(0.1));
    REQUIRE(fortran.getCornerBoundingBox().min_corner() == VecFloat(12.25));
}

TEST_CASE("Fast Marching Utilities", "[fm]")
{
    // load image
    auto img = RawField<float>::load("test_data/volumes/grayimage01_333");
    // calculate fm gradient for location 1,1,1
    auto grad = fastmarching::gradient(img, img.createLocation(VecSize(1, 1, 1)));
    // TODO: this should maybe be negative instead of positive...
    REQUIRE(grad.equals(VecFloat(0.70147467, 0.13889703, 0.76361562)));

    auto singular1 = RawField<float>::load("test_data/volumes/singularx1");
    auto singular2 = RawField<float>::load("test_data/volumes/singularx2");
    auto singular3 = RawField<float>::load("test_data/volumes/singularx3");

    std::size_t counter = 0;
    for (std::size_t z = 0; z < singular1.dims().z(); ++z)
    {
        for (std::size_t y = 0; y < singular1.dims().y(); ++y)
        {
            for (std::size_t x = 0; x < singular1.dims().x(); ++x)
            {
                REQUIRE(std::min(fastmarching::euclidean_field(singular1, singular1.createLocation(VecInt(x, y, z))), singular1.get(VecSize(x, y, z)).value()) == singular2.get(VecSize(x, y, z)).value());
                REQUIRE(std::min(fastmarching::euclidean_field(singular2, singular2.createLocation(VecInt(x, y, z))), singular2.get(VecSize(x, y, z)).value()) == singular3.get(VecSize(x, y, z)).value());
                ++counter;
            }
        }
    }

    auto pot = RawField<float>::load("test_data/volumes/grayimage01_333");
    auto vals = RawField<float>::load("test_data/volumes/grayimage02_333");
    auto eikonal = RawField<float>::load("test_data/volumes/eikonal_01x02");
    auto euclidean = RawField<float>::load("test_data/volumes/euclidean02");

    // check if eikonal and euclidean are the same as in amira
    counter = 0;
    for (std::size_t z = 0; z < pot.dims().z(); ++z)
    {
        for (std::size_t y = 0; y < pot.dims().y(); ++y)
        {
            for (std::size_t x = 0; x < pot.dims().x(); ++x)
            {
                REQUIRE_THAT(fastmarching::euclidean_field(vals, vals.createLocation(VecInt(x, y, z))), WithinULP(euclidean.get(VecSize(x, y, z)).value(), 10));
                REQUIRE_THAT(fastmarching::eikonal_field(pot, 0.0, 1.0, vals, vals.createLocation(VecInt(x, y, z))), WithinULP(eikonal.get(VecSize(x, y, z)).value(), 10));
                ++counter;
            }
        }
    }
}

TEST_CASE("Fast Marching", "[fm]")
{
    progressbar::ProgressbarReportNone progress;

    SECTION("Probobilistic field"){
        // load images
        auto pot = RawField<float>::load("test_data/volumes/10x10x10");
        auto time = RawField<float>::load("test_data/volumes/10x10x10_time");
        auto dist = RawField<float>::load("test_data/volumes/10x10x10_dist");
        // create fm with two RawFields, and use ArrayMappingView (pointer) for now...
        auto time_view = mutil::ArrayMappingView<uint64_t,float>(time.data(), time.getDims().size());
        auto dist_view = mutil::ArrayMappingView<uint64_t,float>(dist.data(), dist.getDims().size());
        fastmarching::FastMarching<fastmarching::ObserverAll<mutil::ArrayMappingView<uint64_t, float>, BigHashSet<uint64_t>>, progressbar::ProgressbarReportNone> fm(progress);
        // start fm from the top left (0,0,0)
        fm.data(&pot, time_view, dist_view);
        fm.thresholds(0.0f,1.0f);
        fm.setStartVoxel(0);
        fm.march();
        // compare with amira result
        auto counter = 0;
        for (std::size_t z = 0; z < pot.dims().z(); ++z)
        {
            for (std::size_t y = 0; y < pot.dims().y(); ++y)
            {
                for (std::size_t x = 0; x < pot.dims().x(); ++x)
                {
                    REQUIRE_THAT(fm.time().get(counter), WithinULP(time.get(VecSize(x, y, z)).value(), 10));
                    REQUIRE_THAT(fm.distance().get(counter), WithinULP(dist.get(VecSize(x, y, z)).value(), 10));
                    ++counter;
                }
            }
        }
    }

    SECTION("Simple Ridge Field"){
        auto pot = RawField<float>::load("test_data/surfaces/simple_ridge");
        // using MappingView = mutil::HashMapMappingView<SmallHashMap<int64_t, float>>;
        // fastmarching::FastMarching<fastmarching::ObserverDistance<MappingView>> m_fm;
        fastmarching::FastMarching<fastmarching::ObserverDistance<mutil::HashMapMappingView<SmallHashMap<int64_t, float>>>, progressbar::ProgressbarReportNone> fm(progress);
        fm.data(&pot);
        fm.thresholds(0.5f,1.5f);
        fm.setStartPoint(VecFloat(5.f));
        fm.march();

        // debug: save
        auto dist = RawField<float>(pot.dims());
        fm.distance(dist.data(), dist.dims().size(), false);
        dist.save("test_data/output/simple_ridge_dist_test");
        
        // compare with amira result
        auto time_ref = RawField<float>::load("test_data/surfaces/simple_ridge_time_p05-p15_nocap");
        auto dist_ref = RawField<float>::load("test_data/surfaces/simple_ridge_dist_p05-p15_nocap");
        auto counter = 0;
        for (std::size_t z = 0; z < pot.dims().z(); ++z)
        {
            for (std::size_t y = 0; y < pot.dims().y(); ++y)
            {
                for (std::size_t x = 0; x < pot.dims().x(); ++x)
                {
                    REQUIRE_THAT(fm.time().get(counter), WithinULP(time_ref.get(VecSize(x, y, z)).value(), 10));
                    REQUIRE_THAT(fm.distance().get(counter), WithinULP(dist_ref.get(VecSize(x, y, z)).value(), 10));
                    ++counter;
                }
            }
        }
    }
}

TEST_CASE("Wavefront Files", "[obj][io]")
{
    SECTION("default singular global patch"){
        auto data = read_wavefront("test_data/surfaces/unit_cube");
        std::vector<VecFloat> points = {
            VecFloat(0.f,0.f,0.f),
            VecFloat(0.f,0.f,1.f),
            VecFloat(0.f,1.f,0.f),
            VecFloat(0.f,1.f,1.f),
            VecFloat(1.f,0.f,0.f),
            VecFloat(1.f,0.f,1.f),
            VecFloat(1.f,1.f,0.f),
            VecFloat(1.f,1.f,1.f),
        };
        REQUIRE_THAT(data.points, UnorderedEquals(points));
        REQUIRE(data.patches.size() == 0);

        auto cube = surface::Surface::load("test_data/surfaces/unit_cube");
        REQUIRE(cube.number_of_points() == 8);
        REQUIRE(cube.number_of_trianlges() == 12);
        REQUIRE(cube.number_of_patches() == 1);
    }

    SECTION("small surface example with explicit patches"){
        auto data = read_wavefront("test_data/surfaces/house");
        std::vector<VecFloat> points = {
            VecFloat(1.f,0.f,1.f),
            VecFloat(-1.f,0.f,1.f),
            VecFloat(-1.f,0.f,-1.f),
            VecFloat(1.f,0.f,-1.f),
            VecFloat(1.f,1.f,1.f),
            VecFloat(-1.f,1.f,1.f),
            VecFloat(1.f,1.f,-1.f),
            VecFloat(-1.f,1.f,-1.f),
            VecFloat(0.f,1.5f,1.f),
            VecFloat(0.f,1.5f,-1.f)
        };
        REQUIRE_THAT(data.points, UnorderedEquals(points));
        std::vector<std::array<std::size_t, 3>> triangles = {
            {0,1,3},
            {1,2,3},
            {4,5,7},
            {5,6,7},
            {0,4,1},
            {4,5,1},
            {1,5,2},
            {5,6,2},
            {2,6,3},
            {6,7,3},
            {3,7,0},
            {7,4,0},
            {4,8,5},
            {6,9,7},
            {4,7,8},
            {7,9,8},
            {5,8,6},
            {6,8,9}
        };
        REQUIRE_THAT(data.triangles, UnorderedEquals(triangles));
        std::vector<std::size_t> patches = {
            2,
            4,
            12
        };
        REQUIRE_THAT(data.patches, UnorderedEquals(patches));
    
        auto surf = surface::Surface::load("test_data/surfaces/house");
        REQUIRE(surf.number_of_points() == 10);
        REQUIRE(surf.number_of_trianlges() == 18);
        REQUIRE(surf.number_of_patches() == 4);
    }
}

TEST_CASE("Permutation", "[perm]")
{
    std::vector<VecFloat> src = {
        VecFloat(0,1,0),
        VecFloat(0,0,0),
        VecFloat(0,0,1),
        VecFloat(0,1,1),
        VecFloat(1,1,1),
        VecFloat(1,0,1),
        VecFloat(1,0,0),
        VecFloat(1,1,0)
    };

    auto sperm = mutil::Permutation::sortPermutation(src,  
        [&](const VecFloat& a, const VecFloat& b) {
            return a.lexicographic_order_epsilon(b);
        }
    );

    std::vector<std::size_t> sorted = {
        1, 2, 0, 3, 6, 5, 7, 4
    };

    for(std::size_t i = 0; i < sperm.size(); ++i){
        REQUIRE(sperm[i] == sorted[i]);
    }

    std::vector<VecFloat> target = {
        VecFloat(0,0,0),
        VecFloat(1,0,0),
        VecFloat(0,0,1),
        VecFloat(1,0,1),
        VecFloat(0,1,0),
        VecFloat(1,1,0),
        VecFloat(0,1,1),
        VecFloat(1,1,1)
    };

    auto mperm = mutil::Permutation::mapPermutation(src, target, 
        [&](const VecFloat& a, const VecFloat& b) {
            return a.lexicographic_order_epsilon(b);
        }
    );
    
    for(std::size_t i = 0; i < mperm.size(); ++i){
        REQUIRE(src[i] == target[mperm[i]]);
    }
}


TEST_CASE("Face to Surface", "[face]")
{
    // generate a unit cube in face form (same start voxel, use all direction)
    Face faces[6] = {
        Face(0, Direction::LEFT),
        Face(0, Direction::RIGHT),
        Face(0, Direction::UP),
        Face(0, Direction::DOWN),
        Face(0, Direction::FORWARD),
        Face(0, Direction::BACKWARD)
    };
    auto builder = FacesToCubicalMesh();
    auto surf = surface::Surface();
    // origin is the center of the first voxel, which is 0.5f
    builder.populateSurface(&surf, Lattice(Dims(1,1,1), VecFloat(0.5f), VecFloat(1.0f)), &faces[0], &faces[6]);
    
    auto ref = surface::Surface::load("test_data/surfaces/unit_cube");
    REQUIRE(surf.number_of_points() == 8);
    REQUIRE(surf.number_of_trianlges() == 12);
    REQUIRE(surf.number_of_patches() == 1);
    REQUIRE(surf.equals(ref));
}

// TEST_CASE("Change Graph", "[graph]")
// {
//     //TODO: Test change graph
// }

// TEST_CASE("Flooding", "[flood]")
// {
//     //TODO test flooding
// }

TEST_CASE("Ridge Surface Finder", "[rsf]")
{
    auto progressbar = progressbar::ProgressbarReportDynamic();
    ridgesurface::CubicalRidgeSurfaceFinder m_finder(progressbar);

    SECTION("plane surface"){
        progressbar.level(progressbar::NoReport);
        auto img = RawField<float>::load("test_data/surfaces/simple_ridge");
        REQUIRE(img.dims() == Dims(10));
        REQUIRE(img.getCornerBoundingBox() == CornerBBox::fromCenterBBox(VecFloat(0.f), VecFloat(9.f), img.dims()));

        m_finder.setInput(&img);
        m_finder.addSeed(ridgesurface::Seed::Seedpoint(VecFloat(5.f), 50.0f));
        m_finder.setThresholds(0.5f,1.5f);
        m_finder.compute();
        
        auto labels = RawField<uint16_t>(img.dims());
        m_finder.castToLabels(labels.data());
        labels.save("test_data/output/simple_ridge_labels");

        // m_finder.patchedSurface().save("test_data/output/simple_ridge_test");
        const auto ref = surface::Surface::load("test_data/surfaces/simple_ridge_p05-p15");
        auto res = surface::Surface();
        m_finder.finalize(&res);
        REQUIRE(res.equals(ref));
        REQUIRE(res.number_of_patches() == 1);
        REQUIRE(m_finder.patchedSurface().equals(ref));
        REQUIRE(m_finder.patchedSurface().number_of_patches() == 1);
    }
}

// Hidden by default, as test needs ~30 seconds
TEST_CASE("Ridge Surface Finder", "[.][rsf][exhaustive]")
{
    auto progressbar = progressbar::ProgressbarReportDynamic();
    ridgesurface::CubicalRidgeSurfaceFinder m_finder(progressbar);

    // TODO: check if file exists, otherwise skip
    RawField<float> img;
    try {
        img = RawField<float>::load("test_data/monkey_saddle/monkey-capped-4");
    }
    catch (const std::ios_base::failure& e){
        SKIP("monkey_saddle file does not exist and test is skipped");
    }

    SECTION("monkey saddle"){
        progressbar.level(progressbar::FullReport);
        auto seedpoints = read_vertexset("test_data/monkey_saddle/seedpoints", 50.0f);

        REQUIRE(img.getCenterBoundingBox() == CenterBBox(VecFloat(-1.0f), VecFloat(1.0f)));
        REQUIRE(img.getVoxelSize() == VecFloat(0.005f));

        m_finder.setInput(&img);
        m_finder.setThresholds(0.f,1.5f);
        // only add first seed yet
        m_finder.addSeed(ridgesurface::Seed::Seedpoint(seedpoints.points[0], seedpoints.data[0]));
        m_finder.compute();

        // add all other seeds
        for(std::size_t i = 1; i < seedpoints.points.size(); ++i){
            m_finder.addSeed(ridgesurface::Seed::Seedpoint(seedpoints.points[i], seedpoints.data[i]));
        }
        m_finder.compute();

        auto labels = RawField<uint16_t>(img.lattice());
        m_finder.castToLabels(labels.data());
        const auto ref_labels = RawField<uint16_t>::load("test_data/monkey_saddle/label-00-p15");

        auto res = surface::Surface();
        m_finder.finalize(&res);
        const auto ref_res = surface::Surface::load("test_data/monkey_saddle/ridge-00-p15");

        const auto ref_patched = surface::Surface::load("test_data/monkey_saddle/patched-00-p15");

        REQUIRE(m_finder.patchedSurface().equals(ref_patched));
        REQUIRE(res.equals(ref_res));

        for (std::size_t z = 0; z < img.dims().z(); ++z)
        {
            for (std::size_t y = 0; y < img.dims().y(); ++y)
            {
                for (std::size_t x = 0; x < img.dims().x(); ++x)
                {
                    auto coord = VecSize(x, y, z);
                    REQUIRE(labels.get(coord) == ref_labels.get(coord));
                }
            }
        }

        // labels.save("test_data/output/monkey_saddle_test_label");
        // m_finder.patchedSurface().save("test_data/output/monkey_saddle_test_patched");
        // res.save("test_data/output/monkey_saddle_test");
    }
}


TEST_CASE("BitFace", "[face]"){
    const auto face1 = BitFace(VecSize(1,2,3), Direction::UP);
    REQUIRE(face1.toIndex() == 0x2000030000200001);
    REQUIRE(face1.startVoxel() == VecSize(1,2,3));
    const auto face2 = BitFace(VecSize(1,2,3), Direction::DOWN);
    REQUIRE(face2.toIndex() == 0xa000030000200001);
    REQUIRE(face2.startVoxel() == VecSize(1,2,3));
    const auto face3 = BitFace(VecSize(1,2,3), Direction::LEFT);
    REQUIRE(face3.toIndex() == 0x9000030000200001);
    REQUIRE(face3.startVoxel() == VecSize(1,2,3));



}