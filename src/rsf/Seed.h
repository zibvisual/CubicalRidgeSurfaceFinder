#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <queue>
#include <optional>
#include <cassert>        // assert 

#include <utils/Vec.hpp>
#include <utils/Dims.hpp>
#include <utils/Mapping.hpp>
#include <utils/Lattice.hpp>
#include <utils/Range.hpp>

namespace ridgesurface
{

    enum SeedEnum
    {
        UNKNOWN,
        POINT,
        LINE
    };

    class SeedPoint
    {
    public:
        VecFloat point;

        SeedPoint(VecFloat point)
            : point(point)
        {
        }

        bool operator==(const SeedPoint& rhs) const
        {
            return (point == rhs.point);
        }
    };

    class SeedLine
    {
    public:
        std::vector<VecFloat> points;

        SeedLine(std::vector<VecFloat> points)
            : points(points)
        {}

        bool operator==(const SeedLine& rhs) const
        {
            return (points == rhs.points);
        }
    };

    class SeedPointIterator
    {
    public:
        SeedPointIterator(VecFloat point, Lattice lattice)
            : point(point)
            , lattice(lattice)
            , inner(true)
        {}

        std::optional<std::size_t>
        next()
        {
            if (inner)
            {
                std::size_t pointId = lattice.c_index(point);
                inner = false;
                return pointId;
            }
            return std::optional<std::size_t>();
        }

        using value_type = std::size_t;
        // ADD_RANGE_TO_RUST_ITERATOR(SeedPointIterator)

    private:
        VecFloat point;
        Lattice lattice;
        bool inner;
    };

    class SeedLineIterator
    {
    public:
        SeedLineIterator(std::vector<VecFloat>& points, Lattice lattice)
            : points(points)
            , index(0)
            , lattice(lattice)
        {}

        std::optional<std::size_t>
        next()
        {
            return index < points.size() ? lattice.c_index(points[index++]) : std::optional<std::size_t>();
        }

        using value_type = std::size_t;
        // ADD_RANGE_TO_RUST_ITERATOR(SeedLineIterator)

        // copy-assignment
        SeedLineIterator&
        operator=(const SeedLineIterator& sli)
        {
            // Guard self assignment
            // if (this == &sli)
            //     return *this;

            // we just copy the reference
            points = sli.points;
            lattice = sli.lattice;

            return *this;
        }

    private:
        std::vector<VecFloat>& points;
        std::size_t index;
        Lattice lattice;
    };

    class SeedIterator;

    class Seed
    {
        friend SeedIterator;
        friend std::hash<Seed>;

    protected:
        SeedEnum id;

        union
        {
            SeedPoint point;
            SeedLine line;
        };

        float max;

        void
        destroy_value()
        {
            switch (id)
            {
                case LINE:
                    (&line)->SeedLine::~SeedLine();
                    break;
                case POINT: // nothing to do
                default:
                    break;
            }
        }

        void
        copy_value(const Seed& seed)
        {
            switch (id)
            {
                case LINE:
                    line = seed.line;
                    break;
                case POINT:
                    point = seed.point;
                    break;
                default:
                    break;
            }
        }

    public:
        // Seed needs a default constructor for hashmap, but it should never be used otherwise
        Seed()
            : id(UNKNOWN)
            , max(0.0f)
        {}

        Seed(SeedPoint seedpoint, float distance)
            : id(POINT)
            , point({ seedpoint })
            , max(distance)
        {}

        Seed(SeedLine seedline, float distance)
            : id(LINE)
            , line({ seedline })
            , max(distance)
        {}

        Seed(const Seed& seed)
            : id(seed.id)
            , max(seed.max)
        {
            copy_value(seed);
        }

        ~Seed()
        {
            destroy_value();
        }

        // copy-assignment
        Seed&
        operator=(const Seed& seed)
        {
            // Guard self assignment
            if (this == &seed)
                return *this;
            destroy_value();
            id = seed.id;
            max = seed.max;
            copy_value(seed);
            return *this;
        }

        static Seed
        Seedpoint(VecFloat point, float distance)
        {
            return Seed(SeedPoint(point), distance);
        }

        static Seed
        Seedline(std::vector<VecFloat> points, float distance)
        {
            return Seed(SeedLine(points), distance);
        }

        SeedIterator getVoxelSources(Lattice lattice);

        float
        getDistance() const
        {
            return max;
        }

        void
        setDistance(float distance)
        {
            max = distance;
        }

        bool
        isPoint() const
        {
            return id == POINT;
        }

        bool
        isLine() const
        {
            return id == LINE;
        }

        // return a point representing the seed (for spatial graph, mostly used in seedpoints)
        VecFloat
        firstPoint() const
        {
            switch (id)
            {
                case LINE:
                    return line.points[0];
                case POINT:
                    return point.point;
                default: // not reachable
                    assert(false);
                    throw std::runtime_error("Not Reachable");
            }
        }

        bool
        operator==(const Seed& rhs) const
        {
            bool val = ((id == rhs.id) & (max == rhs.max));
            if (!val)
            {
                return false;
            }

            switch (id)
            {
                case LINE:
                    return line == rhs.line;
                case POINT:
                    return point == rhs.point;
                default: // not reachable
                    assert(false);
                    throw std::runtime_error("Not Reachable");
            }
        }
    };

    class SeedIterator
    {
    protected:
        SeedEnum id;

        union
        {
            SeedPointIterator point;
            SeedLineIterator line;
        };

        void
        destroy_value()
        {
            switch (id)
            {
                case POINT:
                    (&point)->SeedPointIterator::~SeedPointIterator();
                    break;
                case LINE:
                    (&line)->SeedLineIterator::~SeedLineIterator();
                    break;
                default:
                    break;
            }
        }

    public:
        SeedIterator(Seed seed, Lattice lattice)
            : id(seed.id)
        {
            switch (id)
            {
                case POINT:
                    point = SeedPointIterator(seed.point.point, lattice);
                    break;
                case LINE:
                    line = SeedLineIterator(seed.line.points, lattice);
                    break;
                default:
                    break;
            }
        }

        std::optional<std::size_t>
        next()
        {
            switch (id)
            {
                case POINT:
                    return point.next();
                case LINE:
                    return line.next();
                default:
                    return std::optional<std::size_t>();
            }
        }

        using value_type = std::size_t;

        ADD_RANGE_TO_RUST_ITERATOR(SeedIterator)
    };
}

template<>
struct std::hash<ridgesurface::SeedLine>
{
    std::size_t operator()(const ridgesurface::SeedLine& seed) const
    {
        // TODO
        throw std::runtime_error("Not implemented yet");
    }
};

template<>
struct std::hash<ridgesurface::SeedPoint>
{
    std::size_t operator()(const ridgesurface::SeedPoint& seed) const
    {
        return std::hash<VecFloat>()(seed.point);
    }
};

template<>
struct std::hash<ridgesurface::Seed>
{
    std::size_t operator()(const ridgesurface::Seed& seed) const
    {
        switch (seed.id)
        {
            case ridgesurface::LINE:
                return (std::hash<ridgesurface::SeedEnum>()(seed.id)
                ^ (std::hash<ridgesurface::SeedLine>()(seed.line) << 1));
            case ridgesurface::POINT:
                return (std::hash<ridgesurface::SeedEnum>()(seed.id)
                ^ (std::hash<ridgesurface::SeedPoint>()(seed.point) << 1));
            default: // not reachable
                assert(false);
                throw std::runtime_error("Not Reachable");
        }

    }
};