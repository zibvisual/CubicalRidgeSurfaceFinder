#pragma once

#include <cstdint>
#include <stdexcept>
#include <optional>

#include <utils/Dims.hpp>
#include <utils/Lattice.hpp>
#include <utils/Direction.hpp>
#include <field/RawField.hpp>

#include <boost/functional/hash.hpp>

/**
 * We save 20 bit for each of the 3 dimensions, as well as a sign bit and a bit for the x,y or z direction.
 * That is we have sign bit, z,y,x bit for direction an then 20 bits for start position z, then y, then x
 */
class BitFace 
{
private:
    uint64_t m_index;

    inline static uint64_t convert(VecSize voxel, Direction direction){
        return (static_cast<uint64_t>(!direction.sign()) << 63) | (static_cast<uint64_t>(1) << (60 + direction.abs())) | (voxel.z() << 40) | (voxel.y() << 20) | voxel.x();
    }

    inline static VecSize toPos(uint64_t val){
        return VecSize(
            val & 0xFFFFF,
            (val >> 20) & 0xFFFFF,
            (val >> 40) & 0xFFFFF
        );
    }

    inline static bool dirSign(uint64_t val){
        return (val >> 63) & 1;
    }

    inline static VecInt toDir(uint64_t val){
        bool sign = dirSign(val);
        auto vec = VecInt(
            (val >> 60) & 1,
            (val >> 61) & 1,
            (val >> 62) & 1
        );
        return sign ? -vec : vec;
    }

    inline static VecInt movePos(uint64_t val){
        return static_cast<VecInt>(toPos(val)) + toDir(val);
    }

    BitFace(uint64_t index) : m_index(index){}

public:
    inline bool operator==(const BitFace& other) const
    {
        return (m_index == other.m_index);
    }

    inline bool operator!=(const BitFace& other) const
    {
        return (m_index != other.m_index);
    }

    BitFace()
        : m_index(0)
    {}

    BitFace(VecSize voxel, Direction direction)
        : m_index(convert(voxel, direction))
    {}

    // Face(VecSize voxel, VecSize endVoxel)
    //     : voxel(startVoxel)
    //     , direction(Direction::LEFT)
    // {
    //     auto diff = endVoxel - startVoxel;
    //     if (diff == -1)
    //         direction = Direction::LEFT;
    //     else if (diff == 1)
    //         direction = Direction::RIGHT;
    //     else if (diff == -dims[0])
    //         direction = Direction::DOWN;
    //     else if (diff == dims[0])
    //         direction = Direction::UP;
    //     else if (diff == -dims[0] * dims[1])
    //         direction = Direction::BACKWARD;
    //     else if (diff == dims[0] * dims[1])
    //         direction = Direction::FORWARD;
    //     else
    //         throw std::invalid_argument("startVoxel and endVoxel are not neighbors to each other");
    // }

    inline static BitFace fromIndex(uint64_t index)
    {
        return BitFace(index);
    }

    inline uint64_t toIndex() const
    {
        return m_index;
    }

    inline VecSize startVoxel() const
    {
        return toPos(m_index);
    }

    inline VecInt endVoxel() const
    {
        return movePos(m_index);
    }

    // uint64_t startId() const
    // {
    //     return voxel;
    // }

    // uint64_t endId(Dims dims) const
    // {
    //     // TODO: add/subtract to the correct place
    //     if (direction == Direction::LEFT)
    //         return voxel - 1;
    //     else if (direction == Direction::RIGHT)
    //         return voxel + 1;
    //     else if (direction == Direction::DOWN)
    //         return voxel - dims[0];
    //     else if (direction == Direction::UP)
    //         return voxel + dims[0];
    //     else if (direction == Direction::BACKWARD)
    //         return voxel - dims[0] * dims[1];
    //     else //(direction == Direction::FORWARD)
    //         return voxel + dims[0] * dims[1];
    // }

    // // this is only true for RawField<float>
    // std::optional<uint64_t> validEndId(Dims dims) const
    // {
    //     auto vec = startVoxel(dims);
    //     if (direction == Direction::LEFT)
    //         return vec[0] == 0 ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel - 1);
    //     else if (direction == Direction::RIGHT)
    //         return vec[0] == dims[0] ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel + 1);
    //     else if (direction == Direction::DOWN)
    //         return vec[1] == 0 ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel - dims[0]);
    //     else if (direction == Direction::UP)
    //         return vec[1] == dims[1] ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel + dims[0]);
    //     else if (direction == Direction::BACKWARD)
    //         return vec[2] == 0 ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel - dims[0] * dims[1]);
    //     else //(direction == Direction::FORWARD)
    //         return vec[2] == dims[2] ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel + dims[0] * dims[1]);
    // }

};

class Face
{
public:
    uint64_t voxel;
    // direction is just empty_cell - voxel;
    Direction direction;

    bool operator==(const Face& other) const
    {
        return (voxel == other.voxel
                && direction == other.direction);
    }

    bool operator!=(const Face& other) const
    {
        return (voxel != other.voxel
                || direction != other.direction);
    }

    // Face() = default;
    Face()
        : voxel(0)
        , direction(Direction::LEFT)
    {}

    Face(uint64_t voxel, Direction direction)
        : voxel(voxel)
        , direction(direction)
    {}

    Face(uint64_t startVoxel, uint64_t endVoxel, Dims dims)
        : voxel(startVoxel)
        , direction(Direction::LEFT)
    {
        auto diff = endVoxel - startVoxel;
        if (diff == -1)
            direction = Direction::LEFT;
        else if (diff == 1)
            direction = Direction::RIGHT;
        else if (diff == -dims[0])
            direction = Direction::DOWN;
        else if (diff == dims[0])
            direction = Direction::UP;
        else if (diff == -dims[0] * dims[1])
            direction = Direction::BACKWARD;
        else if (diff == dims[0] * dims[1])
            direction = Direction::FORWARD;
        else
            throw std::invalid_argument("startVoxel and endVoxel are not neighbors to each other");
    }

    static Face fromIndex(uint64_t index)
    {
        return Face(index / 6, static_cast<Direction::Value>(index % 6));
    }

    uint64_t toIndex() const
    {
        return voxel * 6 + static_cast<uint64_t>(direction);
    }

    VecSize startVoxel(Dims dims) const
    {
        return RawField<float>::FieldLocation(voxel).location(dims);
    }

    VecInt endVoxel(Dims dims) const
    {
        const auto start = startVoxel(dims);
        if (direction == Direction::LEFT)
            return VecInt(start[0] - 1, start[1], start[2]);
        else if (direction == Direction::RIGHT)
            return VecInt(start[0] + 1, start[1], start[2]);
        else if (direction == Direction::DOWN)
            return VecInt(start[0], start[1] - 1, start[2]);
        else if (direction == Direction::UP)
            return VecInt(start[0], start[1] + 1, start[2]);
        else if (direction == Direction::BACKWARD)
            return VecInt(start[0], start[1], start[2] - 1);
        else //(direction == Direction::FORWARD)
            return VecInt(start[0], start[1], start[2] + 1);
    }

    uint64_t startId() const
    {
        return voxel;
    }

    // this is only true for RawField<float>
    uint64_t endId(Dims dims) const
    {
        if (direction == Direction::LEFT)
            return voxel - 1;
        else if (direction == Direction::RIGHT)
            return voxel + 1;
        else if (direction == Direction::DOWN)
            return voxel - dims[0];
        else if (direction == Direction::UP)
            return voxel + dims[0];
        else if (direction == Direction::BACKWARD)
            return voxel - dims[0] * dims[1];
        else //(direction == Direction::FORWARD)
            return voxel + dims[0] * dims[1];
    }

    // this is only true for RawField<float>
    std::optional<uint64_t> validEndId(Dims dims) const
    {
        auto vec = startVoxel(dims);
        if (direction == Direction::LEFT)
            return vec[0] == 0 ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel - 1);
        else if (direction == Direction::RIGHT)
            return vec[0] == dims[0] ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel + 1);
        else if (direction == Direction::DOWN)
            return vec[1] == 0 ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel - dims[0]);
        else if (direction == Direction::UP)
            return vec[1] == dims[1] ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel + dims[0]);
        else if (direction == Direction::BACKWARD)
            return vec[2] == 0 ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel - dims[0] * dims[1]);
        else //(direction == Direction::FORWARD)
            return vec[2] == dims[2] ? std::optional<uint64_t>() : std::optional<uint64_t>(voxel + dims[0] * dims[1]);
    }

    Face flipDirection(Dims dims) const {
        return Face {
            voxel + direction.c_offset(dims),
            direction.flip()
        };
    }

    inline std::array<Face, 12> possibleBorderNeighbors(Dims dims, Face none) const
    {
        auto neighbors = std::array<Face, 12>();

        const VecSize voxel = this->startVoxel(dims);
        const std::array<int64_t, 6> direction = {
            1,
            static_cast<int64_t>(dims[0]),
            static_cast<int64_t>(dims[0] * dims[1]),
            -1,
            -static_cast<int64_t>(dims[0]),
            -static_cast<int64_t>(dims[0] * dims[1]),
        };

        const std::array<uint64_t, 6> grid_check = {
            dims[0] - 1,
            dims[1] - 1,
            dims[2] - 1,
            0,
            0,
            0,
        };

        // We follow right hand rule
        std::array<int, 4> borders = { (this->direction + 4) % 6, (this->direction + 2) % 6, (this->direction + 1) % 6, (this->direction + 5) % 6 };

        std::size_t counter = 0;
        for (int border : borders)
        {
            int border_axis = border % 3;
            int face_axis = this->direction % 3;

            // first possible face is always possible -> no grid check
            Direction first_dir = static_cast<Direction::Value>(border);
            neighbors[counter++] = Face(this->voxel, static_cast<Direction::Value>(border));

            // second possible face only possible if the right voxel exists
            neighbors[counter++] = voxel[border_axis] == grid_check[border] ? none : Face(this->voxel + direction[border], this->direction);

            // third possible face only possible if right and upper voxel exists
            bool check = voxel[border_axis] == grid_check[border] || voxel[face_axis] != grid_check[this->direction];
            neighbors[counter++] = check ? none : Face(this->voxel + direction[border] + direction[this->direction], first_dir.flip());
        }

        return neighbors;
    }

    inline std::array<int64_t, 12> possibleBorderNeighborsId(Dims dims, int64_t none) const
    {
        auto neighbors = std::array<int64_t, 12>();
        auto size = dims.size();
        auto counter = 0;
        for(auto neighbor_face : possibleBorderNeighbors(dims, Face(size, Direction::RIGHT))){
            neighbors[counter++] = neighbor_face.startId() == size ? none : neighbor_face.toIndex();
        }
        return neighbors;
    }

    template <typename T>
    std::size_t
    neighbors(T& set, Dims dims, int64_t* neighbors) const
    {
        auto candidates = possibleBorderNeighborsId(dims, -1);
        // go through all candidates, every batch of 3 candidates has at most one face in set (otherwise not a manifold)
        std::size_t neighbors_counter = 0;
        std::size_t counter = 0;
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                // set.contains(candidates[counter])
                if (set.find(candidates[counter]) != set.end())
                {
                    neighbors[neighbors_counter++] = candidates[counter];
                    break;
                }
                ++counter;
            }
        }
        return neighbors_counter;
    }

    template <typename T>
    std::size_t
    neighbors(T& set, Dims dims, Face* neighbors) const
    {
        auto candidates = possibleBorderNeighbors(dims, Face(dims.size(), Direction::RIGHT));
        // go through all candidates, every batch of 3 candidates has at most one face in set (otherwise not a manifold)
        std::size_t neighbors_counter = 0;
        std::size_t counter = 0;
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                // set.contains(candidates[counter])
                if (set.find(candidates[counter]) != set.end())
                {
                    neighbors[neighbors_counter++] = candidates[counter];
                    break;
                }
                ++counter;
            }
        }
        return neighbors_counter;
    }


    /**
     * TODO: instead of uint64_t we want to return FieldLocation....
     * 
     * @brief return the indices of the corners of the face (the indices are the voxels with the corner as their lower left corner)
     *
     * @param dims: Must be the dimensions of the number of voxels
     * 
     * @return std::array<uint64_t,4>
     */
    std::array<uint64_t, 4> cornerIndices(Dims dims) const
    {
        const auto vertDim = dims.extend(1);
        const auto gridPos = Lattice::gridLocationFromCIndex(voxel, dims);
        const auto vertex = RawField<float>::FieldLocation(gridPos, vertDim).index();
        // const auto vertex = gridPositionToIndex(indexToGridPosition(voxel, dims), vertDim);

        std::array<uint64_t, 4> arr;
        if (direction == Direction::LEFT)
        {
            arr[0] = vertex + vertDim[0];
            arr[1] = vertex;
            arr[2] = vertex + vertDim[0] * vertDim[1];
            arr[3] = vertex + vertDim[0] + vertDim[0] * vertDim[1];
        }
        else if (direction == Direction::RIGHT)
        {
            arr[0] = vertex + 1 + vertDim[0] + vertDim[0] * vertDim[1];
            arr[1] = vertex + 1 + vertDim[0] * vertDim[1];
            arr[2] = vertex + 1;
            arr[3] = vertex + 1 + vertDim[0];
        }
        else if (direction == Direction::DOWN)
        {
            arr[0] = vertex + vertDim[0] * vertDim[1];
            arr[1] = vertex;
            arr[2] = vertex + 1;
            arr[3] = vertex + 1 + vertDim[0] * vertDim[1];
        }
        else if (direction == Direction::UP)
        {
            arr[0] = vertex + 1 + vertDim[0] + vertDim[0] * vertDim[1];
            arr[1] = vertex + 1 + vertDim[0];
            arr[2] = vertex + vertDim[0];
            arr[3] = vertex + vertDim[0] + vertDim[0] * vertDim[1];
        }
        else if (direction == Direction::BACKWARD)
        {
            arr[0] = vertex + 1 + vertDim[0];
            arr[1] = vertex + 1;
            arr[2] = vertex;
            arr[3] = vertex + vertDim[0];
        }
        else
        { //(direction == Direction::FORWARD)
            arr[0] = vertex + vertDim[0] + vertDim[0] * vertDim[1];
            arr[1] = vertex + vertDim[0] * vertDim[1];
            arr[2] = vertex + 1 + vertDim[0] * vertDim[1];
            arr[3] = vertex + 1 + vertDim[0] + vertDim[0] * vertDim[1];
        }
        return arr;
    }

    // VecSize corner, Dims dims, VecFloat voxelSize, VecFloat fieldPosition
    std::array<VecFloat, 4> cornerPosition(Lattice lattice){
        const auto gridPos = Lattice::gridLocationFromCIndex(voxel, lattice.dims());
        const auto pos = lattice.cornerPosition(gridPos);
        const auto voxelsize = lattice.voxelsize();

        if (direction == Direction::LEFT)
        {
            return {
                pos + VecInt::UP * voxelsize[1],
                pos,
                pos + VecInt::BACKWARD * voxelsize[2],
                pos + VecInt::UP * voxelsize[1] + VecInt::BACKWARD * voxelsize[2],
            };
        }
        else if (direction == Direction::RIGHT)
        {
            return {
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::UP * voxelsize[1] + VecInt::BACKWARD * voxelsize[2],
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::BACKWARD * voxelsize[2],
                pos + VecInt::RIGHT * voxelsize[0],
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::UP * voxelsize[1],
            };
        }
        else if (direction == Direction::DOWN)
        {
            return {
                pos + VecInt::BACKWARD * voxelsize[2],
                pos,
                pos + VecInt::RIGHT * voxelsize[0],
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::BACKWARD * voxelsize[2],
            };
        }
        else if (direction == Direction::UP)
        {
            return {
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::UP * voxelsize[1] + VecInt::BACKWARD * voxelsize[2],
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::UP * voxelsize[1],
                pos + VecInt::UP * voxelsize[1],
                pos + VecInt::UP * voxelsize[1] + VecInt::BACKWARD * voxelsize[2],
            };
        }
        else if (direction == Direction::BACKWARD)
        {
            return {
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::UP * voxelsize[1],
                pos + VecInt::RIGHT * voxelsize[0],
                pos,
                pos + VecInt::UP * voxelsize[1],
            };
        }
        else
        { //(direction == Direction::FORWARD)
            return {
                pos + VecInt::UP * voxelsize[1] + VecInt::BACKWARD * voxelsize[2],
                pos + VecInt::BACKWARD * voxelsize[2],
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::BACKWARD * voxelsize[2],
                pos + VecInt::RIGHT * voxelsize[0] + VecInt::UP * voxelsize[1] + VecInt::BACKWARD * voxelsize[2],
            };
        }
    }

    /**
     * @brief returns corners of both faces
     *
     * @return number of corners (at most 4)
     */
    std::size_t
    intersect(Face other, Dims dims, uint64_t* corners)
    {
        const VecInt gridPosition = static_cast<VecInt>(RawField<float>::FieldLocation(voxel).location(dims));
        const VecInt otherPosition = static_cast<VecInt>(RawField<float>::FieldLocation(other.voxel).location(dims));
        // const McVec3i gridPosition = indexToGridPosition(voxel, dims);
        // const McVec3i otherPosition = indexToGridPosition(other.voxel, dims);
        const VecInt diff = (gridPosition - otherPosition).abs();

        if (diff[0] > 1 || diff[1] > 1 || diff[2] > 1)
        {
            return 0;
        }

        const auto selfCorners = cornerIndices(dims);
        const auto otherCorners = other.cornerIndices(dims);

        std::size_t cnt = 0;

        for (auto corner : selfCorners)
        {
            for (auto otherCorner : otherCorners)
            {
                if (corner == otherCorner)
                {
                    corners[cnt++] = corner;
                }
            }
        }

        return cnt;
    }

    // /**
    //  * @brief returns an edge, if the intersection is exaclty the edge returned
    //  *
    //  * @param other
    //  * @param dims
    //  * @return ri::Option<Edge>
    //  */
    // std::optional<Edge>
    // intersectionEdge(Face other, Dims dims)
    // {
    //     std::array<uint64_t, 4> corners;
    //     auto cnt = intersect(other, dims, &corners[0]);
    //     if (cnt == 2)
    //     {
    //         // we may unwrap as all possible offsets are valid directions!
    //         return Edge(corners[0], offsetToDirection(Dims(dims[0] + 1, dims[1] + 1, dims[2] + 1), corners[1] - corners[0]).unwrap());
    //     }
    //     else
    //     {
    //         return ri::none;
    //     }
    // }

    // template <typename T>
    // Face(T voxel, T direction)
    //     : voxel(static_cast<uint64_t>(voxel))
    //     , direction(static_cast<uint64_t>(direction))
    // {}

    // hashing function defined at the end
};

template <typename Iter>
class FlipFaceDirectionConverter {
    using iterator_category = std::forward_iterator_tag;
    using value_type = Face;

    public:
        FlipFaceDirectionConverter(Iter iter, Dims dims) : m_iter(iter), m_dims(dims) {}

        Face operator*() const { return (*m_iter).flipDirection(m_dims); }
        Face operator->() { return (*m_iter).flipDirection(m_dims); }
        FlipFaceDirectionConverter& operator++() { ++m_iter; return *this; }
        FlipFaceDirectionConverter operator++(int) {FlipFaceDirectionConverter tmp = *this; ++(*this); return tmp;}

        friend bool operator== (const FlipFaceDirectionConverter& a, const FlipFaceDirectionConverter& b) { return a.m_iter == b.m_iter; }
        friend bool operator!= (const FlipFaceDirectionConverter& a, const FlipFaceDirectionConverter& b) { return a.m_iter != b.m_iter; }

    protected:
        Iter m_iter;
        Dims m_dims;
};

template <typename Iter>
class IndexToFaceConverter{
    using iterator_category = std::forward_iterator_tag;
    using value_type = Face;

    public:
        IndexToFaceConverter(Iter iter) : m_iter(iter) {}

        Face operator*() const { return Face::fromIndex(*m_iter); }
        Face operator->() { return Face::fromIndex(*m_iter); }
        IndexToFaceConverter& operator++() { ++m_iter; return *this; }
        IndexToFaceConverter operator++(int) {IndexToFaceConverter tmp = *this; ++(*this); return tmp;}

        FlipFaceDirectionConverter<IndexToFaceConverter<Iter>> flipDirection(Dims dims) { return FlipFaceDirectionConverter<IndexToFaceConverter<Iter>>(*this, dims); }

        friend bool operator== (const IndexToFaceConverter& a, const IndexToFaceConverter& b) { return a.m_iter == b.m_iter; }
        friend bool operator!= (const IndexToFaceConverter& a, const IndexToFaceConverter& b) { return a.m_iter != b.m_iter; }

    protected:
        Iter m_iter;
};

/**
 * @brief Like borderNeighbor4, but one does not need to input boolean-values. The Classifiert gets voxelIds and shall tell if we are inside or not.
 *
 * @tparam Classifier
 * @param classifier
 * @param dims
 * @param face
 * @param neighbors
 */
template <class Classifier>
void borderNeighbor4Classifier(Classifier classifier, Dims dims, Face face, Face* neighbors)
{
    const VecSize voxel = Lattice::gridLocationFromCIndex(face.voxel, dims);
    // const McVec3i voxel = indexToGridPosition(face.voxel, dims);
    const std::array<int64_t, 6> direction = {
        1,
        static_cast<int64_t>(dims[0]),
        static_cast<int64_t>(dims[0] * dims[1]),
        -1,
        -static_cast<int64_t>(dims[0]),
        -static_cast<int64_t>(dims[0] * dims[1]),
    };

    const std::array<std::size_t, 6> grid_check = {
        dims[0] - 1,
        dims[1] - 1,
        dims[2] - 1,
        0,
        0,
        0,
    };

    // Direction of our edges
    // std::array<int, 4> borders = { (face.direction + 4) % 6, (face.direction + 5) % 6, (face.direction + 1) % 6, (face.direction + 2) % 6 };
    // We follow right hand rule
    std::array<int, 4> borders = { (face.direction + 4) % 6, (face.direction + 2) % 6, (face.direction + 5) % 6, (face.direction + 1) % 6 };

    for (int neighbor_counter = 0; neighbor_counter < 4; ++neighbor_counter)
    {
        int border = borders[neighbor_counter];
        int border_axis = border % 3;
        int face_axis = face.direction % 3;

        if (voxel[border_axis] != grid_check[border]
            && voxel[face_axis] != grid_check[face.direction]
            && classifier(face.voxel + direction[border])
            && classifier(face.voxel + direction[border] + direction[face.direction]))
        {
            neighbors[neighbor_counter].voxel = face.voxel + direction[border] + direction[face.direction];
            neighbors[neighbor_counter].direction = static_cast<Direction::Value>((border + 3) % 6);
        }
        else if (voxel[border_axis] != grid_check[border] && classifier(face.voxel + direction[border]))
        {
            neighbors[neighbor_counter].voxel = face.voxel + direction[border];
            neighbors[neighbor_counter].direction = face.direction;
        }
        else
        {
            neighbors[neighbor_counter].voxel = face.voxel;
            neighbors[neighbor_counter].direction = static_cast<Direction::Value>(border);
        }
    }
}

namespace std
{
    template <>
    struct hash<Face>
    {
        std::size_t operator()(const Face& k) const
        {
            using boost::hash_combine;
            using boost::hash_value;

            // Start with a hash value of 0.
            std::size_t seed = 0;

            // Modify 'seed' by XORing and bit-shifting in
            // one member of 'Face' after the other:
            hash_combine(seed, hash_value(k.voxel));
            hash_combine(seed, hash_value(static_cast<uint8_t>(k.direction)));

            // Return the result.
            return seed;
        }
    };

} // namespace std