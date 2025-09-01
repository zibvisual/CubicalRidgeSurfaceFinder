#pragma once

#include <vector>
#include <array>
#include <limits>
#include <optional>
#include <filesystem>
#include <utils/Dims.hpp>
#include <utils/BBox.hpp>
#include <utils/Lattice.hpp>
#include <utils/Vec.hpp>
#include <io/npy.hpp>
#include <io/metadata.hpp>
#include <utils/Direction.hpp>
#include <io/spatialinfo.hpp>

/**
 * Fields are just values in some 3D space. If no value exist at a spot, we return None or some default value.
 * If we are out of bounds, we also return None or some default.
 * We might want to discern between OutOfBox and NoValue. -> Field only knows NoValue and Value
 * a bounding box has a default background value if pos is inside box.
 * Neighbors are calculated with the help of the bounding box to check if they are still inside the field?
 */


/**
 * Each Field may use its own FieldLocation. Each Field is responsible for indexing, bounds checking and moving locations
 */
template <typename T>
class RawField
{
    // TODO: Tiled Field to allow for better cache locality? One Tile is of size 8x8x8?
public:
    class FieldLocation
    {
    public:
        FieldLocation(VecSize location, Dims dims)
            : m_index(location.z() * dims.x() * dims.y() + location.y() * dims.x() + location.x())
        {
            // check if we are in bounds!
            if(location.x() < 0 || location.x() >= dims.x() || location.y() < 0 || location.y() >= dims.y() || location.z() < 0 || location.z() >= dims.z())
            {
                m_index = std::numeric_limits<uint64_t>::max();
            }else{
                m_index = location.z() * dims.x() * dims.y() + location.y() * dims.x() + location.x();
            }
        }

        FieldLocation(VecInt location, Dims dims)
        : m_index(std::numeric_limits<uint64_t>::max())
        {
            // check if we are in bounds!
            if(location.x() < 0 || location.x() >= dims.x() || location.y() < 0 || location.y() >= dims.y() || location.z() < 0 || location.z() >= dims.z())
            {
                m_index = std::numeric_limits<uint64_t>::max();
            }else{
                m_index = static_cast<uint64_t>(location.z()) * dims.x() * dims.y() + static_cast<uint64_t>(location.y()) * dims.x() + static_cast<uint64_t>(location.x());
            }
        }

        FieldLocation(uint64_t index)
            : m_index(index)
        {
        }

        // Empty FieldLocation (out of bounds)
        FieldLocation()
            : m_index(std::numeric_limits<uint64_t>::max())
        {
        }

        bool operator==(const FieldLocation& rhs) const 
        {
            return m_index == rhs.m_index;
        }

        // THIS IS DANGEROUS (for empty field locations!)
        VecSize location(Dims dims) const
        {
            const uint64_t z = m_index / (dims.y() * dims.x());
            const uint64_t rest = m_index - z * (dims.y() * dims.x());
            const uint64_t y = rest / dims.x();
            const uint64_t x = rest - y * dims.x();
            return VecSize(x, y, z);
        }

        uint64_t index() const
        {
            return m_index;
        }

        bool empty() const
        {
            return m_index == std::numeric_limits<uint64_t>::max();
        }

        bool valid() const
        {
            return m_index != std::numeric_limits<uint64_t>::max();
        }

    private:
        uint64_t m_index;
    };

    using Location = FieldLocation;

    RawField(Dims dims)
    : m_data(std::vector<T>(dims.size(), T()))
    , m_lattice(dims)
    , m_dirty_vals(true)
    , m_min_val(std::numeric_limits<T>::max())
    , m_max_val(std::numeric_limits<T>::lowest())
    , m_undefined_val(std::numeric_limits<T>::max())
    {
    }

    /**
     * data must be contiguous with x being the fastest iteration
     */
    RawField(std::vector<T> data, Dims dims) 
    : m_data(data)
    , m_lattice(dims)
    , m_dirty_vals(true)
    , m_min_val(std::numeric_limits<T>::max())
    , m_max_val(std::numeric_limits<T>::lowest())
    , m_undefined_val(std::numeric_limits<T>::max())
    {
    }

    /**
     * data must be contiguous with x being the fastest iteration
     */
    RawField(std::vector<T> data, Lattice lattice) 
    : m_data(data)
    , m_lattice(lattice)
    , m_dirty_vals(true)
    , m_min_val(std::numeric_limits<T>::max())
    , m_max_val(std::numeric_limits<T>::lowest())
    , m_undefined_val(std::numeric_limits<T>::max())
    {
    }

    void save(std::string output_path){
        npy::npy_data_ptr<T> data;
        data.data_ptr = this->data();
        data.shape = npy::shape_t(m_lattice.dims()[0], m_lattice.dims()[1], m_lattice.dims()[2]);
        // fortran order is also inside of data (false currently...)
        npy::write_npy(output_path, data);

        const auto voxelsize = m_lattice.voxelsize();
        const auto origin = m_lattice.origin();

        metadata_t metadata;
        metadata.keywords.push_back("voxelsize");
        metadata.data.push_back({std::to_string(voxelsize[0]), std::to_string(voxelsize[1]), std::to_string(voxelsize[2])});
        metadata.keywords.push_back("origin");
        metadata.data.push_back({std::to_string(origin[0]), std::to_string(origin[1]), std::to_string(origin[2])});
        write_metadata(output_path, metadata);
    }

    static RawField load(std::string input_path){
        npy::npy_data<T> data = npy::read_npy<T>(input_path);
        if (data.shape.ndims() != 3)
        {
            throw std::invalid_argument("input image must be 3 dimensional");
        }
        auto dims = Dims(data.shape[2], data.shape[1], data.shape[0]);

        Lattice lattice = Lattice(dims);
        try {
            metadata_t metadata = read_metadata(input_path);
            // std::cout << metadata << std::endl;
            lattice = SpatialInformation::fromMetadata(metadata).toLattice(dims);
        }
        catch (const std::ios_base::failure& e)
        {
            // ignore missing metadata file
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            throw std::runtime_error("Error: metadata file parsing failed!");
        }
    
        return RawField(data.data, lattice);
    }
  
    /**
     * Returns index for the given location (might be empty)
     */
    FieldLocation createLocation(VecSize location) const
    {
        if (m_lattice.dims().contains(location))
        {
            return FieldLocation(location, m_lattice.dims());
        }
        else
        {
            return FieldLocation();
        }
    }

    /**
     * Returns index for the given location (might be empty)
     */
    FieldLocation createLocation(VecInt location) const
    {
        if (m_lattice.dims().contains(location))
        {
            return FieldLocation(location, m_lattice.dims());
        }
        else
        {
            return FieldLocation();
        }
    }

    /**
     * Returns index for the given location (might be empty)
     */
    FieldLocation createLocation(VecFloat location) const
    {
        // calculate location
        auto grid = m_lattice.gridLocation(location);
        return createLocation(grid);
    }

    /**
     * Returns moved location (might be empty if out of bounds)
     */
    FieldLocation moveLocation(FieldLocation location, Direction direction) const
    {
        auto loc = location.location(m_lattice.dims());
        switch (direction)
        {
        case Direction::LEFT:
            if (loc[0] == 0)
            {
                return FieldLocation();
            }
            return FieldLocation(location.index() - 1);
        case Direction::RIGHT:
            if (loc[0] + 1 == m_lattice.dims()[0])
            {
                return FieldLocation();
            }
            return FieldLocation(location.index() + 1);
        case Direction::UP:
            if (loc[1] == 0)
            {
                return FieldLocation();
            }
            return FieldLocation(location.index() - m_lattice.dims()[0]);
        case Direction::DOWN:
            if (loc[1] + 1 == m_lattice.dims()[1])
            {
                return FieldLocation();
            }
            return FieldLocation(location.index() + m_lattice.dims()[0]);
        case Direction::FORWARD:
            if (loc[2] == 0)
            {
                return FieldLocation();
            }
            return FieldLocation(location.index() - m_lattice.dims()[0] * m_lattice.dims()[1]);
        case Direction::BACKWARD:
            if (loc[2] + 1 == m_lattice.dims()[2])
            {
                return FieldLocation();
            }
            return FieldLocation(location.index() + m_lattice.dims()[0] * m_lattice.dims()[1]);
        default:
            return FieldLocation();
        }
    }

    std::array<FieldLocation, 6> gridNeighbors6(FieldLocation location) const
    {
        // do all cases of moveLocation
        auto loc = location.location(m_lattice.dims());
        std::array<FieldLocation, 6> res = {FieldLocation(), FieldLocation(), FieldLocation(), FieldLocation(), FieldLocation(), FieldLocation()};
        if (loc[0] != 0)
        {
            res[0] = FieldLocation(location.index() - 1);
        }
        if (loc[0] + 1 != m_lattice.dims()[0])
        {
            res[1] = FieldLocation(location.index() + 1);
        }
        if (loc[1] != 0)
        {
            res[2] = FieldLocation(location.index() - m_lattice.dims()[0]);
        }
        if (loc[1] + 1 != m_lattice.dims()[1])
        {
            res[3] = FieldLocation(location.index() + m_lattice.dims()[0]);
        }
        if (loc[2] != 0)
        {
            res[4] = FieldLocation(location.index() - m_lattice.dims()[0] * m_lattice.dims()[1]);
        }
        if (loc[2] + 1 != m_lattice.dims()[2])
        {
            res[5] = FieldLocation(location.index() + m_lattice.dims()[0] * m_lattice.dims()[1]);
        }
        return res;
    }

    /**
     * Returns underlying value or None if out of bounds.
     */
    std::optional<T> get(VecSize location) const
    {
        if (m_lattice.dims().contains(location))
        {
            return get_unchecked(FieldLocation(location, m_lattice.dims()));
        }
        else
        {
            return std::optional<T>();
        }
    }

    /**
     * Returns underlying value or None if out of bounds.
     */
    std::optional<T> get(VecInt location) const
    {
        if (m_lattice.dims().contains(location))
        {
            return get_unchecked(FieldLocation(location, m_lattice.dims()));
        }
        else
        {
            return std::optional<T>();
        }
    }

    /**
     * Returns underlying value or NONE if out of bounds.
     */
    std::optional<T> get(FieldLocation location) const
    {
        if (location.index() >= m_data.size())
        {
            return std::optional<T>();
        }
        return m_data[location.index()];
    }

    /**
     * Returns underlying value. Does NOT check for out of bounds.
     */
    T get_unchecked(FieldLocation location) const
    {
        return m_data[location.index()];
    }


    /** --------------------------------------------------------------------------- */

    std::pair<T, T> getRange() const
    {
        if(m_dirty_vals){
            calculateRange();
        }
        return std::pair<T,T>(m_min_val, m_max_val);
    }

    Dims getDims() const
    {
        return m_lattice.dims();
    }

    Dims
    dims() const
    {
        return m_lattice.dims();
    }

    VecFloat getVoxelSize() const
    {
        return m_lattice.voxelsize();
    }

    CornerBBox getCornerBoundingBox() const
    {
        return m_lattice.cornerbox();
    }

    CenterBBox getCenterBoundingBox() const
    {
        return m_lattice.centerbox();
    }

    const Lattice& lattice() const
    {
        return m_lattice;
    }

    T* data()
    {
        return m_data.data();
    }

    void setBBox(CornerBBox bbox) {
        m_lattice.setBBox(bbox);
    }

    void setOrigin(VecFloat orig){
        m_lattice.setOrigin(orig);
    }

    void setVoxelSize(VecFloat voxelsize){
        m_lattice.setVoxelSize(voxelsize);
    }
    

private:
    void calculateRange() const
    {
        m_min_val = std::numeric_limits<float>::max();
        m_max_val = std::numeric_limits<float>::lowest();
        for(auto i = 0; i < m_data.size(); ++i){
            if(m_data[i] != m_undefined_val && m_data[i] < m_min_val){
                m_min_val = m_data[i];
            }
            if(m_data[i] != m_undefined_val && m_data[i] > m_max_val){
                m_max_val = m_data[i];
            }
        }
        m_dirty_vals = false;
    }

private:
    std::vector<T> m_data;
    Lattice m_lattice;

    // cached values
    mutable bool m_dirty_vals;
    mutable T m_min_val;
    mutable T m_max_val;
    mutable T m_undefined_val;

};