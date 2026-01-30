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
#include <field/RawFieldView.hpp>

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
    using Location = typename RawFieldView<T>::Location;
    using FieldLocation = typename RawFieldView<T>::FieldLocation;

    RawField()
    : m_data()
    , m_view(m_data.data(), Lattice())
    {}

    RawField(Dims dims)
    : m_data(std::vector<T>(dims.size(), T()))
    , m_view(m_data.data(), dims)
    {
        //TODO: instead of dirty vals we could just update m_min_val and m_max_val
    }

    RawField(Lattice lattice)
    : m_data(std::vector<T>(lattice.dims().size(), T()))
    , m_view(m_data.data(), lattice)
    {
        //TODO: instead of dirty vals we could just update m_min_val and m_max_val
    }

    /**
     * data must be contiguous with x being the fastest iteration
     */
    RawField(std::vector<T> data, Dims dims) 
    : m_data(data)
    , m_view(m_data.data(), dims)
    {
    }

    /**
     * data must be contiguous with x being the fastest iteration
     */
    RawField(std::vector<T> data, Lattice lattice) 
    : m_data(data)
    , m_view(m_data.data(), lattice)
    {
    }

    /**
     * Get a view into RawField
     */
    const RawFieldView<T>& get_view() const{
        return m_view;
    }

    void save(std::string output_path){
        m_view.save(output_path);
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
            metadata_t metadata = read_metadata(input_path, true);
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
    Location createLocation(VecSize location) const
    {
        return m_view.createLocation(location);
    }

    /**
     * Returns index for the given location (might be empty)
     */
    Location createLocation(VecInt location) const
    {
        return m_view.createLocation(location);
    }

    /**
     * Returns index for the given location (might be empty)
     */
    Location createLocation(VecFloat location) const
    {
        return m_view.createLocation(location);
    }

    /**
     * Returns moved location (might be empty if out of bounds)
     */
    Location moveLocation(Location location, Direction direction) const
    {
        return m_view.moveLocation(location, direction);
    }

    std::array<Location, 6> gridNeighbors6(Location location) const
    {
        return m_view.gridNeighbors6(location);
    }

    /**
     * Returns underlying value or None if out of bounds.
     */
    std::optional<T> get(VecSize location) const
    {
        return m_view.get(location);
    }

    /**
     * Returns underlying value or undefined value if out of bounds.
     */
    T get_or(VecSize location) const
    {
        return m_view.get_or(location);
    }

    /**
     * Returns underlying value or given value if out of bounds.
     */
    T get_or(VecSize location, const T value) const
    {
        return m_view.get_or(location, value);
    }

    T get_unchecked(VecSize location) const
    {
        return m_view.get_unchecked(location);
    }

    bool set(VecSize location, const T value)
    {
        return m_view.set(location, value);
    }

    void set_unchecked(VecSize location, const T value)
    {
        m_view.set_unchecked(location, value);
    }

    /**
     * Returns underlying value or None if out of bounds.
     */
    std::optional<T> get(VecInt location) const
    {
        return m_view.get(location);
    }

    /**
     * Returns underlying value or undefined value if out of bounds.
     */
    T get_or(VecInt location) const
    {
        return m_view.get_or(location);
    }

    /**
     * Returns underlying value or undefined value if out of bounds.
     */
    T get_or(VecInt location, const T value) const
    {
        return m_view.get_or(location, value);
    }

    T get_unchecked(VecInt location) const
    {
        return m_view.get_unchecked(location);
    }

    /**
     * Returns underlying value or NONE if out of bounds.
     */
    std::optional<T> get(Location location) const
    {
        return m_view.get(location);
    }

    /**
     * Returns underlying value. Does NOT check for out of bounds.
     */
    T get_unchecked(Location location) const
    {
        return m_view.get_unchecked(location);
    }


    /** --------------------------------------------------------------------------- */

    std::pair<T, T> getRange()
    {
        return m_view.getRange();
    }

    Dims getDims() const
    {
        return m_view.dims();
    }

    Dims
    dims() const
    {
        return m_view.dims();
    }

    /**
     * Resize field and fill with undefined value if necessary
     */
    void resize(Dims dims) {
        m_data.resize(dims.size(), m_view.getUndefinedValue());
        m_view.resizedData(m_data.data(), dims, m_view.getUndefinedValue());    
    }

    void resize(Dims dims, T fillValue) {
        m_data.resize(dims.size(), fillValue);
        m_view.resizedData(m_data.data(), dims, fillValue);
    }

    /**
     * Set all values to undefined Value
     */
    void fill() {
        m_view.fill();
    }

    void fill(T fillValue) {
        m_view.fill(fillValue);
    }

    void setUndefinedValue(T undefinedValue) {
        m_view.setUndefinedValue(undefinedValue);
    }

    T getUndefinedValue(){
        return m_view.getUndefinedValue();
    }

    VecFloat getVoxelSize() const
    {
        return m_view.lattice().voxelsize();
    }

    CornerBBox getCornerBoundingBox() const
    {
        return m_view.lattice().cornerbox();
    }

    CenterBBox getCenterBoundingBox() const
    {
        return m_view.lattice().centerbox();
    }

    const Lattice& lattice() const
    {
        return m_view.lattice();
    }

    T* data()
    {
        return m_data.data();
    }

    const T* data() const
    {
        return m_data.data();
    }

    void setBBox(CornerBBox bbox) {
        m_view.setBBox(bbox);
    }

    void setOrigin(VecFloat orig){
        m_view.setOrigin(orig);
    }

    void setVoxelSize(VecFloat voxelsize){
        m_view.setVoxelSize(voxelsize);
    }

private:
    std::vector<T> m_data;
    RawFieldView<T> m_view;
};