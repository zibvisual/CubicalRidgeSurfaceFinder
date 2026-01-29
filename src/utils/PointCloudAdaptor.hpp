#pragma once

template <typename T>
struct PointCloudAdaptor
{
    using coord_t = T::value_type;

    const std::vector<T>& m_points;

    PointCloudAdaptor(const std::vector<T>& points) : m_points(points){}

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const {
        return m_points.size();
    }

    // Returns the dim'th component of the idx'th point in the class:
    inline coord_t kdtree_get_pt(const size_t idx, const size_t dim) const {
        return m_points[idx][dim];
    }

    // Optional bounding-box computation: return false to default to a standard
    // bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned
    //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
    //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }
};