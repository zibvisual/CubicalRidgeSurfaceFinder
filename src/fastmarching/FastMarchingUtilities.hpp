#pragma once

#include <array>

#include <utils/Dims.hpp>
#include <utils/Vec.hpp>
#include <field/RawField.hpp>
// #include <field/RawField.hpp>

namespace fastmarching
{
    inline float
    fade(float t)
    {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }

    inline std::array<float, 3> sort(float x, float y, float z)
    {
        auto arr = std::array<float, 3>();
        if (x < y)
        {
            if (x < z)
            {
                arr[0] = x;
                if (y < z)
                {
                    arr[1] = y;
                    arr[2] = z;
                }
                else
                {
                    arr[1] = z;
                    arr[2] = y;
                }
            }
            else
            {
                arr[0] = z;
                arr[1] = x;
                arr[2] = y;
            }
        }
        else
        {
            if (x < z)
            {
                arr[0] = y;
                arr[1] = x;
                arr[2] = z;
            }
            else
            {
                arr[2] = x;
                if (y < z)
                {
                    arr[0] = y;
                    arr[1] = z;
                }
                else
                {
                    arr[0] = z;
                    arr[1] = y;
                }
            }
        }
        return arr;
    }

    // Instead of mapping view, we should use field
    template <class MappingView>
    std::array<float, 3>
    findMinDifferential_view(VecInt pos, Dims dims, const MappingView &values)
    {
        // TODO: instead of VecInt we want to use the location (depending on MappingView)
        auto top_val = values.get_optional(RawField<float>::FieldLocation(pos + VecInt::UP, dims).index()).value_or(INFINITY);
        auto bottom_val = values.get_optional(RawField<float>::FieldLocation(pos + VecInt::DOWN, dims).index()).value_or(INFINITY);
        auto x_val = std::min(top_val, bottom_val);

        auto left_val = values.get_optional(RawField<float>::FieldLocation(pos + VecInt::LEFT, dims).index()).value_or(INFINITY);
        auto right_val = values.get_optional(RawField<float>::FieldLocation(pos + VecInt::RIGHT, dims).index()).value_or(INFINITY);
        auto y_val = std::min(left_val, right_val);

        auto front_val = values.get_optional(RawField<float>::FieldLocation(pos + VecInt::FORWARD, dims).index()).value_or(INFINITY);
        auto back_val = values.get_optional(RawField<float>::FieldLocation(pos + VecInt::BACKWARD, dims).index()).value_or(INFINITY);
        auto z_val = std::min(front_val, back_val);

        // sort neighbors
        return sort(x_val, y_val, z_val);
    }

    // Instead of mapping view, we should use field
    template <class Field>
    std::array<float, 3>
    findMinDifferential_field(const Field& field, typename Field::Location pos)
    {
        auto neighbors = field.gridNeighbors6(pos);

        auto top_val = field.get(neighbors[0]).value_or(INFINITY);
        auto bottom_val = field.get(neighbors[1]).value_or(INFINITY);
        auto x_val = std::min(top_val, bottom_val);

        auto left_val = field.get(neighbors[2]).value_or(INFINITY);
        auto right_val = field.get(neighbors[3]).value_or(INFINITY);
        auto y_val = std::min(left_val, right_val);

        auto front_val = field.get(neighbors[4]).value_or(INFINITY);
        auto back_val = field.get(neighbors[5]).value_or(INFINITY);
        auto z_val = std::min(front_val, back_val);

        // sort neighbors
        return sort(x_val, y_val, z_val);
    }

    inline float
    calculate_speed(float pot, float min_threshold, float max_threshold)
    {
        // clang-format off
        float norm =    pot < min_threshold ? 0.0f : 
                        pot > max_threshold ? 1.0f : 
                        (pot - min_threshold) / (max_threshold - min_threshold);
        // clang-format on
        const float speed = (1 - norm) / norm;
        return speed;
    }

    namespace {
        template <class Field>
        float eikonal(const Field &potential, typename Field::Location pos, float min_threshold, float max_threshold, std::array<float, 3> neighbors)
        {
            // can be get_unchecked
            const float speed = calculate_speed(potential.get(pos).value(), min_threshold, max_threshold);

            // 3-dim update
            if (neighbors[2] != INFINITY)
            {
                // to reduce numerical errors we use differences of values to mean
                double mean = (neighbors[0] + neighbors[1] + neighbors[2]) / 3.0;
                double quad = 0;
                for (int i = 0; i < 3; ++i)
                {
                    double cent = neighbors[i] - mean;
                    quad += cent * cent;
                }
                double nom = (speed * speed - quad) / 3.0;
                double solution = nom > 0 ? (mean + sqrt(nom)) : mean;
                // test if solution is bigger than all participants
                bool bigger = true;
                for (int j = 0; j < 3; ++j)
                {
                    if (solution < neighbors[j])
                    {
                        bigger = false;
                    }
                }
                if (bigger)
                {
                    return (float)solution;
                }
            }

            if (neighbors[1] != INFINITY)
            {
                // 2-dim update
                float left = neighbors[0];
                float right = neighbors[1];
                // to reduce numerical errors we use differences of values to mean
                double mean = (left + right) / 2.0;
                float center_left = left - mean;
                float center_right = right - mean;
                double quad = center_left * center_left + center_right * center_right;
                double nom = (speed * speed - quad) / 2.0;
                double solution = nom > 0 ? (mean + sqrt(nom)) : mean;
                if (solution >= left && solution >= right)
                {
                    return solution;
                }
            }

            // 1-dim update
            return neighbors[0] + speed;
        }

        inline float
        euclidean(std::array<float, 3> neighbors)
        {
            // 3-dim update
            if (neighbors[2] != INFINITY)
            {
                // to reduce numerical errors we use differences of values to mean
                double mean = (neighbors[0] + neighbors[1] + neighbors[2]) / 3.0;
                double quad = 0;
                for (int i = 0; i < 3; ++i)
                {
                    double cent = neighbors[i] - mean;
                    quad += cent * cent;
                }
                double nom = (1 - quad) / 3.0;
                double solution = nom > 0 ? (mean + sqrt(nom)) : mean;
                // test if solution is bigger than all participants
                bool bigger = true;
                for (int j = 0; j < 3; ++j)
                {
                    if (solution < neighbors[j])
                    {
                        bigger = false;
                    }
                }
                if (bigger)
                {
                    return (float)solution;
                }
            }

            // 2-dim update
            if (neighbors[1] != INFINITY)
            {
                float left = neighbors[0];
                float right = neighbors[1];
                // to reduce numerical errors we use differences of values to mean
                double mean = (left + right) / 2.0;
                float center_left = left - mean;
                float center_right = right - mean;
                double quad = center_left * center_left + center_right * center_right;
                double nom = (1 - quad) / 2.0;
                double solution = nom > 0 ? (mean + sqrt(nom)) : mean;
                if (solution >= left && solution >= right)
                {
                    return solution;
                }
            }

            // 1-dim update
            return neighbors[0] + 1.0;
        }
    }

    /**
     * @brief
     *
     * @tparam MappingView
     * @param nodeId
     * @param dims
     * @param potential
     * @param min_threshold
     * @param max_threshold
     * @param timeData
     * @return float
     *
     * We use the Godunov’s scheme or  Rouy-Tourin’s scheme to update our field. Both make use of finite differences to estimate the gradient.
     * The gradient estimation of the fime field should coincide with the given scheme. See "Fast Marching Methods - parallel implementation and analysis"
     * for more information on the different schemes.
     *
     * We might want to use a different scheme to improve accuracy.
     * Also the performance of this function should have a big performance influence on the FM algorithm. Optimizing it is probably a good idea.
     *
     * In this function, to improve numerical accuracy, we calculate the mean. That is, instead of calculating the new value as
     * \sum_i c_i + \sqrt((\sum_i c_i)^2) - n \sum_i c_i^2 + nf^2
     * we calculate
     * k + \sqrt(1/n * (f^2 - \sum_i b_i^2)) with k := (\sum_i c_i) / n and b_i := c_i - k
     */
    template <class Field, class MappingView>
    float
    eikonal_view(const Field &potential, typename Field::Location pos, float min_threshold, float max_threshold, const MappingView &timeData)
    {
        auto neighbors = findMinDifferential_view(static_cast<VecInt>(pos.location(potential.dims())), potential.dims(), timeData);
        return eikonal(potential, pos, min_threshold, max_threshold, neighbors);
    }

    template <class Field>
    float
    eikonal_field(const Field &potential, float min_threshold, float max_threshold, const Field &timeData, typename Field::Location pos)
    {
        auto neighbors = findMinDifferential_field(timeData, pos);
        return eikonal(potential, pos, min_threshold, max_threshold, neighbors);
    }

    template <class MappingView>
    float
    euclidean_view(VecInt pos, Dims dims, const MappingView &distanceData)
    {
        auto neighbors = findMinDifferential_view(pos, dims, distanceData);
        return euclidean(neighbors);
    }

    template <class Field>
    float
    euclidean_field(const Field &distanceData, typename Field::Location pos)
    {
        auto neighbors = findMinDifferential_field(distanceData, pos);
        return euclidean(neighbors);
    }

    template <class Field>
    VecFloat
    gradient(const Field &field, typename Field::Location loc)
    {
        // we want to find the 3 smallest neighbor fields of timedata.
        const typename Field::Location offsets[6] = {
            field.moveLocation(loc, Direction::LEFT),
            field.moveLocation(loc, Direction::RIGHT),
            field.moveLocation(loc, Direction::UP),
            field.moveLocation(loc, Direction::DOWN),
            field.moveLocation(loc, Direction::FORWARD),
            field.moveLocation(loc, Direction::BACKWARD),
        };
        VecFloat grad = VecFloat(0.0f, 0.0f, 0.0f);
        auto value_opt = field.get(loc);
        if (!value_opt.has_value())
        {
            return grad;
        }
        float value = value_opt.value();
        for (int i = 0; i < 3; ++i)
        {
            auto bottom = offsets[i * 2];
            auto top = offsets[i * 2 + 1];

            if (top.valid() && field.get(top).value() < value)
            {
                if (bottom.valid() && field.get(bottom).value() < field.get(top).value())
                {
                    grad[i] = field.get(bottom).value() - value;
                }
                else
                {
                    grad[i] = value - field.get(top).value();
                }
            }
            else if (bottom.valid() && field.get(bottom).value() < value)
            {
                grad[i] = field.get(bottom).value() - value;
            }
        }
        return grad;
    }

    template <class Field>
    VecFloat
    normalized_gradient(const Field &field, typename Field::Location loc)
    {
        auto grad = fastmarching::gradient(field, loc);
        grad.normalize();
        return grad;
    }

    // template <class Field>
    // McVec3f
    // gradient(McVec3f point, const MappingView &values, McDim3l dims, McBox3f bbox, McVec3f voxelsize)
    // {
    //     const mcint64 offsets[8] = {
    //         0,
    //         1,
    //         dims[0],
    //         dims[0] + 1,
    //         dims[0] * dims[1],
    //         dims[0] * dims[1] + 1,
    //         dims[0] * dims[1] + dims[0],
    //         dims[0] * dims[1] + dims[0] + 1};

    //     // TODO: first we would need futil::worldToNearestCornerPosition
    //     McVec3f position = futil::worldToCornerPosition(point, voxelsize, bbox);
    //     // it is not clear if all 8 gradients exist. so we check first if the given corners even are all existing
    //     auto pair = mutil::floor_remainder(position);
    //     McVec3i corner = pair.first;
    //     McVec3f pre_weights = pair.second;
    //     mcint64 cornerId = futil::gridPositionToIndex(corner, dims);

    //     float sum_weights = 0.0f;
    //     float pre_weight;
    //     McVec3f sum_val = McVec3f(0.0f, 0.0f, 0.0f);
    //     mcint64 id;
    //     // first corner
    //     id = cornerId + offsets[0];
    //     if (values.get(id) != values.getDefault())
    //     {
    //         pre_weight = fade(std::min({1.0f - pre_weights[0], 1.0f - pre_weights[1], 1.0f - pre_weights[2]}));
    //         sum_val += pre_weight * fastmarching::gradient(id, dims, values);
    //         sum_weights += pre_weight;
    //     }
    //     // second corner
    //     id = cornerId + offsets[1];
    //     if (values.get(id) != values.getDefault())
    //     {
    //         pre_weight = fade(std::min({pre_weights[0], 1.0f - pre_weights[1], 1.0f - pre_weights[2]}));
    //         sum_val += pre_weight * fastmarching::gradient(id, dims, values);
    //         sum_weights += pre_weight;
    //     }
    //     // third corner
    //     id = cornerId + offsets[2];
    //     if (values.get(id) != values.getDefault())
    //     {
    //         pre_weight = fade(std::min({1.0f - pre_weights[0], pre_weights[1], 1.0f - pre_weights[2]}));
    //         sum_val += pre_weight * fastmarching::gradient(id, dims, values);
    //         sum_weights += pre_weight;
    //     }
    //     // fourth corner
    //     id = cornerId + offsets[3];
    //     if (values.get(id) != values.getDefault())
    //     {
    //         pre_weight = fade(std::min({pre_weights[0], pre_weights[1], 1.0f - pre_weights[2]}));
    //         sum_val += pre_weight * fastmarching::gradient(id, dims, values);
    //         sum_weights += pre_weight;
    //     }
    //     // fifth corner
    //     id = cornerId + offsets[4];
    //     if (values.get(id) != values.getDefault())
    //     {
    //         pre_weight = fade(std::min({1.0f - pre_weights[0], 1.0f - pre_weights[1], pre_weights[2]}));
    //         sum_val += pre_weight * fastmarching::gradient(id, dims, values);
    //         sum_weights += pre_weight;
    //     }
    //     // sixth corner
    //     id = cornerId + offsets[5];
    //     if (values.get(id) != values.getDefault())
    //     {
    //         pre_weight = fade(std::min({pre_weights[0], 1.0f - pre_weights[1], pre_weights[2]}));
    //         sum_val += pre_weight * fastmarching::gradient(id, dims, values);
    //         sum_weights += pre_weight;
    //     }
    //     // seventh corner
    //     id = cornerId + offsets[6];
    //     if (values.get(id) != values.getDefault())
    //     {
    //         pre_weight = fade(std::min({1.0f - pre_weights[0], pre_weights[1], pre_weights[2]}));
    //         sum_val += pre_weight * fastmarching::gradient(id, dims, values);
    //         sum_weights += pre_weight;
    //     }
    //     // eigth corner
    //     id = cornerId + offsets[7];
    //     if (values.get(id) != values.getDefault())
    //     {
    //         pre_weight = fade(std::min({pre_weights[0], pre_weights[1], pre_weights[2]}));
    //         sum_val += pre_weight * fastmarching::gradient(id, dims, values);
    //         sum_weights += pre_weight;
    //     }

    //     return sum_weights > 0.0f ? sum_val / sum_weights : sum_val;
    // }
}
