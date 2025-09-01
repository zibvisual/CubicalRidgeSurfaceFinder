#pragma once

#include <vector>
#include <limits>
#include <cmath>
#include <cassert>
#include <algorithm>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Dims.hpp"
#include "Vec.hpp"

inline std::vector<float> gauss_kernel(float sigma, const int radius, const int order)
{
    // Arbitrary orders require recursive computation of Hermite polynomials.
    // Order up to 2 is enough for our applications.
    assert(order == 0 || order == 1 || order == 2);

    sigma = std::abs(sigma);
    sigma = std::max(sigma, std::numeric_limits<float>::epsilon());

    float sigma2 = sigma * sigma;

    std::vector<float> kernel(2 * radius + 1, 0.0f);
    float sum = 0.0f;

    for (int i = 0; i < 2 * radius + 1; ++i)
    {
        float dx = i - radius;
        float value = std::exp(dx * dx / (-2.0f * sigma2));
        sum += value;

        // I flipped some signs here (experimentally) so that the results
        // match with the stencil methods. Make sure I actually compute the
        // convolution above and that the derivatives are correct.
        if (order == 1)
        {
            value = value * (dx / sigma2);
        }
        else if (order == 2)
        {
            value = value * (dx * dx / sigma2 - 1.0f) / sigma2;
        }
        kernel[i] = value;
    }

    // We just normalize so that the 0-th order Gaussian kernel
    // has unit area. Every derivative is thus normalized wrt. to the
    // discretized 0th order kernel.
    for (int i = 0; i < 2 * radius + 1; ++i)
    {
        kernel[i] /= sum;
    }
    return kernel;
}

/**
 * @brief 1d convolution of a specific window
 *
 * @param data
 * @param output
 * @param dims the dimensions of the given vector (is usually bigger than the wanted window)
 * @param startVoxel startVoxel of the window (entrywise smaller than endVoxel)
 * @param endVoxel not inclusive!
 * @param direction for x->1, y->dims[0] and z->dims[0]*dims[1]
 * @param kernel
 */
inline void conv_1d(const float* data, float* output, const Dims dims, const VecInt startVoxel, const VecInt endVoxel, const int direction, std::vector<float> kernel)
{
    const int radius = kernel.size() / 2;

    for (int iz = startVoxel[2]; iz < endVoxel[2]; ++iz)
    {
        for (int iy = startVoxel[1]; iy < endVoxel[1]; ++iy)
        {
            for (int ix = startVoxel[0]; ix < endVoxel[0]; ++ix)
            {
                const int counter = ix + iy * dims[0] + iz * dims[0] * dims[1];
                // convolution at point ix,iy,iz
                float sum = 0.0f;
                for (int i = -radius; i <= radius; ++i)
                {
                    sum += kernel[i + radius] * data[counter + i * direction];
                }
                output[counter] = sum;
            }
        }
    }
}

/**
 * @brief
 *
 * @param data
 * @param dims
 * @param startVoxel where to start the convolution
 * @param endVoxel where to end the convolution
 * @param sigma
 * @param radius
 * @param order
 * @return std::vector<float> a vector with the results.
 * This will be smaller than the window given by start- and endVoxel, as convolution shrinks the size of the window.
 */
inline std::vector<float> gauss_filter(const std::vector<float> data, const Dims dims, VecInt startVoxel, VecInt endVoxel, const float sigma, const int radius, const VecInt order)
{
    auto output = std::vector<float>(dims.size());
    // PERF: one could use a heap array for temp, no vector necessary
    auto temp = std::vector<float>(dims.size());

    auto windowStart = VecInt(radius, 0, 0);
    auto windowEnd = VecInt(dims[0] - radius, dims[1], dims[2]);
    auto kernelx = gauss_kernel(sigma, radius, order[0]);
    conv_1d(&data[0], &output[0], dims, windowStart, windowEnd, 1, kernelx);
    windowStart[1] += radius;
    windowEnd[1] -= radius;
    auto kernely = gauss_kernel(sigma, radius, order[1]);
    conv_1d(&output[0], &temp[0], dims, windowStart, windowEnd, dims[0], kernely);
    windowStart[2] += radius;
    windowEnd[2] -= radius;
    auto kernelz = gauss_kernel(sigma, radius, order[2]);
    conv_1d(&temp[0], &output[0], dims, windowStart, windowEnd, dims[0] * dims[1], kernelz);
    return output;
}

/**
 * @brief convolution at one specific point
 *
 * @tparam T Any View which allows reading values by indexing
 * @param data
 * @param dims
 * @param voxel
 * @param kernel
 * @param radius of the kernel, such that (radius*2+1)^3 == kernel.size()
 * @return float
 */
template <class T>
float
conv_3d(T data, const Dims dims, const VecInt voxel, std::vector<float> kernel, const uint radius)
{
    float sum = 0.0f;
    std::size_t counter = 0;
    // clip to the edge for voxels outside
    for (int iz = -radius; iz <= radius; ++iz)
    {
        int z = std::clamp<int>(voxel[2] + iz, 0, dims[2] - 1);
        for (int iy = -radius; iy <= radius; ++iy)
        {
            int y = std::clamp<int>(voxel[1] + iy, 0, dims[1] - 1);
            for (int ix = -radius; ix <= radius; ++ix)
            {
                int x = std::clamp<int>(voxel[0] + ix, 0, dims[0] - 1);
                int pointer = x + y * dims[0] + z * dims[0] * dims[1];
                sum += kernel[counter++] * static_cast<float>(data[pointer]);
            }
        }
    }
    return sum;
}

/**
 * @brief Convulution of array, the array represents the necessary neighborhood for the convolution
 *
 * @param data gradients of neighborhood of point e.g. x1,y1,z1, x2,y2,z2, ...
 * @param kernel 1d kernel
 * @return float
 */
inline
float structure_conv(std::vector<float> first, std::vector<float> second, Dims dims, VecInt startVoxel, VecInt endVoxel, std::vector<float> kernel)
{
    float sum = 0.0f;
    for (int iz = startVoxel[2]; iz < endVoxel[2]; ++iz)
    {
        float ysum = 0.0f;
        for (int iy = startVoxel[1]; iy < endVoxel[1]; ++iy)
        {
            float xsum = 0.0f;
            for (int ix = startVoxel[0]; ix < endVoxel[0]; ++ix)
            {
                const int x = ix - startVoxel[0];
                const int counter = ix + iy * dims[0] + iz * dims[0] * dims[1];
                xsum += kernel[x] * first[counter] * second[counter];
            }
            const int y = iy - startVoxel[1];
            ysum += kernel[y] * xsum;
        }
        const int z = iz - startVoxel[2];
        sum += kernel[z] * ysum;
    }
    return sum;
}

// copy data into vec (clip edges if necessary)
template <class T>
std::vector<float> clipEdges(T* data, const Dims dims, const VecInt startVoxel, const VecInt endVoxel)
{
    auto size = (endVoxel[0] - startVoxel[0]) * (endVoxel[1] - startVoxel[1]) * (endVoxel[2] - startVoxel[2]);
    std::vector<float> input = std::vector<float>(size);
    int outputCounter = 0;
    for (int iz = startVoxel[2]; iz < endVoxel[2]; ++iz)
    {
        const int z = std::max<int>(0, std::min<int>(dims[2] - 1, iz));
        for (int iy = startVoxel[1]; iy < endVoxel[1]; ++iy)
        {
            const int y = std::max<int>(0, std::min<int>(dims[1] - 1, iy));
            for (int ix = startVoxel[0]; ix < endVoxel[0]; ++ix)
            {
                const int x = std::max<int>(0, std::min<int>(dims[0] - 1, ix));
                const int counter = x + y * dims[0] + z * dims[0] * dims[1];
                input[outputCounter] = static_cast<float>(data[counter]);
                ++outputCounter;
            }
        }
    }
    return input;
}

template <class T>
Eigen::Matrix3f structure_tensor(T* data, const Dims dims, const VecInt voxel, const float gradient_sigma, const float tensor_sigma)
{
    const int gradient_radius = std::ceil(4.0f * gradient_sigma);
    const int tensor_radius = std::ceil(4.0f * tensor_sigma);
    const int radius = gradient_radius + tensor_radius;
    auto window_dims = Dims(static_cast<std::size_t>(radius * 2 + 1));
    auto startVoxel = voxel - radius;
    auto endVoxel = voxel + (radius + 1);

    auto input = clipEdges(data, dims, startVoxel, endVoxel);

    auto gradient_start = VecInt(0, 0, 0);
    auto gradient_end = VecInt(radius * 2 + 1);
    auto gradx = gauss_filter(input, window_dims, gradient_start, gradient_end, gradient_sigma, gradient_radius, VecInt::RIGHT);
    auto grady = gauss_filter(input, window_dims, gradient_start, gradient_end, gradient_sigma, gradient_radius, VecInt::UP);
    auto gradz = gauss_filter(input, window_dims, gradient_start, gradient_end, gradient_sigma, gradient_radius, VecInt::BACKWARD);

    auto tensor_start = VecInt(gradient_radius);
    auto tensor_end = gradient_end - gradient_radius;
    auto kernel = gauss_kernel(tensor_sigma, tensor_radius, 0);

    Eigen::Matrix3f tensor;
    tensor(0, 0) = structure_conv(gradx, gradx, window_dims, tensor_start, tensor_end, kernel);
    tensor(1, 1) = structure_conv(grady, grady, window_dims, tensor_start, tensor_end, kernel);
    tensor(2, 2) = structure_conv(gradz, gradz, window_dims, tensor_start, tensor_end, kernel);

    tensor(0, 1) = tensor(1, 0) = structure_conv(gradx, grady, window_dims, tensor_start, tensor_end, kernel);
    tensor(0, 2) = tensor(2, 0) = structure_conv(gradx, gradz, window_dims, tensor_start, tensor_end, kernel);
    tensor(1, 2) = tensor(2, 1) = structure_conv(grady, gradz, window_dims, tensor_start, tensor_end, kernel);

    return tensor;
}