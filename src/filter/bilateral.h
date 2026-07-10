#pragma once

// STL
#include <vector>
#include <tuple>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include <utils/Vec.hpp>
#include <algorithm>
#include <limits>

#include <surface/StaticSurface.hpp>
#include <filter/dijkstra.h>

template <class Adjacency>
struct BilateralFilterSettings {
    const std::vector<VecFloat>& input_points;
    const std::vector<Adjacency>& neighbors;
    std::vector<VecFloat>& output_points;
    float radius;
    float sigma_d;
    float sigma_n;
};

template <class T>
void set_insert(std::vector<T>& vec, T elem){
    if (std::find(vec.begin(), vec.end(), elem) == vec.end())
        vec.push_back(elem);
}

template <class T>
std::vector<std::vector<T>>
buildAdjacency(const surface::StaticSurface& surf)
{
    auto adj = std::vector<std::vector<T>>();
    const std::size_t n = surf.points().size();
    assert(n < std::numerical_limits<T>::max());
    adj.assign(n, std::vector<T>());

    const std::vector<surface::SimpleTriangle>& tris = surf.triangles();
    const std::size_t m = tris.size();

    for (std::size_t t = 0; t < m; ++t)
    {
        const surface::SimpleTriangle& tri = tris[t];
        const std::array<std::size_t, 3>& indices = tri.toRaw();
        const T v0 = indices[0];
        const T v1 = indices[1];
        const T v2 = indices[2];

        set_insert(adj[v0], v1);
        set_insert(adj[v1], v0);
        set_insert(adj[v0], v2);
        set_insert(adj[v2], v0);
        set_insert(adj[v1], v2);
        set_insert(adj[v2], v1);
    }

    return adj;
}

inline bool performLocalPCA(
    const std::vector<std::tuple<VecFloat, float>>& neighbors,
    VecFloat& barycenter,
    VecFloat& normal,
    float sigma_d)
{
    const size_t m = neighbors.size();
    if (m < 3 || sigma_d <= 0.0f)
    {
        barycenter = VecFloat(0.0f, 0.0f, 0.0f);
        normal = VecFloat(0.0f, 0.0f, 0.0f);
        return false;
    }
    // 1) weighted barycenter mu with spatial gauss weight
    // w = exp(-d^2 / (2*sigma_d^2)), d = ||q-p||
    const double inv2sig2 = 1.0 / (2.0 * double(sigma_d) * double(sigma_d));
    double sumW = 0.0;
    Eigen::Vector3d mu(0.0, 0.0, 0.0);

    for (size_t i = 0; i < m; ++i)
    {
        const VecFloat& q = std::get<0>(neighbors[i]);
        const double d = double(std::get<1>(neighbors[i])); // distance
        const double w = std::exp( -(d*d) * inv2sig2);      // gauss weight

        mu.x() += w * double(q.x());
        mu.y() += w * double(q.y());
        mu.z() += w * double(q.z());
        sumW += w;
    }

    if(sumW <= 0.0 || !std::isfinite(sumW))
    {
        barycenter = VecFloat(0.0f, 0.0f, 0.0f);
        normal = VecFloat(0.0f, 0.0f, 0.0f);
        return false;
    }

    // mean
    mu.x() /= sumW;
    mu.y() /= sumW;
    mu.z() /= sumW;

    // 2) weighted covariance C um mu
    // C = (1/sumW) *sum_i w_i * (q_i - mu)(q_i - mu)^T
    Eigen::Matrix3d C;
    C.setZero();

    for (size_t i = 0; i < m; ++i)
    {
        const VecFloat& q = std::get<0>(neighbors[i]);
        const double d = double(std::get<1>(neighbors[i]));
        const double w = std::exp( -(d*d) *inv2sig2);

        Eigen::Vector3d diff(double(q.x()), double(q.y()), double(q.z()));
        diff.x() -= mu.x();
        diff.y() -= mu.y();
        diff.z() -= mu.z();

        C += w * (diff * diff.transpose());
    }
    C /= sumW; // scaled Eigenvalue

    // 3)Eigenanalyse (smallest Eigenvector = normal)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
    if (es.info() != Eigen::Success)
    {
        // if it failes
        barycenter = VecFloat(float(mu.x()), float(mu.y()), float(mu.z()));
        normal = VecFloat(0.0f, 0.0f, 1.0f);
        return false;
    }
    // SelfAdjointEigenSolver sorted in ascending order:
    // column 0 -> Eigenvector corresponding to the smallest Eigenvalue
    Eigen::Vector3d n = es.eigenvectors().col(0);
    n.normalize();
    barycenter = VecFloat(float(mu.x()), float(mu.y()), float(mu.z()));
    normal = VecFloat(float(n.x()), float(n.y()), float(n.z()));
    return true;
}

template <class Adjacency, class Index>
void bilateral_filter_step(BilateralFilterSettings<Adjacency> settings)
{
    // Start timer
    auto startTime = std::chrono::high_resolution_clock::now();

    const auto n = settings.input_points.size();
    settings.output_points.resize(n);
    
    int64_t sumSteps = 0;

    int numThreads = 1;
#ifdef _OPENMP
    numThreads = omp_get_max_threads();
    #pragma omp parallel for
#endif
    // main filter loop
    for (int i = 0; i < n; ++i)
    {
#ifdef _OPENMP
        if (omp_get_thread_num() == 0){
            #pragma omp critical
#endif
            ++sumSteps;

            if (sumSteps % 100 == 0)
            {
                const float progress = (float)(sumSteps * numThreads) / (float)n;
                
                std::cout << "\rProgress: " << progress * 100 << "% - ";
                std::cout << sumSteps * numThreads << " / " << n << " vertices processed.";
                std::cout.flush();
            }
#ifdef _OPENMP
        }
#endif


        // Dijkstra object is not thread safe, so must be created here
        auto dijkstra = Dijkstra<Adjacency, Index>();
        const auto& neighbors = dijkstra.run(settings.input_points, settings.neighbors, i, settings.radius);

        VecFloat mu, vec;
        bool ok = performLocalPCA(neighbors, mu, vec, settings.sigma_d);
        if (!ok)
        {
            mu = settings.input_points[i];
            vec = VecFloat(0.0f, 0.0f, 1.0f);
        }
        auto nvec = vec.normalize().value_or(VecFloat(0.0f, 0.0f, 1.0f));

        const VecFloat pi = settings.input_points[i];         // current point

        double sumW = 0.0;                      // SIGMA w_ij
        double sumWH = 0.0;                     // SIGMA w_ij * h_ij (delta_p)
        for (size_t j = 0; j < neighbors.size(); ++j)
        {
            const VecFloat qj = std::get<0>(neighbors[j]);    // neighbor point
            const double d = double(std::get<1>(neighbors[j])); // || q - p || (euk. dist)

            // h_ij = <q_j - p_i, n_i> (scale product of the neighboron the normal)
            const double dx = double(qj.x()) - double(pi.x());
            const double dy = double(qj.y()) - double(pi.y());
            const double dz = double(qj.z()) - double(pi.z());
            const double h = dx * double(nvec.x()) + dy * double(nvec.y()) + dz * double(nvec.z());     // d_n

            // weights W_d(d), W_n(h)
            const double wd = std::exp( - (d * d) * (1.0 / ( 2.0 * double(settings.sigma_d) * double(settings.sigma_d))));
            const double wn = std::exp( - (h * h) * (1.0 / (2.0 * double(settings.sigma_n) * double(settings.sigma_n))));     // d_n

            const double w = wd * wn;

            sumWH += w * h;         // counter
            sumW += w;              // denominator
        }
        // new p_i = p_i + (delta * n_i), delta = sumWH / sumW
        // shift into delta, then add to it and output
        const double delta = sumWH / sumW;
        VecFloat pout = pi;
        if (sumW > 0.0 && std::isfinite(sumW))
        {
            pout[0] += float(delta * double(nvec.x()));
            pout[1] += float(delta * double(nvec.y()));
            pout[2] += float(delta * double(nvec.z()));
        }
        settings.output_points[i] = pout;
    }
    // End timer
    auto endTime = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsedTime = endTime - startTime;

    std::cout << "\nBilateral filtering finished in "
              << elapsedTime.count() << " seconds." << std::endl;
}

template <class Index = std::size_t>
surface::StaticSurface bilateral_filter_template(const surface::StaticSurface& surf, float r, float sigma_d, float sigma_n, int nIter) {
    using Adjacency = std::vector<Index>;

    std::vector<VecFloat> inPoints = surf.points();
    std::vector<VecFloat> outPoints;
    outPoints.reserve(inPoints.size());

    // settings
    auto adj = buildAdjacency<Index>(surf);
    auto settings = BilateralFilterSettings<Adjacency> {
        inPoints,
        adj,
        outPoints,
        r,
        sigma_d,
        sigma_n
    };

    for (int iter = 0; iter < nIter; ++iter) {
        std::cout << "\nIteration " << iter + 1 << " out of " << nIter << ":" << std::endl;
        bilateral_filter_step<Adjacency, Index>(settings);
        outPoints.swap(inPoints);
    }
        
    // create a raw wavefront structure
    wavefront_data_t raw;
    raw.points = inPoints;

    // convert triangles from the original surface
    raw.triangles.reserve(surf.triangles().size());
    for (const auto& t : surf.triangles())
    {
        raw.triangles.push_back(t.toRaw());
    }

    // build a new StaticSurface
    surface::StaticSurface filteredSurface =
        surface::StaticSurface::fromRaw(raw);

    return filteredSurface;
}

inline surface::StaticSurface bilateral_filter(const surface::StaticSurface& surf, float r, float sigma_d, float sigma_n, int nIter) {
    std::size_t n = surf.points().size();
    if(n <= std::numeric_limits<uint8_t>::max()){
        return bilateral_filter_template<uint8_t>(surf, r, sigma_d, sigma_n, nIter);
    }else if(n <= std::numeric_limits<uint16_t>::max()){
        return bilateral_filter_template<uint16_t>(surf, r, sigma_d, sigma_n, nIter);
    }else if(n <= std::numeric_limits<uint32_t>::max()){
        return bilateral_filter_template<uint32_t>(surf, r, sigma_d, sigma_n, nIter);
    }else if(n <= std::numeric_limits<uint64_t>::max()){
        return bilateral_filter_template<uint64_t>(surf, r, sigma_d, sigma_n, nIter);
    }else{
        throw new std::invalid_argument("Surface is too big");
    }
}