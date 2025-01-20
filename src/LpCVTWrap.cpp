//
// Created by Hongbo on 12/30/23.
//
#include "CVT/LpCVTWrap.h"
#include <geogram/basic/progress.h>

namespace LpCVT {
    LpCVTWrap::LpCVTWrap(GEO::Mesh *M, GEO::index_t dim)
            : GEO::CentroidalVoronoiTesselation(M, GEO::coord_index_t(dim)) {
    }

    void LpCVTWrap::set_simplex_func(GEO::IntegrationSimplex_var is) {
        simplex_func_ = is;
    }

    void LpCVTWrap::set_constrain_points(std::vector<double> constrain_points) {
        auto nb_samples = points_.size() / dimension();
        auto constrain_size = constrain_points.size()/dimension();
        std::vector<int> constrain_id;

        for(auto i = 0;i<points_.size();i++){
            constrain_points.push_back(points_[i]);
        }

        set_points(constrain_size+nb_samples, constrain_points.data());

        for(auto i=0; i<constrain_size; i++){
            lock_point(i);
        }

    }
}