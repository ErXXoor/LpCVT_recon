//
// Created by Hongbo on 12/30/23.
//
#include "CVT/LpCVTWrap.h"
#include <geogram/basic/progress.h>

namespace LpCVT {
    LpCVTWrap::LpCVTWrap(GEO::Mesh *M, GEO::index_t dim)
            : GEO::CentroidalVoronoiTesselation(M, GEO::coord_index_t(dim)) {
//        this->is_projection_ = false;
        this->nn_search = std::make_shared<Base::NNSearch>(dim);
        this->nn_search->ConstructKDTree(*M);
    }

    void LpCVTWrap::newiteration() {
        if(progress_ != nullptr) {
            progress_->next();
        }
        cur_iter_++;

        if(cur_iter_==200){
            GEO::vector<double> nearest;
            this->nn_search->QueryPointSet(points_,nearest);

            for(auto i=0;i<points_.size();i++){
                points_[i] = 0.8*points_[i]+0.2* nearest[i];
            }
        }
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