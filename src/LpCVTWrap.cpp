//
// Created by Hongbo on 12/30/23.
//
#include "CVT/LpCVTWrap.h"

namespace LpCVT {
    LpCVTWrap::LpCVTWrap(GEO::Mesh *M, GEO::index_t dim)
            : GEO::CentroidalVoronoiTesselation(M, GEO::coord_index_t(dim)) {
    }

    void LpCVTWrap::set_simplex_func(GEO::IntegrationSimplex_var is) {
        simplex_func_ = is;
    }
}