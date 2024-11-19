//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_LPCVTWRAP_H
#define LPCVT_RECON_LPCVTWRAP_H

#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/integration_simplex.h>
#include "Base/NNSearch.h"
#include <memory>
namespace LpCVT {
    class LpCVTWrap : public GEO::CentroidalVoronoiTesselation {
    public:
        LpCVTWrap(GEO::Mesh *M, GEO::index_t dim);

        void set_simplex_func(GEO::IntegrationSimplex_var is);
        void newiteration() override;

        void set_constrain_points(std::vector<double> constrain_points);
    private:
        std::shared_ptr<Base::NNSearch> nn_search;
    };
}
#endif //LPCVT_RECON_LPCVTWRAP_H
