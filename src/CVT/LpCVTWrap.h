//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_LPCVTWRAP_H
#define LPCVT_RECON_LPCVTWRAP_H

#include <geogram/voronoi/CVT.h>
#include <geogram/voronoi/integration_simplex.h>

namespace LpCVT {
    class LpCVTWrap : public GEO::CentroidalVoronoiTesselation {
    public:
        LpCVTWrap(GEO::Mesh *M, GEO::index_t dim);

        void set_simplex_func(GEO::IntegrationSimplex_var is);

        void set_constrain_points(std::vector<double> constrain_points);
    };
}
#endif //LPCVT_RECON_LPCVTWRAP_H
