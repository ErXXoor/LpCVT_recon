//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_REMESH_H
#define LPCVT_RECON_REMESH_H

#include <geogram/mesh/mesh.h>
#include <geogram/voronoi/integration_simplex.h>
#include "CVT/LpCVTWrap.h"

namespace LpCVT {
    class Remesher {
    public:
        enum class RemeshType {
            Lloyd_CVT = 0,
            Newton_CVT,
            LPCVT
        };

        Remesher() = default;

        ~Remesher() = default;

        void Init(const GEO::Mesh &M_in,
                  GEO::coord_index_t dim = 3,
                  RemeshType type = RemeshType::Lloyd_CVT);

        void Remeshing(unsigned int nb_pts = 1000, unsigned int nb_iter = 100);

        void GetRVD(GEO::Mesh &M_out);

        void GetRDT(GEO::Mesh &M_out);


    private:
        std::shared_ptr<LpCVTWrap> m_cvt;
        GEO::IntegrationSimplex_var m_is = nullptr;
        RemeshType m_type;
    };
}
#endif //LPCVT_RECON_REMESH_H
