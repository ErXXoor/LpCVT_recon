//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_REMESH_H
#define LPCVT_RECON_REMESH_H

#include <geogram/mesh/mesh.h>
#include <geogram/voronoi/integration_simplex.h>
#include "CVT/LpCVTWrap.h"
#include <geogram/mesh/mesh_AABB.h>
#include <memory>
#include <geogram/basic/smart_pointer.h>
#include "CVT/LpCVTIS.h"

namespace LpCVT {
    class Remesher {
    public:
        enum class RemeshType {
            Lloyd_CVT = 0,
            Newton_CVT,
            L2CVT,
            L8CVT
        };

        Remesher() = default;

        ~Remesher() = default;

        void Init(GEO::Mesh *M_in,
                  GEO::SmartPointer<LpCVTIS> is,
                  GEO::coord_index_t dim = 3,
                  RemeshType type = RemeshType::Lloyd_CVT);

        void Remeshing(unsigned int nb_pts = 1000,
                       unsigned int nb_iter = 100);

        void GetRVD(GEO::Mesh &M_out);

        void GetRDT(GEO::Mesh &M_out, bool post_process = false);


    private:
        std::shared_ptr<LpCVTWrap> m_cvt;
        GEO::IntegrationSimplex_var m_is = nullptr;
        std::shared_ptr<GEO::MeshFacetsAABB> m_facetsAABB = nullptr;
        RemeshType m_type;
    };
}
#endif //LPCVT_RECON_REMESH_H
