//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_LPCVTIS_H
#define LPCVT_RECON_LPCVTIS_H

#include <geogram/voronoi/integration_simplex.h>
#include <geogram/mesh/mesh_AABB.h>

namespace LpCVT {
    class LpCVTIS : public GEO::IntegrationSimplex {
    public:
        LpCVTIS(const GEO::Mesh &mesh, bool volumetric, unsigned int dim, unsigned int degree);

        ~LpCVTIS() override = default;

        const double *vertex_ptr(GEO::index_t v) const;

        void set_facetsAABB(std::shared_ptr<GEO::MeshFacetsAABB> facetsAABB);

        double eval(
                GEO::index_t center_vertex_index,

                const GEOGen::Vertex &v0,
                const GEOGen::Vertex &v1,
                const GEOGen::Vertex &v2,
                GEO::index_t t,
                GEO::index_t t_adj = GEO::index_t(-1),
                GEO::index_t v_adj = GEO::index_t(-1)
        ) override;

        double grad_tri(const GEO::vec3 &U1, const GEO::vec3 &U2, const GEO::vec3 &U3,
                        GEO::vec3 &dTdU1, GEO::vec3 &dTdU2, GEO::vec3 &dTdU3);

        //Utils
        void vecmul(const double *p1, const double *p2, double *to);

        void vecmul(const double *p1, const double *p2, const double *p3, double *to);

        double vecbar(const double *p1);

        void vecmadd(double s, const double *p1, const double *p2, const double *p3, double *to);

        void vecmadd(double s, const GEO::vec3 &p1, double t, const GEO::vec3 &p2, GEO::vec3 &to);

        void matTvecmul(const GEO::mat3 &M, const GEO::vec3 &U, GEO::vec3 &V);

    private:
        unsigned int m_dim;
        unsigned int m_degree;
        unsigned int nb_coeffs;
        unsigned int nb_dcoeffs;
        std::vector<std::vector<unsigned int>> E_pow;
        std::vector<std::vector<std::vector<unsigned int>>> dE_pow;
        std::shared_ptr<GEO::MeshFacetsAABB> m_facetsAABB = nullptr;
        GEO::Attribute<double> m_face_normal;

    };

}
#endif //LPCVT_RECON_LPCVTIS_H
