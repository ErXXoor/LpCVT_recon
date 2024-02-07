//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_LPCVTIS_H
#define LPCVT_RECON_LPCVTIS_H

#include <stan/math.hpp>
#include <geogram/voronoi/integration_simplex.h>
#include <geogram/mesh/mesh_AABB.h>
#include "Base/Mesh.h"
#include <memory>

namespace LpCVT {
    class LpCVTIS : public GEO::IntegrationSimplex {
    public:
        LpCVTIS(const GEO::Mesh &mesh, bool volumetric, unsigned int dim, unsigned int degree,
                double metric_weight = 6.0);

        ~LpCVTIS() override = default;

        const double *vertex_ptr(GEO::index_t v) const;

        void CalQuadMetric();

        void set_facetsAABB(std::shared_ptr<GEO::MeshFacetsAABB> facetsAABB);

        void set_mesh(std::shared_ptr<Mesh> mesh);

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

        //Don't want to expose Eigen
        double grad_tri(unsigned int dim,
                        const double *U1,
                        const double *U2,
                        const double *U3,
                        std::vector<double> &dTdU1,
                        std::vector<double> &dTdU2,
                        std::vector<double> &dTdU3);

        //Utils
        void vecmul(const double *p1, const double *p2, double *to);

        void vecmul(const double *p1, const double *p2, const double *p3, double *to);

        double vecbar(const double *p1);

        void vecmadd(double s, const double *p1, const double *p2, const double *p3, double *to);

        void vecmadd(double s, const GEO::vec3 &p1, double t, const GEO::vec3 &p2, GEO::vec3 &to);

        void matTvecmul(const GEO::mat3 &M, const GEO::vec3 &U, GEO::vec3 &V);


    private:
        const unsigned int m_dim;
        unsigned int m_degree;
        unsigned int nb_coeffs;
        unsigned int nb_dcoeffs;
        double m_metric_weight;
        std::vector<std::vector<unsigned int>> E_pow;
        std::vector<std::vector<std::vector<unsigned int>>> dE_pow;
        std::shared_ptr<GEO::MeshFacetsAABB> m_facetsAABB = nullptr;
        std::shared_ptr<Mesh> m_mesh;
        std::vector<std::vector<unsigned int>> m_vertices_id;
        std::vector<std::vector<double>> m_vertices;
        std::vector<Eigen::Matrix3d> m_quad_metric;

    };

}
#endif //LPCVT_RECON_LPCVTIS_H
