//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_LPCVTIS_H
#define LPCVT_RECON_LPCVTIS_H

#include <stan/math.hpp>
#include <geogram/voronoi/integration_simplex.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/delaunay/delaunay.h>
#include "Base/Mesh.h"
#include <memory>

namespace LpCVT {
    class LpCVTIS : public GEO::IntegrationSimplex {
    public:
        enum class MetricType {
            Feature = 0,
            Quad
        };

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

        double eval_explicit(const Eigen::VectorXd &p0,
                             const Eigen::VectorXd &p1,
                             const Eigen::VectorXd &p2,
                             const Eigen::VectorXd &p3,
                             const Eigen::MatrixXd &M,
                             Eigen::VectorXd &dFdp1,
                             Eigen::VectorXd &dFdp2,
                             Eigen::VectorXd &dFdp3);

        double eval_ad(const Eigen::VectorXd &p0,
                       const Eigen::VectorXd &p1,
                       const Eigen::VectorXd &p2,
                       const Eigen::VectorXd &p3,
                       const Eigen::MatrixXd &M,
                       Eigen::VectorXd &dFdp1,
                       Eigen::VectorXd &dFdp2,
                       Eigen::VectorXd &dFdp3);

        void grad_site(GEO::index_t site_id,
                       const GEOGen::Vertex &v,
                       const Eigen::Vector3d &N,
                       const Eigen::VectorXd &dFdC,
                       const GEO::index_t t);

        void compute_W(const Eigen::Vector3d &N0,
                       const Eigen::Vector3d &N1,
                       const Eigen::Vector3d &N2,
                       Eigen::Vector3d &W0,
                       Eigen::Vector3d &W1);


        double grad_tri(unsigned int dim,
                        const double *U1,
                        const double *U2,
                        const double *U3,
                        std::vector<double> &dTdU1,
                        std::vector<double> &dTdU2,
                        std::vector<double> &dTdU3);

        //Utils
        void add_grad(Eigen::VectorXd grad, GEO::index_t site_id);

        double grad_tri(const GEO::vec3 &U1, const GEO::vec3 &U2, const GEO::vec3 &U3,
                        GEO::vec3 &dTdU1, GEO::vec3 &dTdU2, GEO::vec3 &dTdU3);


    private:
        const unsigned int m_dim;
        unsigned int m_degree;
        unsigned int nb_coeffs;
        unsigned int nb_dcoeffs;
        double m_metric_weight;
        MetricType m_metric_type;
        std::vector<std::vector<unsigned int>> E_pow;
        std::vector<std::vector<std::vector<unsigned int>>> dE_pow;
        std::shared_ptr<GEO::MeshFacetsAABB> m_facetsAABB = nullptr;
        std::shared_ptr<Mesh> m_mesh;
        GEO::Mesh *m_mesh_geo;
        std::vector<std::vector<unsigned int>> m_vertices_id;
        std::vector<std::vector<double>> m_vertices;
        std::vector<Eigen::Matrix3d> m_quad_metric;

        GEO::Delaunay *m_dly;

    public:
        void set_Delaunay(GEO::Delaunay *delaunay) {
            m_dly = delaunay;
        }

        void set_mesh_geo(GEO::Mesh *mesh_geo) {
            m_mesh_geo = mesh_geo;
        }

        void set_metric_type(MetricType metric_type) {
            m_metric_type = metric_type;
        }

    };

}
#endif //LPCVT_RECON_LPCVTIS_H
