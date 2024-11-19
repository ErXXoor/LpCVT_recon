//
// Created by Hongbo on 12/30/23.
//
#include "CVT/LpCVTIS.h"
#include <geogram/voronoi/generic_RVD.h>
#include <geogram/basic/smart_pointer.h>
#include <geogram/basic/geometry_nd.h>
#include <stan/math/rev/core/chainablestack.hpp>
#include <igl/rotation_matrix_from_directions.h>

namespace LpCVT {
    void point2vec(const double *point,
                   unsigned int dim,
                   Eigen::VectorXd &vec) {
        vec = Eigen::VectorXd::Zero(dim);
        for (unsigned int i = 0; i < dim; i++) {
            vec(i) = point[i];
        }
    }

    void point2vec(const double *point,
                   unsigned int dim,
                   Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> &vec) {
        vec = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>::Zero(dim);
        for (unsigned int i = 0; i < dim; i++) {
            vec(i) = point[i];
        }
    }

    void soft_zero(Eigen::Vector3d &vec, double eps = 1e-10) {
        for (auto i = 0; i < vec.size(); i++) {
            if (::fabs(vec[i]) < eps) {
                vec[i] = 0.0;
            }
        }
    }

    LpCVTIS::LpCVTIS(const GEO::Mesh &mesh,
                     bool volumetric,
                     unsigned int dim,
                     unsigned int degree,
                     double metric_weight) :
            IntegrationSimplex(mesh,
                               volumetric,
                               0,
                               0,
                               nullptr),
            m_dim(dim) {
        m_degree = degree;
        nb_coeffs = ((m_degree + 1) * (m_degree + 2)) / 2;
        nb_dcoeffs = nb_coeffs - (m_degree + 1);
        m_metric_weight = metric_weight;

        // Pre-compute indices for evaluating F_{L_p}^T
        // (see Appendix A)
        {
            for (unsigned int alpha = 0; alpha <= m_degree; alpha++) {
                for (unsigned int beta = 0; beta <= m_degree - alpha; beta++) {
                    unsigned int gamma = m_degree - alpha - beta;
                    std::vector<unsigned int> row_tmp;
                    row_tmp.push_back(alpha);
                    row_tmp.push_back(beta);
                    row_tmp.push_back(gamma);
                    E_pow.emplace_back(row_tmp);
                }
            }
        }

        // Pre-compute indices for evaluating \nabla F_{L_p}^T
        // (see Appendix B.1)
        {
            dE_pow.resize(3);

            for (unsigned int alpha = 0; alpha <= m_degree; alpha++) {
                for (unsigned int beta = 0; beta <= m_degree - alpha; beta++) {
                    unsigned int gamma = m_degree - alpha - beta;
                    if (alpha != 0) {
                        std::vector<unsigned int> dU1_tmp;
                        dU1_tmp.push_back(alpha);
                        dU1_tmp.push_back(beta);
                        dU1_tmp.push_back(gamma);
                        dE_pow[0].emplace_back(dU1_tmp);
                    }
                    if (beta != 0) {
                        std::vector<unsigned int> dU2_tmp;
                        dU2_tmp.push_back(alpha);
                        dU2_tmp.push_back(beta);
                        dU2_tmp.push_back(gamma);
                        dE_pow[1].emplace_back(dU2_tmp);
                    }
                    if (gamma != 0) {
                        std::vector<unsigned int> dU3_tmp;
                        dU3_tmp.push_back(alpha);
                        dU3_tmp.push_back(beta);
                        dU3_tmp.push_back(gamma);
                        dE_pow[2].emplace_back(dU3_tmp);
                    }
                }
            }
        }

    }

    const double *LpCVTIS::vertex_ptr(GEO::index_t i) const {
        return points_ + i * points_stride_;
    }

    void LpCVTIS::set_facetsAABB(std::shared_ptr<GEO::MeshFacetsAABB> facetsAABB) {
        if (facetsAABB == nullptr) {
            std::cerr << "facetsAABB is nullptr" << std::endl;
            return;
        }
        m_facetsAABB = facetsAABB;

        for (auto i = 0; i < m_facetsAABB->mesh()->facets.nb(); i++) {
            unsigned int a = m_facetsAABB->mesh()->facets.vertex(i, 0);
            unsigned int b = m_facetsAABB->mesh()->facets.vertex(i, 1);
            unsigned int c = m_facetsAABB->mesh()->facets.vertex(i, 2);
            m_vertices_id.push_back({a, b, c});
            std::vector<double> vert_a = std::vector<double>(m_facetsAABB->mesh()->vertices.point_ptr(a),
                                                             m_facetsAABB->mesh()->vertices.point_ptr(a) + 3);
            std::vector<double> vert_b = std::vector<double>(m_facetsAABB->mesh()->vertices.point_ptr(b),
                                                             m_facetsAABB->mesh()->vertices.point_ptr(b) + 3);
            std::vector<double> vert_c = std::vector<double>(m_facetsAABB->mesh()->vertices.point_ptr(c),
                                                             m_facetsAABB->mesh()->vertices.point_ptr(c) + 3);
            m_vertices.push_back(vert_a);
            m_vertices.push_back(vert_b);
            m_vertices.push_back(vert_c);
        }
    }

    void LpCVTIS::set_mesh(std::shared_ptr<Mesh> mesh) {
        m_mesh = mesh;
    }

    void LpCVTIS::CalQuadMetric() {
        auto fn = m_mesh->GetFaceNormals();
        auto cross_w = m_mesh->GetCrossFieldW();
        auto cross_v = m_mesh->GetCrossFieldV();
        for (auto i = 0; i < fn.rows(); i++) {
            Eigen::Matrix3d R;
            R.row(0) = cross_w.row(i);
            R.row(1) = cross_v.row(i);
            R.row(2) = fn.row(i);
            Eigen::Matrix3d S = Eigen::Matrix3d::Identity();
            S(2, 2) *= m_metric_weight;
            Eigen::Matrix3d M = R * S * R.transpose();
            m_quad_metric.emplace_back(M);
        }
    }

    double LpCVTIS::eval(
            GEO::index_t v,
            const GEOGen::Vertex &p1,
            const GEOGen::Vertex &p2,
            const GEOGen::Vertex &p3,
            GEO::index_t t,
            GEO::index_t t_adj,
            GEO::index_t v_adj
    ) {


        // Compute normal of the triangle
        auto p0 = vertex_ptr(v);
        Eigen::VectorXd p1_vec, p2_vec, p3_vec, p0_vec;
        p1_vec.resize(m_dim);
        p2_vec.resize(m_dim);
        p3_vec.resize(m_dim);
        p0_vec.resize(m_dim);
        for (auto i = 0; i < m_dim; i++) {
            p1_vec[i] = p1.point()[i];
            p2_vec[i] = p2.point()[i];
            p3_vec[i] = p3.point()[i];
            p0_vec[i] = p0[i];
        }

        auto p1_vec3 = Eigen::Vector3d(p1.point()[0], p1.point()[1], p1.point()[2]);
        auto p2_vec3 = Eigen::Vector3d(p2.point()[0], p2.point()[1], p2.point()[2]);
        auto p3_vec3 = Eigen::Vector3d(p3.point()[0], p3.point()[1], p3.point()[2]);

        Eigen::MatrixXd M = Eigen::MatrixXd::Identity(m_dim, m_dim);
        Eigen::Vector3d N = (p2_vec3 - p1_vec3).cross(p3_vec3 - p1_vec3);
        N.normalize();

        if (m_metric_type == MetricType::Quad) {
//            auto mid_tmp = (p1_vec3 + p2_vec3 + p3_vec3) / 3.0;
//            auto mid = GEO::vec3(mid_tmp[0], mid_tmp[1], mid_tmp[2]);
//            auto f_id = m_facetsAABB->nearest_facet(mid);
//            auto metric = m_quad_metric[t];
//            Eigen::Vector3d normal = m_quad_metric[f_id].col(2);
//            metric.col(2) = N;

//            M.block<3, 3>(0, 0) = m_quad_metric[t];
        } else {
            Eigen::Matrix3d N_mat3 = N * N.transpose();
            N_mat3 = N_mat3 * m_metric_weight;

            M.block<3, 3>(0, 0) += N_mat3;
        }

        Eigen::VectorXd dFdp1, dFdp2, dFdp3;

//        double f = eval_explicit(p0_vec, p1_vec, p2_vec, p3_vec, M, dFdp1, dFdp2, dFdp3);
        stan::math::ChainableStack stack_;
        double f = eval_ad(p0_vec, p1_vec, p2_vec, p3_vec, M, dFdp1, dFdp2, dFdp3);

//        grad_site(v, p1, N, dFdp1);
//        grad_site(v, p2, N, dFdp2);
//        grad_site(v, p3, N, dFdp3);

        for (auto i = 0; i < m_dim; i++) {
            g_[m_dim * v + i] += -dFdp1[i] - dFdp2[i] - dFdp3[i];
        }

        return f;
    }

    double LpCVTIS::eval_explicit(const Eigen::VectorXd &p0_vec,
                                  const Eigen::VectorXd &p1_vec,
                                  const Eigen::VectorXd &p2_vec,
                                  const Eigen::VectorXd &p3_vec,
                                  const Eigen::MatrixXd &M,
                                  Eigen::VectorXd &dFdp1,
                                  Eigen::VectorXd &dFdp2,
                                  Eigen::VectorXd &dFdp3) {
        std::vector<std::vector<Eigen::VectorXd>> U_pow;
        U_pow.resize(3);
        U_pow[0].resize(m_degree + 1);
        U_pow[1].resize(m_degree + 1);
        U_pow[2].resize(m_degree + 1);
        U_pow[0][0] = Eigen::VectorXd::Ones(m_dim);
        U_pow[1][0] = Eigen::VectorXd::Ones(m_dim);
        U_pow[2][0] = Eigen::VectorXd::Ones(m_dim);

        Eigen::VectorXd u0 = p1_vec - p0_vec;
        Eigen::VectorXd u1 = p2_vec - p0_vec;
        Eigen::VectorXd u2 = p3_vec - p0_vec;
        U_pow[0][1] = M * u0;
        U_pow[1][1] = M * u1;
        U_pow[2][1] = M * u2;

        for (unsigned int i = 2; i <= m_degree; i++) {
            U_pow[0][i] = U_pow[0][1].cwiseProduct(U_pow[0][i - 1]);
            U_pow[1][i] = U_pow[1][1].cwiseProduct(U_pow[1][i - 1]);
            U_pow[2][i] = U_pow[2][1].cwiseProduct(U_pow[2][i - 1]);
        }

        // Computation of function value.
        double E = 0.0;
        for (unsigned int i = 0; i < nb_coeffs; i++) {
            Eigen::VectorXd W;
            unsigned int alpha = E_pow[i][0];
            unsigned int beta = E_pow[i][1];
            unsigned int gamma = E_pow[i][2];
            W = U_pow[0][alpha].cwiseProduct(U_pow[1][beta]).cwiseProduct(U_pow[2][gamma]);
            E += W.sum();
        }

        // Computation of gradient
        Eigen::VectorXd dEdU1, dEdU2, dEdU3;
        dEdU1.setZero(m_dim);
        dEdU2.setZero(m_dim);
        dEdU3.setZero(m_dim);
        for (unsigned int i = 0; i < nb_dcoeffs; i++) {
            {
                unsigned int alpha = dE_pow[0][i][0];
                unsigned int beta = dE_pow[0][i][1];
                unsigned int gamma = dE_pow[0][i][2];
                dEdU1 += alpha * U_pow[0][alpha - 1].cwiseProduct(U_pow[1][beta]).cwiseProduct(U_pow[2][gamma]);
            }
            {
                unsigned int alpha = dE_pow[1][i][0];
                unsigned int beta = dE_pow[1][i][1];
                unsigned int gamma = dE_pow[1][i][2];
                dEdU2 += beta * U_pow[0][alpha].cwiseProduct(U_pow[1][beta - 1]).cwiseProduct(U_pow[2][gamma]);

            }
            {
                unsigned int alpha = dE_pow[2][i][0];
                unsigned int beta = dE_pow[2][i][1];
                unsigned int gamma = dE_pow[2][i][2];
                dEdU3 += gamma * U_pow[0][alpha].cwiseProduct(U_pow[1][beta]).cwiseProduct(U_pow[2][gamma - 1]);
            }
        }

        // Compute tet measure and its
        // derivatives relative to U1, U2 and U3.

        std::vector<double> dTdU1_vec, dTdU2_vec, dTdU3_vec;
        dTdU1_vec.resize(m_dim);
        dTdU2_vec.resize(m_dim);
        dTdU3_vec.resize(m_dim);

        double T = grad_tri(m_dim,
                            U_pow[0][1].data(),
                            U_pow[1][1].data(),
                            U_pow[2][1].data(),
                            dTdU1_vec,
                            dTdU2_vec,
                            dTdU3_vec);
        Eigen::Map<Eigen::VectorXd> dTdU1(dTdU1_vec.data(), dTdU1_vec.size());
        Eigen::Map<Eigen::VectorXd> dTdU2(dTdU2_vec.data(), dTdU2_vec.size());
        Eigen::Map<Eigen::VectorXd> dTdU3(dTdU3_vec.data(), dTdU3_vec.size());

        // Assemble dF = E.d|T| + |T|.dE
        // Rem: anisotropy matrix needs to be transposed
        // grad(F(MX)) = J(F(MX))^T = (gradF^T M)^T = M^T grad F
        Eigen::VectorXd dFdU1, dFdU2, dFdU3;
        dFdp1 = M.transpose() * (E * dTdU1 + T * dEdU1);
        dFdp2 = M.transpose() * (E * dTdU2 + T * dEdU2);
        dFdp3 = M.transpose() * (E * dTdU3 + T * dEdU3);

        return T * E;
    }

    double LpCVTIS::eval_ad(const Eigen::VectorXd &p0,
                            const Eigen::VectorXd &p1,
                            const Eigen::VectorXd &p2,
                            const Eigen::VectorXd &p3,
                            const Eigen::MatrixXd &M,
                            Eigen::VectorXd &dFdp1,
                            Eigen::VectorXd &dFdp2,
                            Eigen::VectorXd &dFdp3) {
        Eigen::VectorXd U1 = M * (p1 - p0);
        Eigen::VectorXd U2 = M * (p2 - p0);
        Eigen::VectorXd U3 = M * (p3 - p0);

        // Computation of function value.
        std::vector<std::vector<Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>>> U_pow;
        U_pow.resize(3);
        U_pow[0].resize(m_degree + 1);
        U_pow[1].resize(m_degree + 1);
        U_pow[2].resize(m_degree + 1);
        U_pow[0][0] = Eigen::VectorXd::Ones(m_dim);
        U_pow[1][0] = Eigen::VectorXd::Ones(m_dim);
        U_pow[2][0] = Eigen::VectorXd::Ones(m_dim);

        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U1_E = U1;
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U2_E = U2;
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U3_E = U3;
        U_pow[0][1] = U1_E;
        U_pow[1][1] = U2_E;
        U_pow[2][1] = U3_E;

        for (unsigned int i = 2; i <= m_degree; i++) {
            U_pow[0][i] = U_pow[0][1].cwiseProduct(U_pow[0][i - 1]);
            U_pow[1][i] = U_pow[1][1].cwiseProduct(U_pow[1][i - 1]);
            U_pow[2][i] = U_pow[2][1].cwiseProduct(U_pow[2][i - 1]);
        }

        stan::math::var E = 0.0;
        for (unsigned int i = 0; i < nb_coeffs; i++) {
            Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> W;
            unsigned int alpha = E_pow[i][0];
            unsigned int beta = E_pow[i][1];
            unsigned int gamma = E_pow[i][2];
            W = U_pow[0][alpha].cwiseProduct(U_pow[1][beta]).cwiseProduct(U_pow[2][gamma]);
            E += W.sum();
        }

        Eigen::VectorXd dEdU1, dEdU2, dEdU3;
        dEdU1.setZero(m_dim);
        dEdU2.setZero(m_dim);
        dEdU3.setZero(m_dim);
        E.grad();
        dEdU1 = U1_E.adj();
        dEdU2 = U2_E.adj();
        dEdU3 = U3_E.adj();


        // Computation of area
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U1_T = U1;
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U2_T = U2;
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U3_T = U3;

        stan::math::var a = stan::math::distance(U1_T, U2_T);
        stan::math::var b = stan::math::distance(U2_T, U3_T);
        stan::math::var c = stan::math::distance(U3_T, U1_T);
        stan::math::var s = double(0.5) * (a + b + c);
        stan::math::var A2 = s * (s - a) * (s - b) * (s - c);
        A2 = stan::math::if_else(A2 < 0.0, 0.0, A2);
        stan::math::var T = stan::math::sqrt(A2);

        Eigen::VectorXd dTdU1, dTdU2, dTdU3;
        dTdU1.setZero(m_dim);
        dTdU2.setZero(m_dim);
        dTdU3.setZero(m_dim);

        T.grad();
        if (A2 >= 1e-10) {
            dTdU1 = U1_T.adj();
            dTdU2 = U2_T.adj();
            dTdU3 = U3_T.adj();
        }
        dFdp1 = M.transpose() * (E.val() * dTdU1 + T.val() * dEdU1);
        dFdp2 = M.transpose() * (E.val() * dTdU2 + T.val() * dEdU2);
        dFdp3 = M.transpose() * (E.val() * dTdU3 + T.val() * dEdU3);

        return T.val() * E.val();
    }

    double
    LpCVTIS::grad_tri(const GEO::vec3 &U1,
                      const GEO::vec3 &U2,
                      const GEO::vec3 &U3,
                      GEO::vec3 &dTdU1,
                      GEO::vec3 &dTdU2,
                      GEO::vec3 &dTdU3) {
        GEO::vec3 N = cross(U1 - U3, U2 - U3);
        double T = length(N);
        if (::fabs(T) < 1e-10) {
            dTdU1 = GEO::vec3(0.0, 0.0, 0.0);
            dTdU2 = GEO::vec3(0.0, 0.0, 0.0);
            dTdU3 = GEO::vec3(0.0, 0.0, 0.0);
        } else {
            N = (1.0 / T) * N;
            dTdU1 = cross(N, U3 - U2);
            dTdU2 = cross(N, U1 - U3);
            dTdU3 = cross(N, U2 - U1);
        }
        return T;
    }

    double LpCVTIS::grad_tri(unsigned int dim,
                             const double *U1,
                             const double *U2,
                             const double *U3,
                             std::vector<double> &dTdU1,
                             std::vector<double> &dTdU2,
                             std::vector<double> &dTdU3) {
        stan::math::ChainableStack stack_;
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U1_vec(dim);
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U2_vec(dim);
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> U3_vec(dim);
        U1_vec.setZero();
        U2_vec.setZero();
        U3_vec.setZero();

        for (auto i = 0; i < dim; i++) {
            U1_vec(i) = U1[i];
            U2_vec(i) = U2[i];
            U3_vec(i) = U3[i];
        }
        stan::math::var a = stan::math::distance(U1_vec, U2_vec);
        stan::math::var b = stan::math::distance(U2_vec, U3_vec);
        stan::math::var c = stan::math::distance(U3_vec, U1_vec);

        stan::math::var s = double(0.5) * (a + b + c);
        stan::math::var A2 = s * (s - a) * (s - b) * (s - c);

        if (A2 < 1e-10) {
            for (auto i = 0; i < dim; i++) {
                dTdU1[i] = 0.0;
                dTdU2[i] = 0.0;
                dTdU3[i] = 0.0;
            }
            A2 = stan::math::if_else(A2 < 0.0, 0.0, A2);
            return stan::math::sqrt(A2).val();
        }

        stan::math::var A = stan::math::sqrt(A2);
        Eigen::Matrix<double, 1, Eigen::Dynamic> dAdU1(dim);
        Eigen::Matrix<double, 1, Eigen::Dynamic> dAdU2(dim);
        Eigen::Matrix<double, 1, Eigen::Dynamic> dAdU3(dim);

        A.grad();
        dAdU1 = U1_vec.adj();
        dAdU2 = U2_vec.adj();
        dAdU3 = U3_vec.adj();

        for (auto i = 0; i < dim; i++) {
            dTdU1[i] = dAdU1(i);
            dTdU2[i] = dAdU2(i);
            dTdU3[i] = dAdU3(i);
        }

        return A.val();
    }

    void LpCVTIS::add_grad(Eigen::VectorXd grad, GEO::index_t site_id) {
        for (auto i = 0; i < m_dim; i++) {
            g_[m_dim * site_id + i] += grad[i];
        }
    }

    void LpCVTIS::compute_W(const Eigen::Vector3d &N0,
                            const Eigen::Vector3d &N1,
                            const Eigen::Vector3d &N2,
                            Eigen::Vector3d &W0,
                            Eigen::Vector3d &W1) {
        W0 = N1.cross(N2);
        double Tinv = 1.0;
        if (N0.dot(W0) != 0)
            Tinv = 1.0 / (N0.dot(W0));
        W0 *= Tinv;

        W1 = N2.cross(N0);
        W1 *= Tinv;
    }

    void intersect_geom(Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> &site_a,
                        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> &site_b,
                        Eigen::VectorXd &edg_v1,
                        Eigen::VectorXd &edg_v2,
                        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> &C
    ) {
        Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> n = site_a - site_b;
        stan::math::var d = -n.dot(site_a + site_b);
        stan::math::var l1 = edg_v2.dot(n);
        stan::math::var l2 = edg_v1.dot(n);
        d = 0.5 * d;
        l1 = stan::math::fabs(l1 + d);
        l2 = stan::math::fabs(l2 + d);
        stan::math::var l = l1 + l2;
        l1 = stan::math::if_else(l > 1e-30, l1 / l, 0.5);
        l2 = stan::math::if_else(l > 1e-30, l2 / l, 0.5);
        C = l1 * edg_v1 + l2 * edg_v2;
    }

    void LpCVTIS::grad_site(GEO::index_t site_id,
                            const GEOGen::Vertex &C,
                            const Eigen::Vector3d &N,
                            const Eigen::VectorXd &dFdC) {


        Eigen::VectorXd site_vec;
        point2vec(vertex_ptr(site_id), m_dim, site_vec);

        Eigen::VectorXd C_vec;
        point2vec(C.point(), m_dim, C_vec);

        switch (C.sym().nb_bisectors()) {
            case 1: {
                // config B: The point v is the intersection between
                //   two facets of the surface and one bisector.

                Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> site_ad = site_vec;
                auto b = C.sym().bisector(0);
                auto site_b = m_dly->vertex_ptr(b);
                Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> site_b_ad;
                point2vec(site_b, m_dim, site_b_ad);

                GEO::index_t edg_v1_idx, edg_v2_idx;
                C.sym().get_boundary_edge(edg_v1_idx, edg_v2_idx);

                Eigen::VectorXd edg_v1;
                point2vec(
                        mesh_.vertices.point_ptr(edg_v1_idx), m_dim, edg_v1);
                Eigen::VectorXd edg_v2;
                point2vec(
                        mesh_.vertices.point_ptr(edg_v2_idx), m_dim, edg_v2);

                Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> C_ad;
                intersect_geom(site_ad, site_b_ad, edg_v1, edg_v2, C_ad);

                //site gradient
                Eigen::MatrixXd J_site, J_site_b;
                J_site.resize(C_ad.rows(), site_ad.rows());
                J_site_b.resize(C_ad.rows(), site_b_ad.rows());
                for (auto i = 0; i < C_ad.rows(); i++) {
                    if (i > 0) {
                        stan::math::set_zero_all_adjoints();
                    }
                    C_ad(i).grad();
                    for (auto j = 0; j < site_ad.rows(); j++) {
                        J_site(i, j) = site_ad(j).adj();
                        J_site_b(i, j) = site_b_ad(j).adj();
                    }
                }

                Eigen::VectorXd dFdCJ = dFdC.transpose() * J_site;
                add_grad(dFdCJ, site_id);

                Eigen::VectorXd dFdCJ_b = dFdC.transpose() * J_site_b;
                add_grad(dFdCJ_b, b);

                break;
            }
            case 2: {
                // config C: The point v is the intersection between
                // one facets of the surface and two bisectors.

                //site gradient for a, c
                auto b = C.sym().bisector(0);
                auto site_b = m_dly->vertex_ptr(b);
                Eigen::VectorXd site_b_vec;
                point2vec(site_b, m_dim, site_b_vec);

                auto c = C.sym().bisector(1);
                auto site_c = m_dly->vertex_ptr(c);
                Eigen::VectorXd site_c_vec;
                point2vec(site_c, m_dim, site_c_vec);


                Eigen::Vector3d W0, W1;
                compute_W((site_b_vec - site_vec).head<3>(),
                          (site_c_vec - site_vec).head<3>(),
                          N, W0, W1);

                Eigen::Matrix3d J;
                //============= dC/dp0
                auto d0 = C_vec - site_vec;
                auto W = W0 + W1;
                J = d0 * W.transpose();
                Eigen::Vector3d dFdCJ_ = J * dFdC;
                add_grad(dFdCJ_, site_id);

                //============= dC/dp1
                auto d1 = site_b_vec - C_vec;
                J = d1 * W0.transpose();
                dFdCJ_ = J * dFdC;
                add_grad(dFdCJ_, b);

                //============= dC/dp2
                auto d2 = site_c_vec - C_vec;
                J = d2 * W1.transpose();
                dFdCJ_ = J * dFdC;
                add_grad(dFdCJ_, c);
                break;
            }
            default:
                break;
        }

    }
}