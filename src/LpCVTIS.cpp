//
// Created by Hongbo on 12/30/23.
//
#include "CVT/LpCVTIS.h"
#include <geogram/voronoi/generic_RVD.h>
#include <geogram/basic/smart_pointer.h>
#include <geogram/basic/geometry_nd.h>
#include <stan/math.hpp>
#include <thread>
#include <stan/math/rev/core/chainablestack.hpp>

namespace LpCVT {


    LpCVTIS::LpCVTIS(const GEO::Mesh &mesh,
                     bool volumetric,
                     unsigned int dim,
                     unsigned int degree) :
            IntegrationSimplex(mesh,
                               volumetric,
                               0,
                               0,
                               nullptr),
            m_dim(dim) {
        m_degree = degree;
        nb_coeffs = ((m_degree + 1) * (m_degree + 2)) / 2;
        nb_dcoeffs = nb_coeffs - (m_degree + 1);
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
        m_face_normal.bind_if_is_defined(facetsAABB->mesh()->facets.attributes(), "normal");
        if (!m_face_normal.is_bound() || m_face_normal.dimension() != 3) {
            std::cerr << "normal is not bound" << std::endl;
            return;
        }
        m_facetsAABB = facetsAABB;
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

        Eigen::Vector3d N = (p2_vec3 - p1_vec3).cross(p3_vec3 - p1_vec3);
        N.normalize();

//        if (m_facetsAABB != nullptr) {
//            auto mid = (p1_vec3 + p2_vec3 + p3_vec3) / 3.0;
//            auto f_id = m_facetsAABB->nearest_facet(mid);
//            auto f_n = GEO::vec3(m_face_normal[3 * f_id], m_face_normal[3 * f_id + 1], m_face_normal[3 * f_id + 2]);
//            if (dot(N, f_n) < 0) {
//                N = -N;
//            }
//        }

        Eigen::Matrix3d N_mat3 = N * N.transpose();
        N_mat3 = N_mat3 * 6;

        Eigen::MatrixXd M;
        M = Eigen::MatrixXd::Identity(m_dim, m_dim);
        M.block<3, 3>(0, 0) += N_mat3;


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
        Eigen::VectorXd dFdU1, dFdU2, dFdU3, dFdp1, dFdp2, dFdp3;
        dFdp1 = M * (E * dTdU1 + T * dEdU1);
        dFdp2 = M * (E * dTdU2 + T * dEdU2);
        dFdp3 = M * (E * dTdU3 + T * dEdU3);

        for (auto i = 0; i < m_dim; i++) {
            g_[m_dim * v + i] += -dFdp1[i] - dFdp2[i] - dFdp3[i];
        }

        double f = T * E;
        return f;
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

        if (A2 < 0.0) {
            for (auto i = 0; i < dim; i++) {
                dTdU1[i] = 0.0;
                dTdU2[i] = 0.0;
                dTdU3[i] = 0.0;
            }
            return 0.0;
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

    void LpCVTIS::vecmul(const double *p1, const double *p2, double *to) {
        for (auto i = 0; i < m_dim; i++) {
            to[i] = p1[i] * p2[i];
        }
    }

    void LpCVTIS::vecmul(const double *p1, const double *p2, const double *p3, double *to) {
        for (auto i = 0; i < m_dim; i++) {
            to[i] = p1[i] * p2[i] * p3[i];
        }
    }

    double LpCVTIS::vecbar(const double *p1) {
        double result = 0.0;
        for (auto i = 0; i < m_dim; i++) {
            result += p1[i];
        }
        return result;
    }

    void LpCVTIS::vecmadd(double s, const double *p1, const double *p2, const double *p3, double *to) {
        for (auto i = 0; i < m_dim; i++) {
            to[i] += s * p1[i] * p2[i] * p3[i];
        }
    }

    void LpCVTIS::vecmadd(double s, const GEO::vec3 &p1, double t, const GEO::vec3 &p2, GEO::vec3 &to) {
        to.x = s * p1.x + t * p2.x;
        to.y = s * p1.y + t * p2.y;
        to.z = s * p1.z + t * p2.z;
    }

    void LpCVTIS::matTvecmul(const GEO::mat3 &M, const GEO::vec3 &U, GEO::vec3 &V) {
        V.x = M(0, 0) * U.x + M(1, 0) * U.y + M(2, 0) * U.z;
        V.y = M(0, 1) * U.x + M(1, 1) * U.y + M(2, 1) * U.z;
        V.z = M(0, 2) * U.x + M(1, 2) * U.y + M(2, 2) * U.z;
    }
}