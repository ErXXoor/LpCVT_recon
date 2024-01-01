//
// Created by Hongbo on 12/30/23.
//
#include "CVT/LpCVTIS.h"
#include <geogram/voronoi/generic_RVD.h>
#include <geogram/basic/smart_pointer.h>
#include <geogram/basic/geometry_nd.h>

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
        auto p1_vec8 = GEO::vecng<8, double>(p1.point());
        auto p2_vec8 = GEO::vecng<8, double>(p2.point());
        auto p3_vec8 = GEO::vecng<8, double>(p3.point());

        auto p1_vec3 = GEO::vec3(p1.point()[0], p1.point()[1], p1.point()[2]);
        auto p2_vec3 = GEO::vec3(p2.point()[0], p2.point()[1], p2.point()[2]);
        auto p3_vec3 = GEO::vec3(p3.point()[0], p3.point()[1], p3.point()[2]);

        GEO::vec3 N = normalize(cross(p2_vec3 - p1_vec3, p3_vec3 - p1_vec3));
//        if (m_facetsAABB != nullptr) {
//            auto mid = (p1_vec3 + p2_vec3 + p3_vec3) / 3.0;
//            auto f_id = m_facetsAABB->nearest_facet(mid);
//            auto f_n = GEO::vec3(m_face_normal[3 * f_id], m_face_normal[3 * f_id + 1], m_face_normal[3 * f_id + 2]);
//            if (dot(N, f_n) < 0) {
//                N = -N;
//            }
//        }

        GEO::mat3 N_mat3;
        N_mat3(0, 0) = N.x * N.x;
        N_mat3(0, 1) = N.x * N.y;
        N_mat3(0, 2) = N.x * N.z;
        N_mat3(1, 0) = N.y * N.x;
        N_mat3(1, 1) = N.y * N.y;
        N_mat3(1, 2) = N.y * N.z;
        N_mat3(2, 0) = N.z * N.x;
        N_mat3(2, 1) = N.z * N.y;
        N_mat3(2, 2) = N.z * N.z;
        N_mat3 = N_mat3 * 4;

        GEO::Matrix<8, double> M;
        M.load_identity();
        for (auto i = 0; i < 3; i++) {
            for (auto j = 0; j < 3; j++) {
                M(i, j) += N_mat3(i, j);
            }
        }

        auto p0 = vertex_ptr(v);
        std::vector<std::vector<GEO::vecng<8, double>>> U_pow;
        std::vector<double> init_8d(8, 1.0);
        U_pow.resize(3);
        U_pow[0].resize(m_degree + 1);
        U_pow[1].resize(m_degree + 1);
        U_pow[2].resize(m_degree + 1);
        U_pow[0][0] = GEO::vecng<8, double>(init_8d.data());
        U_pow[1][0] = GEO::vecng<8, double>(init_8d.data());
        U_pow[2][0] = GEO::vecng<8, double>(init_8d.data());

        std::vector<double> u0;
        std::vector<double> u1;
        std::vector<double> u2;
        for (GEO::index_t c = 0; c < m_dim; c++) {
            u0.push_back(p1[c] - p0[c]);
            u1.push_back(p2[c] - p0[c]);
            u2.push_back(p3[c] - p0[c]);
        }

        mult(M, u0.data(), U_pow[0][1].data());
        mult(M, u1.data(), U_pow[1][1].data());
        mult(M, u2.data(), U_pow[2][1].data());

        for (unsigned int i = 2; i <= m_degree; i++) {
            vecmul(U_pow[0][1].data(), U_pow[0][i - 1].data(), U_pow[0][i].data());
            vecmul(U_pow[1][1].data(), U_pow[1][i - 1].data(), U_pow[1][i].data());
            vecmul(U_pow[2][1].data(), U_pow[2][i - 1].data(), U_pow[2][i].data());
        }

        // Computation of function value.
        double E = 0.0;
        for (unsigned int i = 0; i < nb_coeffs; i++) {
            GEO::vecng<8, double> W;
            unsigned int alpha = E_pow[i][0];
            unsigned int beta = E_pow[i][1];
            unsigned int gamma = E_pow[i][2];
            vecmul(U_pow[0][alpha].data(), U_pow[1][beta].data(), U_pow[2][gamma].data(), W.data());
            E += vecbar(W.data());
        }

        // Computation of gradient
        std::vector<double> dE_init(8, 0.0);
        GEO::vecng<8, double> dEdU1(dE_init.data()), dEdU2(dE_init.data()), dEdU3(dE_init.data());
        for (unsigned int i = 0; i < nb_dcoeffs; i++) {
            {
                unsigned int alpha = dE_pow[0][i][0];
                unsigned int beta = dE_pow[0][i][1];
                unsigned int gamma = dE_pow[0][i][2];
                vecmadd(alpha, U_pow[0][alpha - 1].data(), U_pow[1][beta].data(), U_pow[2][gamma].data(), dEdU1.data());
            }
            {
                unsigned int alpha = dE_pow[1][i][0];
                unsigned int beta = dE_pow[1][i][1];
                unsigned int gamma = dE_pow[1][i][2];
                vecmadd(beta, U_pow[0][alpha].data(), U_pow[1][beta - 1].data(), U_pow[2][gamma].data(), dEdU2.data());
            }
            {
                unsigned int alpha = dE_pow[2][i][0];
                unsigned int beta = dE_pow[2][i][1];
                unsigned int gamma = dE_pow[2][i][2];
                vecmadd(gamma, U_pow[0][alpha].data(), U_pow[1][beta].data(), U_pow[2][gamma - 1].data(), dEdU3.data());
            }
        }

        // Compute tet measure and its
        // derivatives relative to U1, U2 and U3.
        GEO::vec3 dTdU1, dTdU2, dTdU3;
//        double T = grad_tri(
//                U_pow[0][1], U_pow[1][1], U_pow[2][1],
//                dTdU1, dTdU2, dTdU3
//        );
        double T = GEO::Geom::triangle_area(U_pow[0][1].data(),
                                            U_pow[1][1].data(),
                                            U_pow[2][1].data(),
                                            m_dim);

        // Assemble dF = E.d|T| + |T|.dE
        GEO::vec3 dFdU1, dFdU2, dFdU3, dFdp1, dFdp2, dFdp3;
        vecmadd(E, dTdU1, T, dEdU1, dFdU1);
        matTvecmul(M, dFdU1, dFdp1);
        vecmadd(E, dTdU2, T, dEdU2, dFdU2);
        matTvecmul(M, dFdU2, dFdp2);
        vecmadd(E, dTdU3, T, dEdU3, dFdU3);
        matTvecmul(M, dFdU3, dFdp3);

        g_[3 * v + 0] += -dFdp1.x - dFdp2.x - dFdp3.x;
        g_[3 * v + 1] += -dFdp1.y - dFdp2.y - dFdp3.y;
        g_[3 * v + 2] += -dFdp1.z - dFdp2.z - dFdp3.z;

        double f = T * E;
//        for(index_t c = 0; c < m_dim; c++) {
//            double Gc = (1.0 / 3.0) * (p1[c] + p2[c] + p3[c]);
//            g_[m_dim * v + c] += (2.0 * t_area) * (p0[c] - Gc);
//        }
        return f;
    }

    double
    LpCVTIS::grad_tri(const GEO::vec3 &U1, const GEO::vec3 &U2, const GEO::vec3 &U3, GEO::vec3 &dTdU1, GEO::vec3 &dTdU2,
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