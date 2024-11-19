//
// Created by hongbo on 18/11/24.
//

#ifndef LPCVT_RECON_NNSEARCH_H
#define LPCVT_RECON_NNSEARCH_H
#include <vector>
#include <geogram/mesh/mesh.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/basic/geometry_nd.h>

namespace Base{
    class NNSearch{
    public:
        NNSearch(int dim);
        ~NNSearch()=default;
        void ConstructKDTree(GEO::Mesh& mesh);

        std::vector<double> GetVertex(unsigned int i);

        template<int DIM>
        void FindNearestPointOnSurface(std::vector<double> query, std::vector<double>& nearest) {
            GEO::index_t v = m_delaunay_nn->nearest_vertex(query.data());
            nearest =GetVertex(v);

            double d2= GEO::Numeric::max_float64();
            GEO::vecng<DIM,double> P(query.data());

            for(GEO::index_t i = 0; i < m_2stars[v].size(); i++) {
                GEO::index_t t = m_2stars[v][i];
                auto p1 = GEO::vecng<DIM,double>(GetVertex(m_triangles[3*t]).data());
                auto p2 = GEO::vecng<DIM,double>(GetVertex(m_triangles[3 * t + 1]).data());
                auto p3 = GEO::vecng<DIM,double>(GetVertex(m_triangles[3 * t + 2]).data());

                double l1, l2, l3;
                GEO::vecng<DIM,double> nearestP;

                double cur_d2 = GEO::Geom::point_triangle_squared_distance(
                        P, p1, p2, p3, nearestP, l1, l2, l3
                );
                if(cur_d2 < d2) {
                    d2 = cur_d2;
//                    const GEO::vec3& p1_R3{p1[0], p1[1], p1[2]};
//                    const GEO::vec3& p2_R3{p2[0], p2[1], p2[2]};
//                    const GEO::vec3& p3_R3{p3[0], p3[1], p3[2]};
//
//                    auto barycenter = l1 * p1_R3 + l2 * p2_R3 + l3 * p3_R3;
//                    for (int j = 0; j < 3; j++) {
//                        nearest[j] = barycenter[j];
//                    }
                    auto barycenter = l1 * p1 + l2 * p2 + l3 * p3;
                    for (int j = 0; j < DIM; j++) {
                        nearest[j] = barycenter[j];
                    }
                }
            }
        }

        void QueryPointSet(GEO::vector<double> query, GEO::vector<double>& nearest);


    private:
        int m_dim;
        GEO::Mesh m_mesh;
        GEO::Delaunay_var m_delaunay_nn;
        std::vector<double> m_vertices;
        std::vector<GEO::index_t> m_triangles;
        GEO::index_t m_nb_triangles;
        std::vector<std::vector<GEO::index_t>> m_1stars;
        std::vector<std::vector<GEO::index_t>> m_2stars;

    public:
        int GetDimension(){
            return m_dim;
        }
    };

}
#endif //LPCVT_RECON_NNSEARCH_H
