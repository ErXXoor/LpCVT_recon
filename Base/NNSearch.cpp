//
// Created by hongbo on 18/11/24.
//
#include "Base/NNSearch.h"
#include <geogram/basic/numeric.h>
namespace Base{
    void sort_unique(std::vector<GEO::index_t> &v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    NNSearch::NNSearch(int dim):m_dim(dim) {}

    void NNSearch::ConstructKDTree(GEO::Mesh& mesh) {
        // Step 1: get triangles
        for(GEO::index_t f = 0; f < mesh.facets.nb(); f++) {
            GEO::index_t i = mesh.facets.corners_begin(f);
            for(GEO::index_t j = i + 1;
                    j + 1 < mesh.facets.corners_end(f); j++
                    ) {
                m_triangles.push_back(mesh.facet_corners.vertex(i));
                m_triangles.push_back(mesh.facet_corners.vertex(j));
                m_triangles.push_back(mesh.facet_corners.vertex(j + 1));
            }
        }
        m_nb_triangles = GEO::index_t(m_triangles.size() / 3);

        // Step 2: get vertices stars
        //    Step 2.1: get one-ring neighborhood
        m_1stars.resize(mesh.vertices.nb());
        for(GEO::index_t t = 0; t < m_nb_triangles; t++) {
            m_1stars[m_triangles[3 * t]].push_back(t);
            m_1stars[m_triangles[3 * t + 1]].push_back(t);
            m_1stars[m_triangles[3 * t + 2]].push_back(t);
        }

        //   Step 2.2: get two-ring neighborhood
        m_2stars.resize(mesh.vertices.nb());
        for(GEO::index_t i = 0; i < m_1stars.size(); i++) {
            std::vector<GEO::index_t> Ni;
            for(GEO::index_t j = 0; j < m_1stars[i].size(); j++) {
                GEO::index_t t = m_1stars[i][j];
                for(GEO::index_t iv = 0; iv < 3; iv++) {
                    GEO::index_t v = m_triangles[3 * t + iv];
                    if(v != i) {
                        Ni.push_back(v);
                    }
                }
            }
            sort_unique(Ni);
            for(GEO::index_t j = 0; j < Ni.size(); j++) {
                GEO::index_t k = Ni[j];
                m_2stars[i].insert(
                        m_2stars[i].end(), m_1stars[k].begin(), m_1stars[k].end()
                );
            }
            sort_unique(m_2stars[i]);
        }

        // Step 3: create search structure
        m_delaunay_nn = GEO::Delaunay::create(m_dim, "NN");
        GEO::index_t nb_vertices = mesh.vertices.nb();

        m_vertices.resize(nb_vertices * m_dim);
        for(GEO::index_t i = 0; i < nb_vertices; i++) {
            for(GEO::index_t coord = 0; coord < m_dim; coord++) {
                m_vertices[i * m_dim + coord] =
                        mesh.vertices.point_ptr(i)[coord];
            }
        }
        m_delaunay_nn->set_vertices(nb_vertices, m_vertices.data());
    }

    std::vector<double> NNSearch::GetVertex(unsigned int i) {
        std::vector<double> vertex(m_dim);
        for(int j=0;j<m_dim;j++){
            vertex[j]=m_vertices[i*m_dim+j];
        }
        return vertex;
    }

    void NNSearch::QueryPointSet(GEO::vector<double> query,
                                 GEO::vector<double> &nearest) {
        int point_num = query.size() / m_dim;
        nearest.resize(point_num * m_dim);
        for(auto i=0;i<point_num;i++){
            std::vector<double> query_ptr(m_dim);
            std::vector<double> nearest_ptr(m_dim);
            for(auto j=0;j<m_dim;j++){
                query_ptr[j]=query[i*m_dim+j];
            }

            if(m_dim==8){
                FindNearestPointOnSurface<8>(query_ptr, nearest_ptr);
            }
            else if(m_dim==3){
                FindNearestPointOnSurface<3>(query_ptr, nearest_ptr);
            }
            else{
                std::cerr<<"The dimension of the point set is not supported!"<<std::endl;
            }
            for(auto j=0;j<m_dim;j++){
                nearest[i*m_dim+j]=nearest_ptr[j];
            }
        }
    }



}