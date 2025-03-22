//
// Created by Hongbo on 12/30/23.
//
#include "CVT/Remesher.h"
#include <geogram/basic/logger.h>
#include <geogram/basic/progress.h>
#include "Base/NNSearch.h"

namespace LpCVT {
    void Remesher::Init(GEO::Mesh *M_in,
                        GEO::SmartPointer<LpCVTIS> is,
                        GEO::coord_index_t dim,
                        RemeshType type
    ) {
        GEO::Logger::div("CVT meshing");

        m_cvt = std::make_shared<LpCVTWrap>(M_in, dim);
        m_cvt->set_volumetric(false);
        is->set_Delaunay(m_cvt->delaunay());
        is->set_mesh_geo(M_in);
        //no cast func for GEO::IntegrationSimplex_var
        m_is = is;
        m_cvt->set_simplex_func(m_is);
        m_type = type;

    }

    void Remesher::Remeshing(unsigned int nb_pts,
                             unsigned int nb_iter,
                             std::vector<double> constrain_points) {
        GEO::Logger::div("Set metric weight");
        GEO::Logger::div("Generate random samples");

        m_cvt->compute_initial_sampling(nb_pts);

        if(constrain_points.size()>0){
            m_cvt->set_constrain_points(constrain_points);
        }


        GEO::Logger::div("Optimize sampling");
        GEO::ProgressTask progress("Optimize sampling", nb_iter);
        m_cvt->set_progress_logger(&progress);

        if (m_type == RemeshType::Lloyd_CVT) {
            m_cvt->Lloyd_iterations(nb_iter);
        } else {
            m_cvt->Newton_iterations(nb_iter);
        }

        m_cvt->set_progress_logger(nullptr);
    }

    void Remesher::GetRVD(GEO::Mesh &M_out) {
        GEO::Logger::div("Generate RVD");
        m_cvt->RVD()->set_exact_predicates(true);
        m_cvt->RVD()->compute_RVD(M_out, 3);
    }

    void Remesher::GetRDT(GEO::Mesh &M_out, bool post_process) {
        GEO::Logger::div("Generate RDT");

        if (post_process) {
            m_cvt->compute_surface(&M_out, true);
//            GEO::vector<double> points;
//            for (auto i = 0; i < M_out.vertices.nb(); i++) {
//                for (auto j = 0; j < M_out.vertices.dimension(); j++) {
//                    points.push_back(M_out.vertices.point_ptr(i)[j]);
//                }
//            }

//            Base::NNSearch search(3);
//            search.ConstructKDTree(*(m_cvt->mesh()));
//
//            GEO::vector<double> nearest;
//            search.QueryPointSet(points, nearest);
//
//            for(auto i=0;i<M_out.vertices.nb();i++){
//                for(auto j=0;j<M_out.vertices.dimension();j++){
//                    M_out.vertices.point_ptr(i)[j] = nearest[i*M_out.vertices.dimension()+j];
//                }
//            }

        } else {
            m_cvt->RVD()->set_exact_predicates(true);
            auto mode_i = GEO::RestrictedVoronoiDiagram::RDT_RVC_CENTROIDS |
                          GEO::RestrictedVoronoiDiagram::RDT_PREFER_SEEDS;
            auto mode = GEO::RestrictedVoronoiDiagram::RDTMode(mode_i);
            m_cvt->RVD()->compute_RDT(M_out, mode);
        }
    }

    void Remesher::GetHDRDT(GEO::Mesh &M_out) {
        GEO::vector<GEO::index_t> simplices;
        GEO::vector<double> embedding;
        m_cvt->RVD()->compute_RDT(simplices, embedding);

        M_out.vertices.set_double_precision();
        M_out.vertices.set_dimension(8);

        for(int i=0;i<embedding.size()-1;i+=8){
            std::vector<double> p(8);
            for(int j=0;j<8;j++){
                p[j] = embedding[i+j];
            }
            M_out.vertices.create_vertex(p.data());
        }

        for(auto i=0;i<simplices.size();i+=3) {
            M_out.facets.create_triangle(simplices[i], simplices[i + 1], simplices[i + 2]);
        }

        M_out.facets.connect();

    }

}