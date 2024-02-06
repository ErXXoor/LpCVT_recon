//
// Created by Hongbo on 12/30/23.
//
#include "CVT/Remesher.h"
#include <geogram/basic/logger.h>
#include "CVT/LpCVTIS.h"
#include <geogram/basic/progress.h>
#include <geogram/basic/smart_pointer.h>

namespace LpCVT {
    void Remesher::Init(const GEO::Mesh &M_in,
                        GEO::coord_index_t dim,
                        unsigned int degree,
                        double metric_weight,
                        RemeshType type) {
        GEO::Logger::div("CVT meshing");
        GEO::Mesh *mesh = new GEO::Mesh();
        mesh->copy(M_in);

        m_type = type;
        if (type == RemeshType::LPCVT_NORMAL) {
            m_facetsAABB = std::make_shared<GEO::MeshFacetsAABB>();
            m_facetsAABB->initialize(*mesh);
        }

        m_cvt = std::make_shared<LpCVTWrap>(mesh, dim);
        m_cvt->set_volumetric(false);

        if (type <= RemeshType::LPCVT) {
            m_is = new LpCVTIS(*mesh, false, dim, degree, metric_weight);
            m_cvt->set_simplex_func(m_is);
        } else if (type == RemeshType::LPCVT_NORMAL) {
            GEO::SmartPointer<LpCVTIS> is = new LpCVTIS(*mesh, false, dim, degree, metric_weight);
            is->set_facetsAABB(m_facetsAABB);
            m_is = is;
            m_cvt->set_simplex_func(m_is);
        }

    }

    void Remesher::Remeshing(unsigned int nb_pts,
                             unsigned int nb_iter) {
        GEO::Logger::div("Set metric weight");
        GEO::Logger::div("Generate random samples");

        m_cvt->compute_initial_sampling(nb_pts);

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
        } else {
            m_cvt->RVD()->set_exact_predicates(true);
            auto mode_i = GEO::RestrictedVoronoiDiagram::RDT_RVC_CENTROIDS |
                          GEO::RestrictedVoronoiDiagram::RDT_PREFER_SEEDS;
            auto mode = GEO::RestrictedVoronoiDiagram::RDTMode(mode_i);
            m_cvt->RVD()->compute_RDT(M_out, mode);
        }
    }

}