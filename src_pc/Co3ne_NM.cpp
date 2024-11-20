//
// Created by hongbo on 15/11/24.
//
#include "CVT_PC/Co3ne_NM.h"
#include "Base/MeshAdaptor.h"
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include "CVT_PC/Co3neT.h"
#include <geogram/points/kd_tree.h>


namespace NMCVT_PC{

    void Co3ne_NM::Reconstruct(Eigen::MatrixXd &points,unsigned int nb_pts,unsigned int nb_iter) {
        m_mesh = std::make_shared<Mesh>();
        LpCVT::MeshAdaptor::Convert(points, *m_mesh);
        double R = GEO::bbox_diagonal(*m_mesh);
        GEO::mesh_repair(*m_mesh, GEO::MESH_REPAIR_COLOCATE, 1e-6*R);

        auto radius = 0.01*R*5.0;
        m_cvt_core = std::make_shared<Co3neCVT>(m_mesh);
        m_is = new LpCVT::LpCVTIS(*m_mesh,false, 3,2,4);

        m_cvt_core->SurfaceRecon(radius);
//        m_cvt_core->InitRemeshing(m_is);
//        m_cvt_core->NMRemeshing(nb_pts, nb_iter);
    }

    void Co3ne_NM::Smooth(Mesh& M, int nb_neighbors, int nb_iterations) {
        Co3neT co3ne(M);
        try {
            ProgressTask progress("Smoothing", nb_iterations);
            for(int i = 0; i < nb_iterations; i++) {
                co3ne.smooth(nb_neighbors);
                if(i != nb_iterations - 1) {
                    co3ne.RVD().update();
                }
                progress.next();
            }
            co3ne.end_smooth();
        }
        catch(const TaskCanceled&) {
        }
    }

    void Co3ne_NM::GetRDT(GEO::Mesh &M_out, bool post_process) {
        m_cvt_core->GetRDT(M_out, post_process);
    }

    void Co3ne_NM::GetReconstructedMesh(LpCVT::Mesh &M_out) {
        LpCVT::MeshAdaptor::Convert(m_mesh, M_out);
    }

}