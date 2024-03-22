//
// Created by Hongbo on 12/30/23.
//
#include "Base/Mesh.h"
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/avg_edge_length.h>
#include <igl/principal_curvature.h>
#include <igl/copyleft/comiso/frame_field.h>
#include <igl/frame_to_cross_field.h>
#include <igl/copyleft/comiso/nrosy.h>
#include <igl/rotate_vectors.h>


namespace LpCVT {
    static Eigen::RowVector3d RED(0.8, 0.2, 0.2);
    static Eigen::RowVector3d BLUE(0.2, 0.2, 0.8);

    bool Mesh::LoadMesh(const std::string &filename) {
        igl::read_triangle_mesh(filename, m_v, m_f);
        return true;
    }

    bool Mesh::CalculateCurvature() {
        igl::per_vertex_normals(m_v, m_f, m_vn);
        igl::per_face_normals(m_v, m_f, m_fn);
        igl::principal_curvature(m_v, m_f, m_max_pd, m_min_pd, m_max_pv, m_min_pv);
        return true;
    }

    void Mesh::CalculateCrossField() {
        if (m_min_pd.size() == 0) {
            CalculateCurvature();
        }

        Eigen::VectorXi constrain_idx;
        constrain_idx.resize(1);
        constrain_idx(0) = 0;

//        constrain_idx.resize(m_f.rows());
//        for (auto i = 0; i < m_f.rows(); i++) {
//            constrain_idx(i) = i;
//        }

//        Eigen::MatrixXd f_min_pd, f_max_pd;
//        igl::barycenter(m_min_pd, m_f, f_min_pd);
//        igl::barycenter(m_max_pd, m_f, f_max_pd);


//        Eigen::MatrixXd frame_w, frame_v;
//        igl::copyleft::comiso::frame_field(m_v,
//                                           m_f,
//                                           constrain_idx,
//                                           f_min_pd,
//                                           f_max_pd,
//                                           frame_w,
//                                           frame_v);
//
//        igl::frame_to_cross_field(m_v, m_f,
//                                  frame_w, frame_v,
//                                  m_cross_field_w);

        Eigen::MatrixXd bc_x(constrain_idx.size(), 3);
//        for (unsigned i = 0; i < constrain_idx.size(); ++i)
//            bc_x.row(i) = f_min_pd.row(constrain_idx(i));
        bc_x.row(0) = Eigen::RowVector3d(1, 0, 0);

        Eigen::VectorXd S;
        igl::copyleft::comiso::nrosy(
                m_v,
                m_f,
                constrain_idx,
                bc_x,
                4,
                m_cross_field_w,
                S);

        Eigen::MatrixXd B1, B2, B3;
        igl::local_basis(m_v, m_f, B1, B2, B3);
        m_cross_field_v =
                igl::rotate_vectors(m_cross_field_w,
                                    Eigen::VectorXd::Constant(1, igl::PI / 2),
                                    B1, B2);
    }

    bool Mesh::NormalizeMesh() {
        Eigen::RowVector3d min_corner, max_corner;
        min_corner = m_v.colwise().minCoeff();
        max_corner = m_v.colwise().maxCoeff();
        auto diag_size = (max_corner - min_corner).norm();
        if (diag_size > 1.2 || diag_size < 0.8) {
            auto mid_point = (max_corner + min_corner) / 2;
            auto scale = 1.0 / diag_size;
            m_v.rowwise() -= mid_point;
            m_v *= scale;
        }
        return true;
    }

    void Mesh::PlotMesh() {
        igl::opengl::glfw::Viewer viewer;
        viewer.data().set_mesh(m_v, m_f);

        const double avg = igl::avg_edge_length(m_v, m_f);
        // Draw a red segment parallel to the maximal curvature direction
        viewer.data().add_edges(m_v, m_v + m_min_pd * avg, RED);

        // Draw a blue segment parallel to the minimal curvature direction
        viewer.data().add_edges(m_v, m_v + m_max_pd * avg, BLUE);

        // Hide wireframe
        viewer.data().show_lines = false;

        viewer.launch();
    }
}