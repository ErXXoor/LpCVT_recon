//
// Created by Hongbo on 12/30/23.
//
#include "Base/Mesh.h"
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/avg_edge_length.h>
#include <igl/principal_curvature.h>
#include <igl/rotate_vectors.h>


namespace LpCVT {
    static Eigen::RowVector3d RED(0.8, 0.2, 0.2);
    static Eigen::RowVector3d BLUE(0.2, 0.2, 0.8);

    bool Mesh::LoadMesh(const std::string &filename) {
        igl::read_triangle_mesh(filename, m_v, m_f);
        return true;
    }

    bool Mesh::CalculateNormal() {
        igl::per_vertex_normals(m_v, m_f, m_vn);
        return true;
    }

    Eigen::MatrixXd Mesh::GetVertexNormals(){
        if(m_vn.rows() == 0)
            CalculateNormal();
        return m_vn;
    }

    bool Mesh::CalculateCurvature() {
        igl::per_vertex_normals(m_v, m_f, m_vn);
        igl::per_face_normals(m_v, m_f, m_fn);
        igl::principal_curvature(m_v, m_f, m_max_pd, m_min_pd, m_max_pv, m_min_pv);
        return true;
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