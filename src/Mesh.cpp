//
// Created by Hongbo on 12/30/23.
//
#include "Base/Mesh.h"
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/avg_edge_length.h>
namespace LpCVT {
    bool Mesh::LoadMesh(const std::string &filename) {
        igl::read_triangle_mesh(filename, m_v, m_f);
        return true;
    }

    bool Mesh::CalculateFaceNormal() {
        igl::per_face_normals(m_v, m_f, m_fn);
        return true;
    }

    void Mesh::PlotMesh() {
        igl::opengl::glfw::Viewer viewer;
        viewer.data().set_mesh(m_v, m_f);
        viewer.data().set_normals(m_fn);

        Eigen::MatrixXd centroids(m_f.rows(), 3);

        for (int i = 0; i < m_f.rows(); ++i)
        {
            // Get vertices of the current face
            Eigen::RowVector3d v1 = m_v.row(m_f(i, 0));
            Eigen::RowVector3d v2 = m_v.row(m_f(i, 1));
            Eigen::RowVector3d v3 = m_v.row(m_f(i, 2));

            // Calculate centroid as 1/3 of the sum of vertices
            Eigen::RowVector3d centroid = (v1 + v2 + v3) / 3.0;

            // Store centroid in the centroids matrix
            centroids.row(i) = centroid;
        }

        const double avg = igl::avg_edge_length(m_v, m_f);
        const Eigen::RowVector3d red(0.8, 0.2, 0.2);
        viewer.data().add_edges(centroids, centroids + m_fn * avg, red);

        viewer.launch();
    }
}