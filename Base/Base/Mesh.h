//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_MESH_H
#define LPCVT_RECON_MESH_H

#include <igl/read_triangle_mesh.h>

namespace LpCVT {
    class MeshAdaptor;

    class Mesh {
    public:
        Mesh() = default;

        ~Mesh() = default;

        bool LoadMesh(const std::string &filename);

        bool CalculateCurvature();

        bool NormalizeMesh();

        void PlotMesh();


        friend class MeshAdaptor;

    private:
        Eigen::MatrixXd m_v;
        Eigen::MatrixXi m_f;
        Eigen::MatrixXd m_vn;
        Eigen::MatrixXd m_fn;
        Eigen::MatrixXd m_min_pd;
        Eigen::MatrixXd m_max_pd;
        Eigen::VectorXd m_min_pv;
        Eigen::VectorXd m_max_pv;

        Eigen::MatrixXd m_cross_field_w;
        Eigen::MatrixXd m_cross_field_v;

    public:
        Eigen::MatrixXd GetVertices() const {
            return m_v;
        }

        Eigen::MatrixXi GetFaces() const {
            return m_f;
        }

        Eigen::MatrixXd GetMinPD() const {
            return m_min_pd;
        }

        Eigen::MatrixXd GetMaxPD() const {
            return m_max_pd;
        }

        Eigen::MatrixXd GetCrossFieldW() const {
            return m_cross_field_w;
        }

        Eigen::MatrixXd GetCrossFieldV() const {
            return m_cross_field_v;
        }

        Eigen::MatrixXd GetVertexNormals() const {
            return m_vn;
        }

        Eigen::MatrixXd GetFaceNormals() const {
            return m_fn;
        }

        std::vector<int> GetFace(int i) const {
            return {m_f(i, 0), m_f(i, 1), m_f(i, 2)};
        }
    };
}
#endif //LPCVT_RECON_MESH_H
