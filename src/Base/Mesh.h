//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_MESH_H
#define LPCVT_RECON_MESH_H
#include <igl/read_triangle_mesh.h>

namespace LpCVT{
class Mesh {
public:
    Mesh() = default;

    ~Mesh() = default;

    bool LoadMesh(const std::string &filename);

    bool CalculateFaceNormal();

    void PlotMesh();

private:
    Eigen::MatrixXd m_v;
    Eigen::MatrixXi m_f;
    Eigen::MatrixXd m_fn;

};
}
#endif //LPCVT_RECON_MESH_H
