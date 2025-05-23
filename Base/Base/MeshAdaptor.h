//
// Created by Hongbo on 12/30/23.
//

#ifndef LPCVT_RECON_MESHADAPTOR_H
#define LPCVT_RECON_MESHADAPTOR_H

#include "Base/Mesh.h"
#include <geogram/mesh/mesh.h>

namespace LpCVT {
    class MeshAdaptor {
    public:
        static void Convert(const Mesh &mesh_in, GEO::Mesh &M_out);
        static void Convert6D(const Mesh &mesh_in, GEO::Mesh &M_out,float aniso_ratio=0.3);

        static void SaveGEOMesh(const std::string &filepath, const GEO::Mesh &M_out);

        static void HdMeshLoad(const std::string &filepath, GEO::Mesh &mesh_out, int dim);

        static void HdMeshSave(const std::string &filepath, const GEO::Mesh &mesh_out);

        static void AttachAttributeFacet(Eigen::MatrixXd attr,
                                         const std::string &attr_name,
                                         int dim,
                                         GEO::Mesh &M);

        static void SaveVoronoiID(const GEO::Mesh &M_out,const std::string &filepath);

    };
}
#endif //LPCVT_RECON_MESHADAPTOR_H
