//
// Created by Hongbo on 12/30/23.
//
#include "Base/MeshAdaptor.h"
#include <geogram/mesh/mesh_io.h>

namespace LpCVT {
    void MeshAdaptor::Convert(const Mesh mesh_in, GEO::Mesh &M_out) {
        M_out.clear();
        M_out.vertices.set_double_precision();
        M_out.vertices.set_dimension(3);
        GEO::Attribute<double> face_normal;
        face_normal.create_vector_attribute(M_out.facets.attributes(), "normal", 3);

        for (auto i = 0; i < mesh_in.m_v.rows(); i++) {
            std::vector<double> p{mesh_in.m_v(i, 0), mesh_in.m_v(i, 1), mesh_in.m_v(i, 2)};
            M_out.vertices.create_vertex(p.data());
        }

        for (auto i = 0; i < mesh_in.m_f.rows(); i++) {
            M_out.facets.create_triangle(mesh_in.m_f(i, 0),
                                         mesh_in.m_f(i, 1),
                                         mesh_in.m_f(i, 2));
            if (mesh_in.m_fn.rows() > 0) {
                face_normal[3 * i] = mesh_in.m_fn(i, 0);
                face_normal[3 * i + 1] = mesh_in.m_fn(i, 1);
                face_normal[3 * i + 2] = mesh_in.m_fn(i, 2);
            }
        }
        M_out.facets.connect();
    }

    void MeshAdaptor::SaveGEOMesh(const std::string &filepath, const GEO::Mesh &M_out) {
        GEO::MeshIOFlags flags;
        flags.set_attribute(GEO::MESH_ALL_ATTRIBUTES);
        GEO::mesh_save(M_out, filepath, flags);
    }
}