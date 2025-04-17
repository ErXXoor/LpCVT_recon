//
// Created by Hongbo on 12/30/23.
//
#include "Base/MeshAdaptor.h"
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/line_stream.h>
#include <Eigen/Dense>

namespace LpCVT {
    void MeshAdaptor::Convert(const Mesh &mesh_in, GEO::Mesh &M_out) {
        M_out.clear();
        M_out.vertices.set_double_precision();
        M_out.vertices.set_dimension(3);

        for (auto i = 0; i < mesh_in.m_v.rows(); i++) {
            std::vector<double> p{mesh_in.m_v(i, 0), mesh_in.m_v(i, 1), mesh_in.m_v(i, 2)};
            M_out.vertices.create_vertex(p.data());
        }

        for (auto i = 0; i < mesh_in.m_f.rows(); i++) {
            M_out.facets.create_triangle(mesh_in.m_f(i, 0),
                                         mesh_in.m_f(i, 1),
                                         mesh_in.m_f(i, 2));
        }

        M_out.facets.connect();
    }

    void MeshAdaptor::Convert6D(const Mesh &mesh_in, GEO::Mesh &M_out,float aniso_ratio) {
        M_out.clear();
        M_out.vertices.set_double_precision();
        M_out.vertices.set_dimension(6);

        for (auto i = 0; i < mesh_in.m_v.rows(); i++) {
            std::vector<double> p{mesh_in.m_v(i, 0), mesh_in.m_v(i, 1), mesh_in.m_v(i, 2),
                                  aniso_ratio*mesh_in.m_vn(i, 0), aniso_ratio*mesh_in.m_vn(i, 1), aniso_ratio*mesh_in.m_vn(i, 2)};
            M_out.vertices.create_vertex(p.data());
        }

        for (auto i = 0; i < mesh_in.m_f.rows(); i++) {
            M_out.facets.create_triangle(mesh_in.m_f(i, 0),
                                         mesh_in.m_f(i, 1),
                                         mesh_in.m_f(i, 2));
        }

        M_out.facets.connect();

    }



    void MeshAdaptor::SaveGEOMesh(const std::string &filepath, const GEO::Mesh &M_out) {
        GEO::MeshIOFlags flags;
        flags.set_attribute(GEO::MESH_ALL_ATTRIBUTES);
        flags.set_dimension(8);
        GEO::mesh_save(M_out, filepath, flags);
    }

    void MeshAdaptor::HdMeshLoad(const std::string &filepath, GEO::Mesh &mesh_out, int dim) {
        mesh_out.clear();
        mesh_out.vertices.set_double_precision();
        mesh_out.vertices.set_dimension(dim);
        GEO::LineInput in(filepath);
        if (!in.OK()) {
            std::cerr << "Error: cannot open file " << filepath << std::endl;
            return;
        }
        std::vector<double> P(dim);
        std::vector<GEO::index_t> facet_vertices;
        while (!in.eof() && in.get_line()) {
            in.get_fields();
            if (in.nb_fields() >= 1) {
                if (in.field_matches(0, "v")) {
                    for (GEO::coord_index_t c = 0; c < dim; c++) {
                        if (GEO::index_t(c + 1) < in.nb_fields()) {
                            P[c] = in.field_as_double(GEO::index_t(c + 1));
                        } else {
                            P[c] = 0.0;
                        }
                    }
                    GEO::index_t v = mesh_out.vertices.create_vertex();
                    double *p = mesh_out.vertices.point_ptr(v);
                    for (GEO::index_t c = 0; c < dim; ++c) {
                        p[c] = P[c];
                    }
                } else if (in.field_matches(0, "f")) {
                    if (in.nb_fields() < 3) {
                        GEO::Logger::err("I/O")
                                << "Line " << in.line_number()
                                << ": facet only has " << in.nb_fields()
                                << " corners (at least 3 required)"
                                << std::endl;
                        return;
                    }

                    facet_vertices.resize(0);
                    for (GEO::index_t i = 1; i < in.nb_fields(); i++) {
                        // In .obj files,
                        // negative vertex index means
                        // nb_vertices - vertex index
                        int s_vertex_index = in.field_as_int(i);
                        GEO::index_t vertex_index = 0;
                        if (s_vertex_index < 0) {
                            vertex_index = GEO::index_t(
                                    1 + int(mesh_out.vertices.nb()) + s_vertex_index
                            );
                        } else {
                            vertex_index = GEO::index_t(s_vertex_index);
                        }
                        facet_vertices.push_back(vertex_index - 1);
                    }

                    GEO::index_t f = mesh_out.facets.create_polygon(
                            facet_vertices.size()
                    );
                    for (GEO::index_t lv = 0; lv < facet_vertices.size(); ++lv) {
                        mesh_out.facets.set_vertex(f, lv, facet_vertices[lv]);
                    }
                }
            }
        }
        mesh_out.facets.connect();
    }

    void MeshAdaptor::HdMeshSave(const std::string &filepath, const GEO::Mesh &mesh_out) {
        auto dim = mesh_out.vertices.dimension();
        std::ofstream out(filepath);
        if (!out.is_open()) {
            std::cerr << "Error: cannot open file " << filepath << std::endl;
            return;
        }
        for (GEO::index_t v = 0; v < mesh_out.vertices.nb(); v++) {
            out << "v";
            for (GEO::index_t c = 0; c < dim; c++) {
                out << " " << mesh_out.vertices.point_ptr(v)[c];
            }
            out << std::endl;
        }

        for (GEO::index_t f = 0; f < mesh_out.facets.nb(); f++) {
            out << "f";
            for (GEO::index_t lv = 0; lv < mesh_out.facets.nb_vertices(f); lv++) {
                out << " " << mesh_out.facets.vertex(f, lv) + 1;
            }
            out << std::endl;
        }
        out.close();
    }

    void MeshAdaptor::AttachAttributeFacet(Eigen::MatrixXd attr,
                                           const std::string &attr_name,
                                           int dim,
                                           GEO::Mesh &M) {
        GEO::Attribute<double> facet_attr;
        facet_attr.create_vector_attribute(M.facets.attributes(), attr_name, dim);
        for (auto i = 0; i < M.facets.nb(); i++) {
            for (auto j = 0; j < dim; j++) {
                facet_attr[i * dim + j] = attr(i, j);
            }
        }
    }

    void MeshAdaptor::SaveVoronoiID(const GEO::Mesh &M_out,const std::string &filepath) {
        GEO::Logger::div("Get Voronoi ID");
        GEO::Attribute<GEO::index_t> facet_region;
        facet_region.bind_if_is_defined(M_out.facets.attributes(), "region");
        if(!facet_region.is_bound()){
            GEO::Logger::err("Vornoi_idx") << "Missing \"region\" facets attribute"
                                           << std::endl;
            return;
        }

        std::ofstream out(filepath);
        if (!out.is_open()) {
            std::cerr << "Error: cannot open file " << filepath << std::endl;
            return;
        }

        out<<"face_id,region"<<std::endl;
        for (auto i= 0; i < M_out.facets.nb(); i++) {
            out <<i<<","<< facet_region[i] << std::endl;
        }

    }
}