#include <stan/math.hpp> // Should be included first
#include <iostream>
#include <geogram/mesh/mesh.h>
#include "Base/MeshAdaptor.h"
#include "Base/Mesh.h"
#include "CVT/Remesher.h"
#include "CVT/LpCVTIS.h"
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/basic/smart_pointer.h>

int main(int argc, char **argv) {
    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("remesh");
    GEO::CmdLine::import_arg_group("post");
    GEO::CmdLine::import_arg_group("poly");

    std::string input_filename;
    std::string output_filename;
    std::string input_ori_filename;
    int iteration = 400;
    int nb_pts = 13435;
    int dim = 8;
    double metric_weight = 6.0;
    int degree = 2;
    bool post_process = false;

    if (argc == 10) {
        input_filename = argv[1];
        output_filename = argv[2];
        dim = std::stoi(argv[3]);
        iteration = std::stoi(argv[4]);
        nb_pts = std::stoi(argv[5]);
        degree = std::stoi(argv[6]);
        metric_weight = std::stod(argv[7]);
        post_process = std::string(argv[8]) == "true";
        input_ori_filename = argv[9];
    } else {
        std::cout << "Parameter length error" << std::endl;
        return 0;
    }


    GEO::Mesh M;
    std::shared_ptr<LpCVT::Mesh> mesh;
    std::shared_ptr<LpCVT::Mesh> mesh_ori;
    if (dim == 3) {
        mesh = std::make_shared<LpCVT::Mesh>();
        mesh->LoadMesh(input_filename);
        LpCVT::MeshAdaptor::Convert(*mesh, M);
    } else {
        LpCVT::MeshAdaptor::HdMeshLoad(input_filename, M, dim);
        mesh_ori = std::make_shared<LpCVT::Mesh>();
        mesh_ori->LoadMesh(input_ori_filename);
    }

    GEO::SmartPointer<LpCVT::LpCVTIS> is;
    auto remesh_type = LpCVT::Remesher::RemeshType::LPCVT_NORMAL;
    if (remesh_type == LpCVT::Remesher::RemeshType::LPCVT) {
        is = new LpCVT::LpCVTIS(M, false, dim, degree, metric_weight);
    } else if (remesh_type == LpCVT::Remesher::RemeshType::LPCVT_NORMAL) {
        is = new LpCVT::LpCVTIS(M, false, dim, degree, metric_weight);
        auto facet_aabb = std::make_shared<GEO::MeshFacetsAABB>();
        facet_aabb->initialize(M, false);


        if(dim>3){
            mesh_ori->CalculateCurvature();
            is->set_mesh(mesh_ori);
        }
        else{
            mesh->CalculateCurvature();
            is->set_mesh(mesh);
        }
        is->CalQuadMetric();
        is->set_facetsAABB(facet_aabb);
    }

    LpCVT::Remesher remesher;
    remesher.Init(M, is, dim);
    remesher.Remeshing(nb_pts, iteration);
    GEO::Mesh M_out;
//    remesher.GetRDT(M_out, post_process);
    remesher.GetRVD(M_out);

    LpCVT::MeshAdaptor::SaveGEOMesh(output_filename, M_out);

}
