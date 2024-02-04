#include <stan/math.hpp> // Should be included first
#include <iostream>
#include <geogram/mesh/mesh.h>
#include "Base/MeshAdaptor.h"
#include "Base/Mesh.h"
#include "CVT/Remesher.h"
#include <geogram/basic/command_line_args.h>


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
    int iteration = 400;
    int nb_pts = 1564;
    int dim = 8;
    bool post_process = false;

    if (argc == 7) {
        input_filename = argv[1];
        output_filename = argv[2];
        dim = std::stoi(argv[3]);
        iteration = std::stoi(argv[4]);
        nb_pts = std::stoi(argv[5]);
        post_process = std::string(argv[6]) == "true";
    } else {
        input_filename = "/media/hongbo/45ad552c-e83b-4f01-9864-7d87cfa1377e/hongbo/Thing10k_tetwild/hd_output_opt/237739/237739_emb.obj";
        output_filename = "/home/hongbo/Desktop/code/LpCVT_recon/tmp/237739_tri.obj";
//        std::cout<<"Parameter length error"<<std::endl;
//        return 0;
    }

    GEO::Mesh M;
    if (dim == 3) {
        std::shared_ptr<LpCVT::Mesh> mesh = std::make_shared<LpCVT::Mesh>();
        mesh->LoadMesh(input_filename);
        LpCVT::MeshAdaptor::Convert(*mesh, M);
    } else {
        LpCVT::MeshAdaptor::HdMeshLoad(input_filename, M, dim);
    }

    LpCVT::Remesher remesher;
    remesher.Init(M, dim, LpCVT::Remesher::RemeshType::LPCVT);
    remesher.Remeshing(nb_pts, iteration);
    GEO::Mesh M_out;
    remesher.GetRDT(M_out, post_process);

    LpCVT::MeshAdaptor::SaveGEOMesh(output_filename, M_out);
}
