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

    std::string input_filename = "";
    std::string output_filename = "";
    int iteration = 400;
    int nb_pts = 7000;
    int dim = 3;
    bool post_process = false;

    if (argc == 7) {
        input_filename = argv[1];
        output_filename = argv[2];
        dim = std::stoi(argv[3]);
        iteration = std::stoi(argv[4]);
        nb_pts = std::stoi(argv[5]);
        post_process = std::string(argv[6]) == "true";
    } else {
        input_filename = "/Users/lihongbo/Desktop/code/LpCVT_recon/tmp/octa_flower_input.obj";
        output_filename = "/Users/lihongbo/Desktop/code/LpCVT_recon/result/octa_flower_tri.obj";
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
