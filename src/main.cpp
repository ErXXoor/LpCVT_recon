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
    int nb_pts = 13435;
    int dim = 8;
    double metric_weight = 6.0;
    int degree = 2;
    bool post_process = false;

    if (argc == 9) {
        input_filename = argv[1];
        output_filename = argv[2];
        dim = std::stoi(argv[3]);
        iteration = std::stoi(argv[4]);
        nb_pts = std::stoi(argv[5]);
        degree = std::stoi(argv[6]);
        metric_weight = std::stod(argv[7]);
        post_process = std::string(argv[8]) == "true";
    } else {
        std::cout << "Parameter length error" << std::endl;
        return 0;
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
    remesher.Init(M, dim, degree, metric_weight, LpCVT::Remesher::RemeshType::LPCVT_NORMAL);
    remesher.Remeshing(nb_pts, iteration);
    GEO::Mesh M_out;
    remesher.GetRDT(M_out, post_process);

    LpCVT::MeshAdaptor::SaveGEOMesh(output_filename, M_out);
}
