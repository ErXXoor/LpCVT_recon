#include <iostream>
#include <geogram/mesh/mesh.h>
#include "Base/MeshAdaptor.h"
#include "Base/Mesh.h"
#include "CVT/Remesher.h"
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

int main() {
    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("remesh");
    GEO::CmdLine::import_arg_group("post");
    GEO::CmdLine::import_arg_group("poly");

    std::string input_filename = "/Users/lihongbo/Desktop/code/LpCVT_recon/tmp/datatech_1_input.obj";

    std::shared_ptr<LpCVT::Mesh> mesh = std::make_shared<LpCVT::Mesh>();
    mesh->LoadMesh(input_filename);
    mesh->CalculateFaceNormal();

    GEO::Mesh M;
    LpCVT::MeshAdaptor::Convert(*mesh, M);
    LpCVT::Remesher remesher;
    remesher.Init(M, 3, LpCVT::Remesher::RemeshType::LPCVT);
    remesher.Remeshing(7051, 100);
    GEO::Mesh M_out;
    remesher.GetRDT(M_out);

    std::string output_filename = "/Users/lihongbo/Desktop/code/LpCVT_recon/result/test.obj";
    LpCVT::MeshAdaptor::SaveGEOMesh(output_filename, M_out);
}
