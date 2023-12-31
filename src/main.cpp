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

    std::string input_filename = "/Users/lihongbo/Desktop/code/LpCVT_recon/tmp/126660_sf_norm.obj";
    std::string hd_filename = "/Users/lihongbo/Desktop/code/LpCVT_recon/tmp/126660_emb.obj";
    std::string output_filename = "/Users/lihongbo/Desktop/code/LpCVT_recon/result/126660_tri.obj";
    const int dim = 8;

    std::shared_ptr<LpCVT::Mesh> mesh = std::make_shared<LpCVT::Mesh>();
    mesh->LoadMesh(input_filename);
    mesh->CalculateFaceNormal();

    GEO::Mesh M_3d;
    LpCVT::MeshAdaptor::Convert(*mesh, M_3d);

    GEO::Mesh M_hd;
    LpCVT::MeshAdaptor::HdMeshLoad(hd_filename, M_hd, dim);

    LpCVT::Remesher remesher;
    remesher.Init(M_hd, dim, LpCVT::Remesher::RemeshType::LPCVT);
    remesher.Remeshing(20000, 100);
    GEO::Mesh M_out;
    remesher.GetRDT(M_out);

    LpCVT::MeshAdaptor::SaveGEOMesh(output_filename, M_out);
}
