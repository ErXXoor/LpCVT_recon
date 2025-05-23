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
    std::string constrain_id;
    std::string vor_path;
    auto vor_region_path="";


    if (argc == 12) {
        input_filename = argv[1];
        output_filename = argv[2];
        dim = std::stoi(argv[3]);
        iteration = std::stoi(argv[4]);
        nb_pts = std::stoi(argv[5]);
        degree = std::stoi(argv[6]);
        metric_weight = std::stod(argv[7]);
        post_process = std::string(argv[8]) == "true";
        constrain_id = argv[9];
        vor_path = argv[10];
        vor_region_path = argv[11];
    } else {
        std::cout << "Parameter length error" << std::endl;
        return 1;
    }

    GEO::Mesh M;
    std::shared_ptr<LpCVT::Mesh> mesh;
    std::shared_ptr<LpCVT::Mesh> mesh_ori;
    if (dim == 3) {
        mesh = std::make_shared<LpCVT::Mesh>();
        mesh->LoadMesh(input_filename);
        LpCVT::MeshAdaptor::Convert(*mesh, M);

    }
    else if (dim==6){
        mesh = std::make_shared<LpCVT::Mesh>();
        mesh->LoadMesh(input_filename);
        mesh->CalculateNormal();
        LpCVT::MeshAdaptor::Convert6D(*mesh, M,0.1);
    }
    else {
        LpCVT::MeshAdaptor::HdMeshLoad(input_filename, M, dim);
    }

    std::istringstream stream(constrain_id);
    std::vector<int> constrain_ids;
    if (!constrain_id.empty()) {
        int id;
        while (stream >> id) {
            constrain_ids.push_back(id);
        }
    }

    std::vector<double> constrain_verts;
    for (int i = 0; i < constrain_ids.size(); i++) {
        for (int j = 0; j < dim; j++) {
            constrain_verts.push_back(M.vertices.point(constrain_ids[i])[j]);
        }
    }

    auto remesh_type = LpCVT::Remesher::RemeshType::L2CVT;
    GEO::SmartPointer<LpCVT::LpCVTIS> is;
    is = new LpCVT::LpCVTIS(M, false, dim, degree, metric_weight);

//    int ori_verts = M.vertices.nb();
//    nb_pts = ori_verts+3000;


    LpCVT::Remesher remesher;
    remesher.Init(&M, is, dim, remesh_type);

    remesher.Remeshing(nb_pts, iteration, constrain_verts);

    GEO::Mesh M_out;

//    remesher.GetHDRDT(M_out);
//    LpCVT::MeshAdaptor::HdMeshSave(output_filename, M_out);

    remesher.GetRDT(M_out, post_process);
    LpCVT::MeshAdaptor::SaveGEOMesh(output_filename, M_out);
//
//    GEO::Mesh M_vor;
//    if(!vor_path.empty()) {
//
//        remesher.GetRVD(M_vor);
//        LpCVT::MeshAdaptor::SaveGEOMesh(vor_path, M_vor);
//        LpCVT::MeshAdaptor::SaveVoronoiID(M_vor, vor_region_path);
//    }
}
