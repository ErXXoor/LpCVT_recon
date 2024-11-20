//
// Created by hongbo on 14/11/24.
//
#include <stan/math.hpp> // Should be included first
#include <iostream>
#include <pcl/point_types.h>
#include <geogram/basic/command_line_args.h>
#include "Base/PointCloud.h"
#include "CVT_PC/Co3ne_NM.h"
#include "Base/MeshAdaptor.h"
#include "Base/Mesh.h"
int main(int argc, char **argv){
    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("remesh");
    GEO::CmdLine::import_arg_group("post");
    GEO::CmdLine::import_arg_group("poly");
    GEO::CmdLine::import_arg_group("co3ne");


    auto pc_path = "/home/hongbo/Desktop/code/tmp/LpCVT_recon/tmp/pc.xyz";
    auto output_filepath = "/home/hongbo/Desktop/code/tmp/LpCVT_recon/tmp/guitar.obj";
    NMCVT_PC::PointCloud pc;
    pc.LoadFromXYZ(pc_path);

    auto points = pc.GetPoints();
    NMCVT_PC::Co3ne_NM co3ne_nm;
    co3ne_nm.Reconstruct(points,10000,400);
    GEO::Mesh M;
//    co3ne_nm.GetRDT(M);
    LpCVT::Mesh m_out;
    co3ne_nm.GetReconstructedMesh(m_out);
//    LpCVT::MeshAdaptor::SaveGEOMesh(output_filepath,m.get());
    LpCVT::MeshAdaptor::SaveBaseMesh(output_filepath,m_out);
    return 0;
}