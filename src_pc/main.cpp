//
// Created by hongbo on 14/11/24.
//
#include <iostream>
#include <pcl/point_types.h>
#include <geogram/basic/command_line_args.h>
#include "Base/PointCloud.h"
#include "CVT_PC/Co3ne_NM.h"
#include "Base/MeshAdaptor.h"
int main(int argc, char **argv){
    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("remesh");
    GEO::CmdLine::import_arg_group("post");
    GEO::CmdLine::import_arg_group("poly");
    GEO::CmdLine::import_arg_group("co3ne");


    auto pc_path = "/home/hongbo/Desktop/code/tmp/LpCVT_recon/res/01_82-block.xyz";
    auto output_filepath = "/home/hongbo/Desktop/code/tmp/LpCVT_recon/tmp/block.obj";
    NMCVT_PC::PointCloud pc;
    pc.LoadFromXYZ(pc_path);

    auto points = pc.GetPoints();
    NMCVT_PC::Co3ne_NM co3ne_nm;
    GEO::Mesh M;
    co3ne_nm.Reconstruct(points,M);
    LpCVT::MeshAdaptor::SaveGEOMesh(output_filepath,M);
    return 0;
}