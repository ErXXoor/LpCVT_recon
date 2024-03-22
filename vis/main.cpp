//
// Created by Hongbo on 2/17/24.
//

#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/frame_field_deformer.h>
#include <igl/frame_to_cross_field.h>
#include <igl/jet.h>
#include <igl/local_basis.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/rotate_vectors.h>
#include <igl/copyleft/comiso/nrosy.h>
#include <igl/copyleft/comiso/miq.h>
#include <igl/copyleft/comiso/frame_field.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/PI.h>
#include "Base/Mesh.h"

static Eigen::RowVector3d RED(0.8, 0.2, 0.2);
static Eigen::RowVector3d BLUE(0.2, 0.2, 0.8);

int main(int argc, char *argv[]) {
    using namespace Eigen;

    LpCVT::Mesh mesh;
    mesh.LoadMesh("/Users/lihongbo/Downloads/research/gt_dataset/hd_output_opt/85537/85537_tmp.obj");
    mesh.CalculateCurvature();
//    mesh.CalculateCrossField();
//
    auto V = mesh.GetVertices();
    auto F = mesh.GetFaces();
//
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
//
    const double avg = igl::avg_edge_length(V, F);
//
//    auto cross_w = mesh.GetCrossFieldW();
//    auto cross_v = mesh.GetCrossFieldV();
//
//    //Get center of each faces
//    MatrixXd C;
//    igl::barycenter(V, F, C);
//    MatrixXd f_pd_min, f_pd_max;
//    igl::barycenter(mesh.GetMinPD(), F, f_pd_min);
//    igl::barycenter(mesh.GetMaxPD(), F, f_pd_max);
//
//    Eigen::VectorXi constrain_idx;
//    constrain_idx.resize(V.rows());
//    for (auto i = 0; i < V.rows(); i++) {
//        constrain_idx(i) = i;
//    }
//
////    Eigen::MatrixXd frame_w, frame_v;
////    igl::copyleft::comiso::frame_field(V,
////                                       F,
////                                       constrain_idx,
////                                       f_pd_min,
////                                       f_pd_max,
////                                       frame_w,
////                                       frame_v);
//
//    Eigen::MatrixXd cross_w;
////    igl::frame_to_cross_field(V, F, frame_w, frame_v, cross_w);
//
//    MatrixXd bc_x(constrain_idx.size(), 3);
//    for (unsigned i = 0; i < constrain_idx.size(); ++i)
////        bc_x.row(i) = cross_w.row(constrain_idx(i));
//        bc_x.row(i) = f_pd_min.row(constrain_idx(i));
//
//    VectorXd S;
//    igl::copyleft::comiso::nrosy(
//            V,
//            F,
//            constrain_idx,
//            bc_x,
//            VectorXi(),
//            VectorXd(),
//            MatrixXd(),
//            4,
//            0,
//            cross_w,
//            S);
//
//    MatrixXd B1, B2, B3;
//    igl::local_basis(V, F, B1, B2, B3);
//    auto cross_v =
//            igl::rotate_vectors(cross_w, VectorXd::Constant(1, igl::PI / 2), B1, B2);

//    // Draw a red segment parallel to the maximal curvature direction
//    viewer.data().add_edges(C, C + cross_w * avg, RED);
//
//    // Draw a blue segment parallel to the minimal curvature direction
//    viewer.data().add_edges(C, C + cross_v * avg, BLUE);
//
//    // Hide wireframe
    viewer.data().show_lines = false;

    viewer.launch();

}