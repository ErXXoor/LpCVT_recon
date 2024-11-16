//
// Created by hongbo on 16/11/24.
//
#include "CVT_PC/Co3neThread.h"
#include "CVT_PC/Co3neT.h"
namespace NMCVT_PC{
Co3neThread::Co3neThread(NMCVT_PC::Co3neT *master, GEO::index_t from, GEO::index_t to):
    master_(master), from_(from), to_(to)
{
//    mode_ = Co3neT::CO3NE_NONE;
}

void Co3neThread::run_reconstruct() {
    Co3NeRestrictedVoronoiDiagram& RVD = master_->RVD();
    vector<index_t> neigh(100);
    vector<double> sq_dist(100);
    Co3NeRestrictedVoronoiDiagram::Polygon P(100);
    Co3NeRestrictedVoronoiDiagram::Polygon Q(100);

    for(index_t i = from_; i < to_; i++) {
        RVD.get_RVC(i, P, Q, neigh, sq_dist);
        for(index_t v1 = 0; v1 < P.nb_vertices(); v1++) {
            index_t v2 = P.next_vertex(v1);
            signed_index_t j = P.vertex(v1).adjacent_seed();
            signed_index_t k = P.vertex(v2).adjacent_seed();
            if(
                j >= 0 && k >= 0 && j != k
            ) {
                triangles_.push_back(i);
                triangles_.push_back(index_t(j));
                triangles_.push_back(index_t(k));
            }
        }
    }
}

void Co3neThread::run_normals_and_reconstruct() {

    Attribute<double> normal;
//    if(CmdLine::get_arg_bool("co3ne:use_normals")) {
        Process::enter_critical_section();
        normal.bind_if_is_defined(
            master_->mesh().vertices.attributes(), "normal"
        );
        Process::leave_critical_section();
//    }

    std::ofstream RVD_file;
    bool debug_RVD = false;
    //if(
    //    CmdLine::get_arg_bool("dbg:co3neRVD")
    //) {
    //    if(CmdLine::get_arg_bool("sys:multithread")) {
    //        Logger::warn("Co3Ne")
    //            << "dbg:Co3NeRVD cannot work in multithread mode"
    //            << std::endl;
    //        Logger::warn("Co3Ne")
    //            << "use sys:multithread=false"
    //            << std::endl;
    //    } else {
    //        Logger::out("Co3Ne") << "Saving RVD in co3neRVD.obj"
    //                             << std::endl;
    //        RVD_file.open("co3neRVD.obj");
    //        debug_RVD=true;
    //    }
    //}
    index_t cur_v = 0;
    index_t tcount = 0;

    Co3NeRestrictedVoronoiDiagram& RVD = master_->RVD();
    index_t nb_neigh = RVD.nb_neighbors();
    vector<index_t> neigh(100);
    vector<double> sq_dist(100);
    Co3NeRestrictedVoronoiDiagram::Polygon P(100);
    Co3NeRestrictedVoronoiDiagram::Polygon Q(100);

    for(index_t i = from_; i < to_; i++) {

        vec3 N;
        if(normal.is_bound()) {
            RVD.get_neighbors(
                i, neigh, sq_dist, nb_neigh
            );
            N = vec3(normal[3*i], normal[3*i+1], normal[3*i+2]);
        } else {
            RVD.get_neighbors(
                i, neigh, sq_dist, nb_neigh
            );
            least_squares_normal_.begin();
            for(index_t jj = 0; jj < neigh.size(); jj++) {
                least_squares_normal_.add_point(RVD.point(neigh[jj]));
            }
            least_squares_normal_.end();
            N = least_squares_normal_.normal();
        }

        RVD.get_RVC(i, N, P, Q, neigh, sq_dist);
        if(debug_RVD) {
            for(index_t v = 0; v < P.nb_vertices(); ++v) {
                RVD_file << "v "
                         << P.vertex(v).point().x
                         << " "
                         << P.vertex(v).point().y
                         << " "
                         << P.vertex(v).point().z
                         << std::endl;
            }
            RVD_file << "f ";
            for(index_t v = 0; v < P.nb_vertices(); ++v) {
                ++cur_v;
                RVD_file << cur_v << " ";
            }
            RVD_file << std::endl;
            RVD_file << "#" << i << " ";
            for(index_t v1 = 0; v1 < P.nb_vertices(); ++v1) {
                RVD_file << P.vertex(v1).adjacent_seed() << " ";
            }
            RVD_file << std::endl;
        }
        for(index_t v1 = 0; v1 < P.nb_vertices(); v1++) {
            index_t v2 = P.next_vertex(v1);
            signed_index_t j = P.vertex(v1).adjacent_seed();
            signed_index_t k = P.vertex(v2).adjacent_seed();
            if(
                j >= 0 && k >= 0 && j != k
            ) {
                ++tcount;
                triangles_.push_back(i);
                triangles_.push_back(index_t(j));
                triangles_.push_back(index_t(k));
            }
        }
    }

    if(normal.is_bound()) {
        Process::enter_critical_section();
        normal.unbind();
        Process::leave_critical_section();
    }
}
}
