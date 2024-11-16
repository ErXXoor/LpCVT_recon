//
// Created by hongbo on 16/11/24.
//

#ifndef LPCVT_RECON_CO3NET_H
#define LPCVT_RECON_CO3NET_H
#include <geogram/mesh/mesh.h>
#include <queue>
#include "Co3neThread.h"
#include "Co3neRVD.h"
#include "Co3neManifoldExtraction.h"
#include <geogram/mesh/mesh_repair.h>
#include "../CVT_UTILS/Utils.h"
#include <geogram/basic/progress.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_topology.h>

namespace NMCVT_PC{
using namespace GEO;
class Co3neT {
  public:
    enum Co3NeMode {
        CO3NE_NONE,    /**< uninitialized */
        CO3NE_NORMALS, /**< estimate normals in pointset */
        CO3NE_SMOOTH,  /**< smooth the pointset */
        CO3NE_RECONSTRUCT, /**< reconstruct the triangles */
        CO3NE_NORMALS_AND_RECONSTRUCT
        /**< combined normal estimation and reconstruction */
    };
    /**
         * \brief Constructs a new Co3Ne.
         * \param[in] M the pointset
     */
    Co3neT(Mesh& M) :
                     mesh_(M) {
        // TODO: interlace threads (more cache friendly)
        RVD_.init(mesh_);
        index_t nb = Process::maximum_concurrent_threads();
        thread_.clear();
        index_t batch_size = RVD_.nb_points() / nb;
        index_t cur = 0;
        index_t remaining = RVD_.nb_points();
        for(index_t i = 0; i < nb; i++) {
            index_t this_batch_size = batch_size;
            if(i == nb - 1) {
                this_batch_size = remaining;
            }
            thread_.push_back(
                new Co3neThread(
                    this, cur, cur + this_batch_size
                    )
            );
            cur += this_batch_size;
            remaining -= this_batch_size;
        }
        geo_assert(remaining == 0);

        // TODO: pass it as an argument and let Vorpaline's main.cpp
        // communicate with CmdLine.
        double alpha = 60.0;
        alpha = alpha * M_PI / 180.0;
        set_max_angle(alpha);
    }

    /**
	 * \brief Runs the threads.
     */
    void run_threads() {
        Process::run_threads(thread_);
    }

    /**
         * \brief Estimates the normals of the point set.
         * \details They are stored in the "normal" vertex attribute.
         * \param[in] nb_neighbors number of neighbors to be
         *  used for normal estimation
     */
    void compute_normals(index_t nb_neighbors) {
        Attribute<double> normals;
        normals.bind_if_is_defined(mesh_.vertices.attributes(), "normal");
        if(!normals.is_bound()) {
            normals.create_vector_attribute(
                mesh_.vertices.attributes(), "normal", 3
            );
        }
        RVD_.init(mesh_);
        RVD_.set_nb_neighbors(nb_neighbors);
//        for(index_t t = 0; t < thread_.size(); t++) {
//            thread_[t]->set_mode(CO3NE_NORMALS);
//        }
        run_threads();
    }

    static inline double cos_angle(
        Attribute<double>& normal, index_t v1, index_t v2
    ) {
        vec3 V1(normal[3*v1],normal[3*v1+1],normal[3*v1+2]);
        vec3 V2(normal[3*v2],normal[3*v2+1],normal[3*v2+2]);
        return Geom::cos_angle(V1,V2);
    }

    static inline void flip(Attribute<double>& normal, index_t v) {
        normal[3*v]   = -normal[3*v];
        normal[3*v+1] = -normal[3*v+1];
        normal[3*v+2] = -normal[3*v+2];
    }

    /**
	 * \brief Tentatively enforces a coherent orientation of normals
	 *  using a breadth-first traveral of the K-nearest-neighbor graph.
	 * \retval true if normals where computed
	 * \retval false otherwise (when the user pushes the cancel button).
     */
    bool reorient_normals() {
        Attribute<double> normal;
        normal.bind_if_is_defined(mesh_.vertices.attributes(), "normal");
        geo_assert(normal.is_bound());

        //  To resist noisy inputs, propagation is prioritized to the points
        // that have smallest normal deviations.

        std::priority_queue<OrientNormal> S;
        vector<index_t> neighbors(RVD_.nb_neighbors());
        vector<double> dist(RVD_.nb_neighbors());

        index_t nb=0;
        ProgressTask progress("Reorient");

        try {
            std::vector<bool> visited(mesh_.vertices.nb(), false);
            for(index_t v=0; v<mesh_.vertices.nb(); ++v) {
                if(!visited[v]) {
                    S.push(OrientNormal(v,0.0));
                    visited[v] = true;
                    while(!S.empty()) {
                        OrientNormal top = S.top();
                        ++nb;
                        progress.progress(nb*100/mesh_.vertices.nb());
                        S.pop();
                        if(top.dot < 0.0) {
                            flip(normal,top.v);
                        }
                        RVD_.get_neighbors(
                            top.v,
                            neighbors.data(),dist.data(),RVD_.nb_neighbors()
                        );
                        for(index_t i=0; i<RVD_.nb_neighbors(); ++i) {
                            index_t neigh = neighbors[i];
                            if(!visited[neigh]) {
                                visited[neigh] = true;
                                double dot =
                                    cos_angle(normal, top.v, neigh);
                                S.push(OrientNormal(neigh,dot));
                            }
                        }
                    }
                }
            }
        } catch(const TaskCanceled&) {
            return false;
        }
        return true;
    }

    /**
         * \brief Smoothes a point set by projection
         * onto the nearest neighbors best
         * approximating planes.
         * \param[in] nb_neighbors number of neighbors to be
         *  used for best approximating plane estimation
     */
    void smooth(index_t nb_neighbors) {
        new_vertices_.resize(mesh_.vertices.nb() * 3);
        RVD_.set_nb_neighbors(nb_neighbors);
//        for(index_t t = 0; t < thread_.size(); t++) {
//            thread_[t]->set_mode(CO3NE_SMOOTH);
//        }
        run_threads();
        /*
          // TODO: once 'steal-arg' mode works for vertices,
          // we can use this one.
        if(RVD_.p_stride_ == 3) {
            MeshMutator::vertices(mesh_).swap(new_vertices_);
        } else */ {
            index_t idx = 0;
            for(index_t i = 0; i < mesh_.vertices.nb(); i++) {
                double* p = mesh_.vertices.point_ptr(i);
                for(coord_index_t c = 0; c < 3; c++) {
                    p[c] = new_vertices_[idx];
                    idx++;
                }
            }
        }
    }

    /**
         * \brief This function needs to be called after the
         * last iteration of smoothing.
         * \details Deallocates the temporary
         *  variables used for smoothing.
     */
    void end_smooth() {
        new_vertices_.clear();
    }

    /**
         * \brief Reconstructs a mesh from a point set.
         * \details If the mesh has a "normal" vertex attribute,
         *  then the existing normals are used, else normals are estimated.
         * \param[in] r maximum distance used to determine
         *  points adjacencies.
     */
    void reconstruct(double r) {
        bool has_normals = false;
        {
            Attribute<double> normal;
            normal.bind_if_is_defined(mesh_.vertices.attributes(),"normal");
            has_normals = (
                normal.is_bound() && normal.dimension() == 3
            );
        }

        ProgressTask progress("reconstruct",100);

        if(has_normals) {
            Stopwatch W("Co3Ne recons");
            RVD_.set_circles_radius(r);
            for(index_t t = 0; t < thread_.size(); t++) {
//                thread_[t]->set_mode(CO3NE_RECONSTRUCT);
                thread_[t]->triangles().clear();
            }
            progress.progress(1);
            run_threads();
            progress.progress(50);
        } else {
            Stopwatch W("Co3Ne recons");
            Logger::out("Co3Ne")
                << "using combined \'normals and reconstruct\'"
                << std::endl;
            RVD_.set_nb_neighbors(
                30
            );
            RVD_.set_circles_radius(r);
            for(index_t t = 0; t < thread_.size(); t++) {
//                thread_[t]->set_mode(CO3NE_NORMALS_AND_RECONSTRUCT);
                thread_[t]->triangles().clear();
            }
            progress.progress(1);
            run_threads();
            progress.progress(50);
        }

        {
            Stopwatch W("Co3Ne manif.");
            RVD_.clear();  // reclaim memory used by ANN

            index_t nb_triangles = 0;
            for(index_t t = 0; t < thread_.size(); t++) {
                nb_triangles += thread_[t]->nb_triangles();
            }

            Logger::out("Co3Ne") << "Raw triangles: "
                                 << nb_triangles
                                 << std::endl;

            vector<index_t> raw_triangles;
            raw_triangles.reserve(nb_triangles * 3);
            for(index_t th = 0; th < thread_.size(); th++) {
                vector<index_t>& triangles = thread_[th]->triangles();
                raw_triangles.insert(
                    raw_triangles.end(),
                    triangles.begin(), triangles.end()
                );
                thread_[th]->triangles().clear();
            }

            vector<index_t> good_triangles;
            vector<index_t> not_so_good_triangles;
            co3ne_split_triangles_list(
                raw_triangles, good_triangles, not_so_good_triangles
            );


//            if(CmdLine::get_arg_bool("dbg:co3ne")) {
//                Logger::out("Co3Ne") << ">> co3ne_T3.geogram"
//                                     << std::endl;
//                Mesh M;
//                M.vertices.assign_points(
//                    mesh_.vertices.point_ptr(0),
//                    mesh_.vertices.dimension(),
//                    mesh_.vertices.nb()
//                );
//                M.facets.assign_triangle_mesh(good_triangles, false);
//                M.vertices.set_dimension(3);
//                mesh_save(M, "co3ne_T3.geogram");
//            }
//
//            if(CmdLine::get_arg_bool("dbg:co3ne")) {
//                Logger::out("Co3Ne") << ">> co3ne_T12.geogram"
//                                     << std::endl;
//                Mesh M;
//                M.vertices.assign_points(
//                    mesh_.vertices.point_ptr(0),
//                    mesh_.vertices.dimension(),
//                    mesh_.vertices.nb()
//                );
//                M.facets.assign_triangle_mesh(not_so_good_triangles, false);
//                M.vertices.set_dimension(3);
//                mesh_save(M, "co3ne_T12.geogram");
//            }

            progress.progress(53);

            Co3NeManifoldExtraction manifold_extraction(
                mesh_, good_triangles
            );

            progress.progress(55);

//            if(CmdLine::get_arg_bool("co3ne:T12")) {
                manifold_extraction.add_triangles(not_so_good_triangles);
//            }

            progress.progress(57);

            mesh_reorient(mesh_);

            progress.progress(60);

//            if(CmdLine::get_arg_bool("dbg:co3ne")) {
//                Logger::out("Co3Ne") << ">> co3ne_manif.geogram"
//                                     << std::endl;
//                mesh_save(mesh_, "co3ne_manif.geogram");
//            }
        }

//        if(CmdLine::get_arg_bool("co3ne:repair")) {
//            Stopwatch W("Co3Ne post.");
//            mesh_repair(mesh_,
//                        MeshRepairMode(
//                            MESH_REPAIR_DEFAULT | MESH_REPAIR_RECONSTRUCT
//                            )
//            );
//            if(CmdLine::get_arg_bool("dbg:co3ne")) {
//                Logger::out("Co3Ne") << ">> co3ne_post.geogram"
//                                     << std::endl;
//                mesh_save(mesh_, "co3ne_post.geogram");
//            }
//        }

        progress.progress(100);

        Logger::out("Topology")
            << "nb components=" << mesh_nb_connected_components(mesh_)
            << " nb borders=" <<  mesh_nb_borders(mesh_)
            << std::endl;

    }


//    void Co3Ne_smooth_and_reconstruct(
//        Mesh& M, index_t nb_neighbors, index_t nb_iterations, double radius
//    ) {
//        Stopwatch W("Co3Ne total");
//
////        if(CmdLine::get_arg_bool("co3ne:use_normals")) {
//            Attribute<double> normal;
//            normal.bind_if_is_defined(M.vertices.attributes(), "normal");
//            if(normal.is_bound() && normal.dimension() == 3) {
//                Logger::out("Co3Ne") << "Using existing normal attribute"
//                                     << std::endl;
//            } else {
//                Logger::out("Co3Ne") << "No \'normal\' vertex attribute found"
//                                     << std::endl;
//                Logger::out("Co3Ne") << "(estimating normals)"
//                                     << std::endl;
//            }
////        }
//
//
//        Co3neT co3ne(M);
//        if(nb_iterations != 0) {
//            try {
//                co3ne.RVD().set_exact(false);
//                ProgressTask progress("Co3Ne smooth", nb_iterations);
//                for(index_t i = 0; i < nb_iterations; i++) {
//                    co3ne.smooth(nb_neighbors);
//                    co3ne.RVD().update();
//                    progress.next();
//                }
//                co3ne.end_smooth();
//            }
//            catch(const TaskCanceled&) {
//                // TODO_CANCEL
//            }
//        }
//        Logger::out("Co3Ne") << "Reconstruction..." << std::endl;
//        co3ne.RVD().set_exact(true);
//        co3ne.reconstruct(radius);
//    }


    /**
         * \brief Gets the Co3NeRestrictedVoronoiDiagram associated
         *  with this Co3Ne.
         * \return a reference to the Co3NeRestrictedVoronoiDiagram
     */
    Co3NeRestrictedVoronoiDiagram& RVD() {
        return RVD_;
    }

    /**
         * \brief Sets a point
         * \param[in] i the index of the point
         * \param[in] P the new geometric location at the point
     */
    void set_point(index_t i, const vec3& P) {
        geo_debug_assert(new_vertices_.size() > 3 * i + 2);
        new_vertices_[3 * i] = P.x;
        new_vertices_[3 * i + 1] = P.y;
        new_vertices_[3 * i + 2] = P.z;
    }

    /**
         * \brief Sets a normal vector
         * \param[in] i the index of the point
         * \param[in] N the new normal vector associated with the point
     */
    void set_normal(index_t i, const vec3& N) {
        RVD_.set_normal(i, N);
    }

    /**
         * \brief Sets the maximum angle for determining admissible triangles.
         * \details Admissible triangles have a deviation between their normals
         *  and the normals estimated in the pointset smaller than a given
         *  threshold \p alpha.
         * \param[in] alpha the maximum normal angle deviation
     */
    void set_max_angle(double alpha) {
        min_cos_angle_ = ::cos(alpha);
    }

    Mesh& mesh() {
        return mesh_;
    }

  private:
    Mesh& mesh_;
    vector<double> new_vertices_;
    Co3NeRestrictedVoronoiDiagram RVD_;
    TypedThreadGroup<Co3neThread> thread_;
    double min_cos_angle_;
};

}
#endif // LPCVT_RECON_CO3NET_H
