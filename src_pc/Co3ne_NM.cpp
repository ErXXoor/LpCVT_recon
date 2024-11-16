//
// Created by hongbo on 15/11/24.
//
#include "CVT_PC/Co3ne_NM.h"
#include "Base/MeshAdaptor.h"
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include "CVT_PC/Co3neT.h"
#include "CVT_PC/Co3neRVD.h"

namespace NMCVT_PC{

    void Co3ne_NM::Reconstruct(Eigen::MatrixXd &points, GEO::Mesh &M) {
        LpCVT::MeshAdaptor::Convert(points, M);
        double R = GEO::bbox_diagonal(M);
        GEO::mesh_repair(M, GEO::MESH_REPAIR_COLOCATE, 1e-6*R);

//        Smooth(M,)

        auto radius = 0.01*R*5.0;
//        GEO::Co3Ne_reconstruct(M, radius);
        Co3neT co3ne(M);
        co3ne.reconstruct(radius);
    }
    void Co3ne_NM::Smooth(Mesh& M, int nb_neighbors, int nb_iterations) {
        Co3neT co3ne(M);
        try {
            ProgressTask progress("Smoothing", nb_iterations);
            for(int i = 0; i < nb_iterations; i++) {
                co3ne.smooth(nb_neighbors);
                if(i != nb_iterations - 1) {
                    co3ne.RVD().update();
                }
                progress.next();
            }
            co3ne.end_smooth();
        }
        catch(const TaskCanceled&) {
        }
    }

}