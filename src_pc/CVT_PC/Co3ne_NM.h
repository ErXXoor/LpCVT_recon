//
// Created by hongbo on 15/11/24.
//

#ifndef LPCVT_RECON_CO3NE_NM_H
#define LPCVT_RECON_CO3NE_NM_H
#include <geogram/points/co3ne.h>
#include <geogram/points/kd_tree.h>
#include <geogram/mesh/mesh.h>
#include <Eigen/Dense>
namespace NMCVT_PC{
    class Co3ne_NM{
      public:
        Co3ne_NM() = default;
        ~Co3ne_NM() = default;

        void Reconstruct(Eigen::MatrixXd &points, GEO::Mesh &M);
        void Smooth(GEO::Mesh &M,int nb_neighbors,int nb_iterations);
    };
}


#endif // LPCVT_RECON_CO3NE_NM_H
