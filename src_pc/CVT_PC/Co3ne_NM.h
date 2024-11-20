//
// Created by hongbo on 15/11/24.
//

#ifndef LPCVT_RECON_CO3NE_NM_H
#define LPCVT_RECON_CO3NE_NM_H
#include <geogram/points/co3ne.h>
#include <geogram/mesh/mesh.h>
#include "Co3neCVT.h"
#include "NMCVT/LpCVTIS.h"
#include <geogram/basic/smart_pointer.h>
namespace NMCVT_PC{
    class Co3ne_NM{
      public:
        Co3ne_NM() = default;
        ~Co3ne_NM() = default;

        void Reconstruct(Eigen::MatrixXd &points,
                         unsigned int nb_pts,
                         unsigned int nb_iter);
        void Smooth(GEO::Mesh &M,int nb_neighbors,int nb_iterations);
        void GetRDT(GEO::Mesh &M_out, bool post_process=true);
        void GetReconstructedMesh(LpCVT::Mesh &M_out);
      private:
        std::shared_ptr<GEO::Mesh> m_mesh;
        std::shared_ptr<Co3neCVT> m_cvt_core;
        GEO::SmartPointer<LpCVT::LpCVTIS> m_is;

      public:

    };
}


#endif // LPCVT_RECON_CO3NE_NM_H
