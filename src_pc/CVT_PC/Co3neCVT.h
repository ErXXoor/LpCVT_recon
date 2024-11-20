//
// Created by hongbo on 16/11/24.
//

#ifndef LPCVT_RECON_CO3NECVT_H
#define LPCVT_RECON_CO3NECVT_H
#include "Co3neT.h"
#include "NMCVT/LpCVTWrap.h"
#include "NMCVT/LpCVTIS.h"
#include <geogram/basic/smart_pointer.h>
#include <memory>

namespace NMCVT_PC{
class Co3neCVT {
  public:
    Co3neCVT(std::shared_ptr<GEO::Mesh> M);
    void InitialSampling(int num_samples, const Eigen::MatrixXd &input_pc);
    void SurfaceRecon(double radius);
    void InitRemeshing(GEO::SmartPointer<LpCVT::LpCVTIS> is);

    void NMRemeshing(unsigned int nb_pts,
                     unsigned int nb_iter);

    void GetRDT(GEO::Mesh &M_out, bool post_process=true);

  private:
    std::shared_ptr<Co3neT> m_co3;
    std::shared_ptr<LpCVT::LpCVTWrap> m_cvt;
    std::shared_ptr<GEO::Mesh> m_mesh;
    GEO::IntegrationSimplex_var m_is_var;
    Eigen::MatrixXd m_sites;
};
}


#endif // LPCVT_RECON_CO3NECVT_H
