//
// Created by hongbo on 15/11/24.
//

#ifndef LPCVT_RECON_CO3NETHREAD_H
#define LPCVT_RECON_CO3NETHREAD_H
#include <geogram/basic/process.h>
#include <geogram/points/principal_axes.h>
#include <geogram/basic/numeric.h>
namespace NMCVT_PC{
using namespace GEO;
class Co3neT;
class Co3neThread : public Thread {
  public:
    /**
         * \brief Creates a new Co3NeThread
         * \param[in] master the Co3Ne this thread depends on
         * \param[in] from index of the first point to process
         * \param[in] to one position past the index of the last point
     */
    Co3neThread(Co3neT* master,index_t from, index_t to);

    /**
         * \brief Sets the mode of this thread
         * \param[in] m the mode, that determines whether normal computation,
         *  smoothing or reconstruction is performed
     */
//    void set_mode(Co3neT::Co3NeMode m) {
//        mode_ = m;
//    }

    /**
         * \brief Does the actual computation of this thread.
         * \details The actual computation is determined by set_mode().
     */
    void run() override {
//        switch(mode_) {
//        case CO3NE_NORMALS:
//            run_normals();
//            break;
//        case CO3NE_SMOOTH:
//            run_smooth();
//            break;
//        case CO3NE_RECONSTRUCT:
//            run_reconstruct();
//            break;
//        case CO3NE_NORMALS_AND_RECONSTRUCT:
            run_normals_and_reconstruct();
//            break;
//        case CO3NE_NONE:
//            break;
//        }
    }

    /**
         * \brief Gets the reconstructed triangles.
         * \return a reference to a vector of indices
     */
    vector<index_t>& triangles() {
        return triangles_;
    }


    /**
         * \brief Gets the number of reconstructed triangles.
         * \return the number of reconstructed triangles
     */
    index_t nb_triangles() const {
        return triangles_.size()/3;
    }

  protected:
    /**
         * \brief Estimates the normals in the pointset.
     */
    void run_normals();

    /**
         * \brief Smoothes the pointset.
     */
    void run_smooth();

    /**
         * \brief Reconstructs the triangles.
     */
    void run_reconstruct();

    /**
         * \brief Estimates the normals and reconstructs the triangles.
     */
    void run_normals_and_reconstruct();

  private:
    Co3neT* master_;
    index_t from_;
    index_t to_;
//    Co3NeMode mode_;
    PrincipalAxes3d least_squares_normal_;
    vector<index_t> triangles_;
};
}
#endif // LPCVT_RECON_CO3NETHREAD_H
