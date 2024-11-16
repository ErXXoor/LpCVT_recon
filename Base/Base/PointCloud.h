//
// Created by hongbo on 15/11/24.
//

#ifndef LPCVT_RECON_POINTS_H
#define LPCVT_RECON_POINTS_H

#endif // LPCVT_RECON_POINTS_H
#include <string>
#include <pcl/point_types.h>
namespace NMCVT_PC{
    class PointCloud{
      public:
        PointCloud() = default;
        ~PointCloud() = default;
        size_t size();

        bool LoadFromXYZ(const std::string &filename);
      private:
        Eigen::MatrixXd m_points;

      public:
        Eigen::MatrixXd GetPoints(){
            return m_points;
        }
    };

}