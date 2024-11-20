//
// Created by hongbo on 16/11/24.
//
#include "CVT_PC/Co3neCVT.h"
#include <pcl/filters/random_sample.h>
#include <pcl/point_types.h>

namespace NMCVT_PC{
Co3neCVT::Co3neCVT(std::shared_ptr<GEO::Mesh> M) {
    m_mesh = M;
    m_co3 = std::make_shared<Co3neT>(*m_mesh);
}

void Co3neCVT::SurfaceRecon(double radius) {
    m_co3->reconstruct(radius);
}

void Co3neCVT::InitRemeshing(GEO::SmartPointer<LpCVT::LpCVTIS> is){
    m_cvt = std::make_shared<LpCVT::LpCVTWrap>(m_mesh.get(),3);
    m_cvt->set_volumetric(false);
    is->set_Delaunay(m_cvt->delaunay());
    m_is_var = is;
    m_cvt->set_simplex_func(m_is_var);
}

void Co3neCVT::NMRemeshing(unsigned int nb_pts, unsigned int nb_iter) {
    m_cvt->compute_initial_sampling(nb_pts);
    GEO::Logger::div("Optimize sampling");
    GEO::ProgressTask progress("Optimize sampling", nb_iter);
    m_cvt->set_progress_logger(&progress);

//    if (m_type == RemeshType::Lloyd_CVT) {
        m_cvt->Lloyd_iterations(nb_iter);
//    } else {
//        m_cvt->Newton_iterations(nb_iter);
//    }

    m_cvt->set_progress_logger(nullptr);
}

void Co3neCVT::GetRDT(GEO::Mesh &M_out, bool post_process) {
    GEO::Logger::div("Generate RDT");

    if (post_process) {
        m_cvt->compute_surface(&M_out, true);
    } else {
        m_cvt->RVD()->set_exact_predicates(true);
        auto mode_i = GEO::RestrictedVoronoiDiagram::RDT_RVC_CENTROIDS |
                      GEO::RestrictedVoronoiDiagram::RDT_PREFER_SEEDS;
        auto mode = GEO::RestrictedVoronoiDiagram::RDTMode(mode_i);
        m_cvt->RVD()->compute_RDT(M_out, mode);
    }

}


void Co3neCVT::InitialSampling(int num_samples, const Eigen::MatrixXd &input_pc) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    cloud->width = input_pc.rows();
    cloud->height = 1;
    cloud->is_dense = true;
    cloud->points.resize(cloud->width * cloud->height);

    for(int i = 0; i < input_pc.rows(); i++) {
        cloud->points[i].x = input_pc(i, 0);
        cloud->points[i].y = input_pc(i, 1);
        cloud->points[i].z = input_pc(i, 2);
    }

    pcl::RandomSample<pcl::PointXYZ> randomSample;
    randomSample.setInputCloud(cloud);
    randomSample.setSample(num_samples);
    pcl::PointCloud<pcl::PointXYZ>::Ptr sampled_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    randomSample.filter(*cloud);

    m_sites.resize(sampled_cloud->size(), 3);
    for(int i = 0; i < sampled_cloud->size(); i++) {
        m_sites(i, 0) = sampled_cloud->points[i].x;
        m_sites(i, 1) = sampled_cloud->points[i].y;
        m_sites(i, 2) = sampled_cloud->points[i].z;
    }
}
}