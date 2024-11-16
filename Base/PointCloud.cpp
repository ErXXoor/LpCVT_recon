//
// Created by hongbo on 15/11/24.
//

#include "Base/PointCloud.h"
#include <fstream>
#include <iostream>
namespace NMCVT_PC{

size_t PointCloud::size(){
    return m_points.rows();
}

bool PointCloud::LoadFromXYZ(const std::string &filename) {
    std::ifstream infile(filename.c_str());
    if(!infile.is_open()){
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return false;
    }

    std::vector<Eigen::VectorXd> points;
    std::string line;
    size_t num_columns = 0;
    while(std::getline(infile, line)){
        std::istringstream iss(line);
        std::vector<double> vals;
        double val;
        while(iss >> val){
            vals.push_back(val);
        }

        if (num_columns == 0) {
            num_columns = vals.size();
            if (num_columns != 3 && num_columns != 6) {
                throw std::runtime_error("File format not supported: each line must have 3 (x, y, z) or 6 (x, y, z, nx, ny, nz) values.");
            }
        }

        // Check for consistent formatting across lines
        if (vals.size() != num_columns) {
            throw std::runtime_error("Inconsistent line format in file: each line must have the same number of values.");
        }

        points.emplace_back(Eigen::VectorXd::Map(vals.data(), vals.size()));
    }
    infile.close();

    m_points.resize(points.size(),num_columns);
    for(size_t i = 0; i < points.size(); i++){
        m_points.row(i) = points[i];
    }
    return true;
}
}
