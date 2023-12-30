#include <iostream>
#include <geogram/mesh/mesh.h>
#include "geogram/basic/line_stream.h"
#include "Base/Mesh.h"

int main() {
    GEO::Mesh M_in;
    std::string input_filename = "/Users/lihongbo/Desktop/code/LpCVT_recon/tmp/datatech_1_input.obj";

    std::shared_ptr<LpCVT::Mesh> mesh = std::make_shared<LpCVT::Mesh>();
    mesh->LoadMesh(input_filename);
    mesh->CalculateFaceNormal();
    mesh->PlotMesh();
}
