#pragma once
#include <vector>
#include <Eigen/Dense>
#include <fstream>

using std::vector;
using Eigen::Vector3d;

void export_xyz(std::ofstream& xyz_stream,  double t,
    const vector<Vector3d>& R,
    const vector<Vector3d>& V);
