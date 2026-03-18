#pragma once
#include <string_view>
#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <gsl/gsl_spline.h>

using std::string;
using std::vector;
using Eigen::Vector3d;
using Eigen::VectorXd;

inline double deg2rad(double deg) { return deg * M_PI / 180.; }

vector<std::string_view> split_view(std::string_view s, char delim, bool keep_empty=false) {
    vector<std::string_view> out;
    size_t start = 0;
    while (true) {
        size_t pos = s.find(delim, start);
        std::string_view token = (pos == std::string_view::npos)
                               ? s.substr(start)
                               : s.substr(start, pos - start);
        if (keep_empty || !token.empty()) out.push_back(token);
        if (pos == std::string_view::npos) break;
        start = pos + 1;
    }
    return out;
};

// Read a file with a table and build a cubic spline from it
// Input: path - file name
// Output: acc, spl - spline
void import_tablefun(const string& path, gsl_interp_accel **acc, gsl_spline **spl);

// Построить сетку из атомов
void getCellsCoords(double a, double b, double c,
                           double alpha_deg, double beta_deg, double gamma_deg,
                           char cent, int n, 
                           vector<Vector3d> &Atoms, int &vertices_count);

/*double findLatticeConstant(double b1, double b2,
                           gsl_interp_accel* den_acc, gsl_spline* den_spl,
                           gsl_interp_accel* embed_acc, gsl_spline* embed_spl,
                           gsl_interp_accel* pair_acc, gsl_spline* pair_spl);*/
