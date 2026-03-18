#pragma once
#include <vector>
#include <string>
#include <list>
#include <Eigen/Dense>
#include <gsl/gsl_spline.h>

using std::vector;
using std::string;
using std::list;
using Eigen::Vector3d;

void VerletSolve(double t1, double t2, int niter,
            const vector<Vector3d>& r0,
            const vector<Vector3d>& v0,
            //vector<double>& Time,
           // vector<vector<Vector3d>>& R,
            //vector<vector<Vector3d>>& V,
            Vector3d& L,
            gsl_interp_accel* den_acc, gsl_spline* den_spl,
            gsl_interp_accel* embed_acc, gsl_spline* embed_spl,
            gsl_interp_accel* pair_acc, gsl_spline* pair_spl,
            const string& respath);

void ConstructNeighborList(vector<list<int>>& neighbor_list, const vector<Vector3d>& R, const Vector3d& L);
