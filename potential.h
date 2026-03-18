#pragma once
#include <Eigen/Dense>
#include <vector>
#include <list>
#include <gsl/gsl_spline.h>

using Eigen::Vector3d;
using std::vector;
using std::list;

Vector3d mic_disp(const Vector3d& dr, const Vector3d& L);

double calcrij(int i, int j, const vector<Vector3d>& Atoms, const Vector3d& L);

double calcEtot(const vector<Vector3d>& Atoms, const Vector3d& L,
                const vector<list<int>>& neighbor_list,
                gsl_interp_accel* den_acc, gsl_spline* den_spl,
                gsl_interp_accel* embed_acc, gsl_spline* embed_spl,
                gsl_interp_accel* pair_acc, gsl_spline* pair_spl);

void calcRho(const vector<Vector3d>& Atoms, const Vector3d& L,
             const vector<list<int>>& neighbor_list,
             gsl_interp_accel* den_acc,   gsl_spline* den_spl,
             gsl_interp_accel* embed_acc, gsl_spline* embed_spl);

// Force acting on the k-th particle
Vector3d calcdEtot(int k, const vector<Vector3d>& Atoms, const Vector3d& L,
                const vector<list<int>>& neighbor_list,
                gsl_interp_accel* den_acc, gsl_spline* den_spl,
                gsl_interp_accel* embed_acc, gsl_spline* embed_spl,
                gsl_interp_accel* pair_acc, gsl_spline* pair_spl);

Vector3d calcLennardJones(int k, const vector<Vector3d>& Atoms, const Vector3d& L);

Vector3d calcCoulomb(int k, const vector<Vector3d>& Atoms, const Vector3d& L);

double spline_eval(gsl_spline* spl, double x, gsl_interp_accel* acc, const double a, const double b);

double spline_eval_deriv(gsl_spline* spl, double x, gsl_interp_accel* acc, const double a, const double b);

void calcFij(int k, const vector<list<int>>& neighbor_list, vector<Vector3d>& Fij);
