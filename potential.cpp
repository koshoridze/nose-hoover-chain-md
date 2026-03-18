#include "potential.h"
#include "consts.h"
#include <iostream>

Vector3d mic_disp(const Vector3d& dr, const Vector3d& L)
{
    if (isPeriodicBoundary)
    {
        Vector3d out = dr;
        /*for (int a = 0; a < 3; a++)
            out(a) -= L(a)*std::nearbyint(out(a)/L(a));
        return out;*/
        for (int i = 0; i < 3; i++)
        {
            if (out[i] >= L(i)/2)
                out[i] -= L(i);
            else if (out[i] < -L(i)/2)
                out[i] += L(i);
        }
        return out;
    }
    else
        return dr;
}

double calcrij(int i, int j, const vector<Vector3d>& Atoms, const Vector3d& L)
{
    Vector3d ri = Atoms[i];
    Vector3d rj = Atoms[j];
    Vector3d rij = mic_disp(rj - ri, L);
    return rij.norm();
}

double spline_eval(gsl_spline* spl, double x, gsl_interp_accel* acc, const double a, const double b)
{
    if (x < a)
        x = a;
    else if (x > b)
        x = b;
    return gsl_spline_eval(spl, x, acc);
}

double spline_eval_deriv(gsl_spline* spl, double x, gsl_interp_accel* acc, const double a, const double b)
{
    if (x < a)
        x = a;
    else if (x > b)
        x = b;
    return gsl_spline_eval_deriv(spl, x, acc);
}

double calcEtot(const vector<Vector3d>& Atoms, const Vector3d& L,
                const vector<list<int>>& neighbor_list,
                gsl_interp_accel* den_acc, gsl_spline* den_spl,
                gsl_interp_accel* embed_acc, gsl_spline* embed_spl,
                gsl_interp_accel* pair_acc, gsl_spline* pair_spl)
{
    double Etot = 0;
    for (int i = 0; i < N; i++)
    {
        double rho = 0;
        for (int j : neighbor_list[i])
        {
            double rij = calcrij(i, j, Atoms, L);
            Etot += 0.5*spline_eval(pair_spl, rij, pair_acc, a_pair, b_pair); // 0.5 ????????????????????????????????????
            rho += spline_eval(den_spl, rij, den_acc, a_den, b_den);
        }
        Etot += spline_eval(embed_spl, rho, embed_acc, a_embed, b_embed);
    }
    return Etot;
}

vector<list<Vector3d>> Fij_lists(N);
vector<double> Rho(N, 0.0);
vector<double> dF(N, 0.0);

void calcRho(const vector<Vector3d>& Atoms, const Vector3d& L,
             const vector<list<int>>& neighbor_list,
             gsl_interp_accel* den_acc,   gsl_spline* den_spl,
             gsl_interp_accel* embed_acc, gsl_spline* embed_spl)
{
    // Calc all \overline{rho}_i
    for (int i = 0; i < N; i++){
        double rho = 0;
        for (int j : neighbor_list[i])
        {
            double rij = calcrij(i, j, Atoms, L);
            rho += spline_eval(den_spl, rij, den_acc, a_den, b_den);
        }
        Rho[i] = rho;
        dF[i] = spline_eval_deriv(embed_spl, Rho[i], embed_acc, a_embed, b_embed);
    }
}

// Force on the k-th particle (EAM)
Vector3d calcdEtot(
    int k, const vector<Vector3d>& Atoms, const Vector3d& L,
    const vector<list<int>>& neighbor_list,
    gsl_interp_accel* den_acc,   gsl_spline* den_spl,   
    gsl_interp_accel* /*embed_acc*/, gsl_spline* /*embed_spl*/,
    gsl_interp_accel* pair_acc,  gsl_spline* pair_spl   
){
    Fij_lists[k].clear();    

    Vector3d F = Vector3d::Zero();
    for (int j : neighbor_list[k])
    {
        Vector3d ekj = mic_disp(Atoms[k]-Atoms[j], L);
        double rkj = ekj.norm();
        ekj.normalize();

        double dphi = spline_eval_deriv(pair_spl, rkj, pair_acc, a_pair, b_pair); 
        double df   = spline_eval_deriv(den_spl,  rkj, den_acc, a_den, b_den);  

        Vector3d Fkj = -(dphi + (dF[k] + dF[j]) * df) * ekj;
        Fij_lists[k].push_back(Fkj);
        F += Fkj;
    }
    return F;
}

void calcFij(int k, const vector<list<int>>& neighbor_list, vector<Vector3d>& Fij)
{
    /*Vector3d sum = Vector3d::Zero();
    for (const auto& lst : Fij_lists)
      for (const auto& f : lst)
        sum += f;

    if (sum.norm() > 1e-10) {
      std::cerr << "Sum of internal forces != 0, norm = " << sum.norm() << "\n";
    }*/

    Fij.resize(neighbor_list[k].size());
    int ind = 0;
    auto iter = Fij_lists[k].begin();
    for (int j : neighbor_list[k])
    {
        Fij[ind++] = *iter++;
    }

}
