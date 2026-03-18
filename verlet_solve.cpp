#include "verlet_solve.h"
#include "potential.h"
#include "consts.h"
#include "export_xyz.h"
#include <vector>
#include <list>
#include <iostream>
#include <fstream>

// tau/2
void L_xi(vector<double>& Xi, const vector<double>& V_xi, double tau)
{
    for (int k = 0; k < M; k++)
    {
        Xi[k] += V_xi[k]*tau;
    }
}

void L_v(vector<Vector3d>& V, const vector<double>& V_xi, const double& v_eps, double tau)
{
    const double k = V_xi[0] + v_eps*(1+(3./Nf));
    for (int i = 0; i < N; i++)
    {
        V[i] *= exp(-k*tau);
    }
}

// tau/4
void L_G1(const vector<Vector3d>& V, vector<double>& V_xi, const double& v_eps, double tau)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum += m*V[i].squaredNorm();
    }
    const double G1 = (sum + W*v_eps*v_eps - (Nf+1)*kb*Text)/Q[0];
    V_xi[0] += G1*tau;
}

void L_Gk(int k, vector<double>& V_xi, double tau) // k > 1
{
    const double Gk = (Q[k-1]*V_xi[k-1]*V_xi[k-1] - kb*Text)/Q[k];
    V_xi[k] += Gk*tau;
}

// tau/8
void L_v_xik(int k, vector<double>& V_xi, double tau)
{
    V_xi[k-1] *= exp( -V_xi[k]*tau );
}

void L_G_eps(const vector<Vector3d>& V, const double& eps,
             double& v_eps, double Wi, double tau)
{
    const double V0 = exp(3*eps);
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum += (1+(3./Nf))*m*V[i].squaredNorm();
    }
    double Geps = (1./W)*(sum + Wi - 3*V0*Pext);
    v_eps += Geps*tau;
}

void L_v_eps(const vector<double>& V_xi, double& v_eps, double tau)
{
    v_eps *= exp(-V_xi[0]*tau);
}

void L_NHCP(vector<Vector3d>& V,
           vector<double>& Xi, vector<double>& V_xi,
           double& eps, double& v_eps, double Wi, double tau)
{
    L_Gk(M-1, V_xi, tau/4);
    for (int k = M-1; k > 1; k--)
    {
        L_v_xik(k, V_xi, tau/8);
        L_Gk(k-1, V_xi, tau/4);
        L_v_xik(k, V_xi, tau/8);
    }
    L_v_xik(1, V_xi, tau/8);
    L_G1(V, V_xi, v_eps, tau/4);
    L_v_xik(1, V_xi, tau/8);

    L_v_eps(V_xi, v_eps, tau/8);
    L_G_eps(V, eps, v_eps, Wi, tau/4);
    L_v_eps(V_xi, v_eps, tau/8);

    L_v(V, V_xi, v_eps, tau/2);
    L_xi(Xi, V_xi, tau/2);

    L_v_eps(V_xi, v_eps, tau/8);
    L_G_eps(V, eps, v_eps, Wi, tau/4);
    L_v_eps(V_xi, v_eps, tau/8);

    L_v_xik(1, V_xi, tau/8);
    L_G1(V, V_xi, v_eps, tau/4);
    L_v_xik(1, V_xi, tau/8);
    for (int k = M-1; k > 1; k--)
    {
        L_v_xik(k, V_xi, tau/8);
        L_Gk(k-1, V_xi, tau/4);
        L_v_xik(k, V_xi, tau/8);
    }
    L_Gk(M-1, V_xi, tau/4);

    /*L_G2(V_xi, tau/4);
    L_v_xi1(V_xi, tau/8);
    L_G1(V, V_xi, v_eps, tau/4);
    L_v_xi1(V_xi, tau/8);

    L_v_eps(V_xi, v_eps, tau/8);
    L_G_eps(r, V, F, eps, v_eps, Wi, tau/4);
    L_v_eps(V_xi, v_eps, tau/8);

    L_v(V, V_xi, v_eps, tau/2);
    L_xi(Xi, V_xi, tau/2);

    L_v_eps(V_xi, v_eps, tau/8);
    L_G_eps(r, V, F, eps, v_eps, Wi, tau/4);
    L_v_eps(V_xi, v_eps, tau/8);

    L_v_xi1(V_xi, tau/8);
    L_G1(V, V_xi, v_eps, tau/4);
    L_v_xi1(V_xi, tau/8);
    L_G2(V_xi, tau/4);*/

}

void L1(vector<Vector3d>& V, const vector<Vector3d>& F, double tau)
{
    for (int i = 0; i < N; i++)
        V[i] += F[i]*tau/m;
}

void L2(vector<Vector3d>& R, const vector<Vector3d>& V,
        double& eps, const double& v_eps, Vector3d& L, double tau)
{
    const double s = std::exp(v_eps*tau);
    const double f = (std::abs(v_eps) > 1e-12) ? (s - 1.0)/v_eps : tau;
    for (int i = 0; i < N; i++)
    {
        R[i] = s*R[i] + f*V[i];
    }
    L *= s;

    eps += v_eps*tau;
}

void nullifyVcm(vector<Vector3d>& V)
{
    Vector3d vcm = Vector3d::Zero();
    const int N = V.size();
    for (int i = 0; i < N; ++i)
        vcm += m * V[i];
    vcm /= (m * N);
    for (int i = 0; i < N; ++i)
        V[i] -= vcm;
}

void VerletSolve(double t1, double t2, int niter,
            const vector<Vector3d>& R0,
            const vector<Vector3d>& V0,
            Vector3d& L,
            gsl_interp_accel* den_acc, gsl_spline* den_spl,
            gsl_interp_accel* embed_acc, gsl_spline* embed_spl,
            gsl_interp_accel* pair_acc, gsl_spline* pair_spl,
            const string& respath)
{
    std::ofstream xyz_stream(respath);    // Result file
    int N = R0.size();
    
    vector<Vector3d> R(N), V(N), F(N);
    vector<double> Xi(M, 0.), V_xi(M, 0.);
    double eps(0.), v_eps(0.);
    double Ve = L(0)*L(1)*L(2);
    eps = log(Ve) / 3.;

    double tau = (t2-t1) / static_cast<double>(niter-1);
    double time = 0;
    double K(0), T, Wi(0), P, sumT(0), sumP(0), H(0);

    vector<list<int>> neighbor_list(N);
    ConstructNeighborList(neighbor_list, R0, L);
    
    R = R0;
    V = V0;
    calcRho(R, L, neighbor_list, den_acc, den_spl, embed_acc, embed_spl);
    for (int k = 0; k < N; k++)
    {
        F[k] = calcdEtot(k, R, L, neighbor_list, den_acc, den_spl,
                           embed_acc, embed_spl, pair_acc, pair_spl);
    }
    
    export_xyz(xyz_stream, time, R, V);
    for (int t = 0; t < niter-1; t++)
    {
        time = (t != niter-2) ? time+tau : t2;

        // kinetic energy and virial
        for (int k = 0; k < N; k++)
        {
            vector<Vector3d> Fkj;
            calcFij(k, neighbor_list, Fkj);
            K += 0.5*m*V[k].squaredNorm();
            auto nblistiter = neighbor_list[k].cbegin();
            for (int j = 0; j < (int)Fkj.size(); j++)
            {
                Wi += 0.5*mic_disp(R[k]-R[*nblistiter++], L).dot(Fkj[j]);
            }

        }
        T = 2*K / (kb*Nf);
        H = K + calcEtot(R, L, neighbor_list, den_acc, den_spl, embed_acc, embed_spl, pair_acc, pair_spl)
                + 0.5*W*v_eps*v_eps + 0.5*Q[0]*V_xi[0]*V_xi[0] +
                 (Nf+1)*kb*T*Xi[0] + Pext*exp(3*eps);
        for (int k = 1; k < M; k++)
            H +=  0.5*Q[k]*V_xi[k]*V_xi[k] + kb*T*Xi[k];
        P = (N*kb*T + (1./3.)*Wi) / exp(3*eps);
        sumT += T;
        sumP += P;
        printf("%d / %d\n", t, niter);
        std::cout << T << std::endl;
        std::cout << P << "\t" << exp(3*eps) << "\t" << H << std::endl << std::endl;

        // integration step
        L_NHCP(V, Xi, V_xi, eps, v_eps, Wi, tau);
        L1(V, F, tau/2);
        L2(R, V, eps, v_eps, L, tau);

        if (isPeriodicBoundary)
        {
            for (int k = 0; k < N; k++)
            {
                for (int i = 0; i < 3; i++)
                {
                    if (R[k][i] >= L(i)/2)
                        R[k][i] -= L(i);
                    else if (R[k][i] < -L(i)/2)
                        R[k][i] += L(i);
                }
            }
        }
        calcRho(R, L, neighbor_list, den_acc, den_spl, embed_acc, embed_spl);
        for (int k = 0; k < N; k++)
        {
            F[k] = calcdEtot(k, R, L, neighbor_list, den_acc, den_spl,
                       embed_acc, embed_spl, pair_acc, pair_spl);
        }
        Wi = 0;
        for (int k = 0; k < N; k++)
        {
            vector<Vector3d> Fkj;
            calcFij(k, neighbor_list, Fkj);
            auto nblistiter = neighbor_list[k].cbegin();
            for (int j = 0; j < (int)Fkj.size(); j++)
            {
                Wi += 0.5*mic_disp(R[k]-R[*nblistiter++], L).dot(Fkj[j]);
            }

        }
        L1(V, F, tau/2);
        L_NHCP(V, Xi, V_xi, eps, v_eps, Wi, tau);


        export_xyz(xyz_stream, time, R, V);
        if ((t+1) % Nm == 0)
        {
            ConstructNeighborList(neighbor_list, R, L);
        }
        K = 0;
        Wi = 0;
    }

    printf("Итого:\n");
    printf("\t <T>=%.4f\n", sumT/niter);
    printf("\t <P>=%.9f\n", sumP/niter);
}

void ConstructNeighborList(vector<list<int>>& neighbor_list, const vector<Vector3d>& R, const Vector3d& L)
{
    int N = R.size();
    for (int i = 0; i < N; i++)
    {
        neighbor_list[i].clear();
        for (int j = 0; j < N; j++)
        {
            if ( (i != j)  && (calcrij(i, j, R, L) < rnei) )
                neighbor_list[i].push_back(j);
        }
    }
}

