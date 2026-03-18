#include "eam.h"
#include "verlet_solve.h"
#include "potential.h"
#include "consts.h"
#include <Eigen/Dense>
#include <gsl/gsl_spline.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <cmath>

void import_tablefun(const string& path, gsl_interp_accel **acc, gsl_spline **spl) {
    std::ifstream file(path);
    string line;
    
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
    }
    
    auto strarr = split_view(line, ' ');
    int N = std::stoi(string(strarr[1]));
    
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
    }
    
    vector<double> X(N), Y(N);
    
    for (int i = 0; i < N; ++i) {
        std::getline(file, line);
        auto strarr = split_view(line, ' ');
        X[i] = std::stod(string(strarr[0]));
        Y[i] = std::stod(string(strarr[1]));
    }
    
    
    *acc = gsl_interp_accel_alloc();
    *spl = gsl_spline_alloc(gsl_interp_cspline, X.size());
    gsl_spline_init(*spl, X.data(), Y.data(), X.size());
    
}


void getCellsCoords(double a, double b, double c,
                           double alpha_deg, double beta_deg, double gamma_deg,
                           char cent, int n,
                           vector<Vector3d> &Atoms, int &vertices_count)
{
    Atoms.clear();
    const double ca = std::cos(deg2rad(alpha_deg));
    const double sa = std::sin(deg2rad(alpha_deg));
    const double cb = std::cos(deg2rad(beta_deg));
    const double cg = std::cos(deg2rad(gamma_deg));

    // Basis vectors
    Vector3d v1(a, 0.0, 0.0);
    Vector3d v2(b * ca, b * sa, 0.0);
    Vector3d v3;
    v3.x() = c * cb;
    v3.y() = c * (cg - ca * cb) / sa;
    const double v3xy2 = v3.x() * v3.x() + v3.y() * v3.y();
    const double z2 = c * c - v3xy2;
    v3.z() = std::sqrt(std::max(0.0, z2));

    std::vector<double> grid(n);
    const double start = -static_cast<double>(n) / 2.0;
    for (int i = 1; i < n+1; ++i)
        grid[i-1] = start + i;

    Atoms.reserve(static_cast<std::size_t>(n) * (n) * (n) * (cent == 'F' ? 4 : 1));

    for (double ia : grid)
        for (double jb : grid)
            for (double kc : grid)
                Atoms.emplace_back(ia * v1 + jb * v2 + kc * v3);

    vertices_count = Atoms.size();

    if (cent == 'F') //&& n >= 2)
    {
        std::vector<double> mid;
        mid.reserve(n);
        for (int i = 0; i < n; ++i)
            mid.push_back(grid[i] + 0.5);

        for (double ia : grid)
            for (double jb : mid)
                for (double kc : mid)
                    Atoms.emplace_back(ia * v1 + jb * v2 + kc * v3);

        for (double ia : mid)
            for (double jb : grid)
                for (double kc : mid)
                    Atoms.emplace_back(ia * v1 + jb * v2 + kc * v3);

        for (double ia : mid)
            for (double jb : mid)
                for (double kc : grid)
                    Atoms.emplace_back(ia * v1 + jb * v2 + kc * v3);
    }
}

void getCellsCoords_nonper(double a, double b, double c,
                           double alpha_deg, double beta_deg, double gamma_deg,
                           char cent, int n,
                           vector<Vector3d> &Atoms, int &vertices_count)
{
    Atoms.clear();
    const double ca = std::cos(deg2rad(alpha_deg));
    const double sa = std::sin(deg2rad(alpha_deg));
    const double cb = std::cos(deg2rad(beta_deg));
    const double cg = std::cos(deg2rad(gamma_deg));

    Vector3d v1(a, 0.0, 0.0);
    Vector3d v2(b * ca, b * sa, 0.0);

    Vector3d v3;
    v3.x() = c * cb;
    v3.y() = c * (cg - ca * cb) / sa;
    const double v3xy2 = v3.x() * v3.x() + v3.y() * v3.y();
    const double z2 = c * c - v3xy2;
    v3.z() = std::sqrt(std::max(0.0, z2));

    std::vector<double> grid(n+1);
    const double start = -static_cast<double>(n) / 2.0;
    for (int i = 0; i < n+1; ++i)
        grid[i] = start + i;

    Atoms.reserve(static_cast<std::size_t>(n+1) * (n+1) * (n+1) * (cent == 'F' ? 4 : 1));

    for (double ia : grid)
        for (double jb : grid)
            for (double kc : grid)
                Atoms.emplace_back(ia * v1 + jb * v2 + kc * v3);

    vertices_count = Atoms.size();

    if (cent == 'F') //&& n >= 2)
    {
        std::vector<double> mid;
        mid.reserve(n);
        for (int i = 0; i < n; ++i)
            mid.push_back(grid[i] + 0.5);

        for (double ia : grid)
            for (double jb : mid)
                for (double kc : mid)
                    Atoms.emplace_back(ia * v1 + jb * v2 + kc * v3);

        for (double ia : mid)
            for (double jb : grid)
                for (double kc : mid)
                    Atoms.emplace_back(ia * v1 + jb * v2 + kc * v3);

        for (double ia : mid)
            for (double jb : mid)
                for (double kc : grid)
                    Atoms.emplace_back(ia * v1 + jb * v2 + kc * v3);
    }
}

void do_ConstructNeighborList(vector<list<int>>& neighbor_list, const vector<Vector3d>& R, const Vector3d& L)
{
    int N = R.size();
    for (int i = 0; i < N; i++)
    {
        neighbor_list[i].clear();
        for (int j = 0; j < N; j++)
        {
            if ( (i != j) && (calcrij(i, j, R, L) < rnei) )
                neighbor_list[i].push_back(j);
        }
    }
}

double findLatticeConstant(double bst, double bend,
                           gsl_interp_accel* den_acc,   gsl_spline* den_spl,
                           gsl_interp_accel* embed_acc, gsl_spline* embed_spl,
                           gsl_interp_accel* pair_acc,  gsl_spline* pair_spl,
                           int n, char bravais = 'F')
{
    const double tau = 0.5*(std::sqrt(5.0)-1.0); // ≈ 0.618
    auto evalE = [&](double a)->double {
        std::vector<Eigen::Vector3d> R;
        int Ncells = 0;
        getCellsCoords(a, a, a, 90, 90, 90, bravais, n, R, Ncells);

        const double delta = 0;//1e-4;
        Eigen::Vector3d L(n*a+delta, n*a+delta, n*a+delta);

        for (int k = 0; k < (int)R.size(); k++)
        {
            R[k] += L/2;
        }


        // список соседей (симметричный) под PBC
        std::vector<std::list<int>> nbr(R.size());
        do_ConstructNeighborList(nbr, R, L);

        double E = calcEtot(R, L, nbr, den_acc, den_spl,
                            embed_acc, embed_spl, pair_acc, pair_spl);
        return E / R.size(); // total energy per atom
    };

    // initial points of the golden ratio
    double b1 = bend - tau*(bend - bst);
    double b2 = bst  + tau*(bend - bst);
    double E1 = evalE(b1);
    double E2 = evalE(b2);

    while ((bend - bst) > 1e-7) {
        if (E1 > E2) {
            bst = b1;
            b1  = b2;
            E1  = E2;
            b2  = bst + tau*(bend - bst);
            E2  = evalE(b2);
        } else {
            bend = b2;
            b2   = b1;
            E2   = E1;
            b1   = bend - tau*(bend - bst);
            E1   = evalE(b1);
        }
    }
    return 0.5*(bst + bend);
}



int main()
{
    gsl_interp_accel *pair_acc, *embed_acc, *den_acc;
    gsl_spline *pair_spl, *embed_spl, *den_spl;
    
    import_tablefun("/home/georgiy/Desktop/mishin1999(al)/Nose-Hoover/al_den_eam_5_1_mishin_1999_al.txt", &den_acc, &den_spl);
    import_tablefun("/home/georgiy/Desktop/mishin1999(al)/Nose-Hoover/al_embed_eam_5_1_mishin_1999_al.txt", &embed_acc, &embed_spl);
    import_tablefun("/home/georgiy/Desktop/mishin1999(al)/Nose-Hoover/al_pair_eam_5_1_mishin_1999_al.txt", &pair_acc, &pair_spl);

    /*double ares = findLatticeConstant(4, 4.5, den_acc, den_spl,
                                      embed_acc, embed_spl, pair_acc, pair_spl, n);
    printf("%.7f\n", ares);
    return 1;
    omp_set_num_threads(1);*/
    
    vector<Vector3d> r0;
    int l;
    Vector3d L;
    const double delta = 0;//a/2; //a/2;//0.005;
    if (isPeriodicBoundary)
        getCellsCoords(a, a, a, 90, 90, 90, 'F', n, r0, l);
    else
        getCellsCoords_nonper(a, a, a, 90, 90, 90, 'F', n, r0, l);

    std::cout << r0.size() << std::endl;
    L(0) = n*a; L(1) = n*a; L(2) = n*a;
    for (int k = 0; k < N; k++)
    {
        r0[k] += L/2;
    }
    L(0) += delta;
    L(1) += delta;
    L(2) += delta;

    if (rnei >= 0.5*L(0))
    {
        std::cout << "The box size is too small" << std::endl;
        return 1;
    }
    
    // generation of initial velocities (Boltzmann distribution)
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::normal_distribution<> boltz_dist(0., sqrt(kb*Text/m));
    vector<Vector3d> v0(N);
    Vector3d vcm; vcm.setZero();
    for (int k = 0; k < N; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            v0[k](i) = boltz_dist(gen);
            vcm(i) += v0[k][i];
        }
    }
    vcm = (1./N)*vcm;
    // remove the movement as a whole
    for (int k = 0; k < N; k++)
    {
        v0[k] = v0[k] - vcm;
    }
        
    double k = 10;
    double t1 = 0;
    double t2 = FromPicoSec*0.8*k;
    int niter = k*400;
    if (isPeriodicBoundary)
    {
        for (int k = 0; k < N; k++)
        {
            for (int i = 0; i < 3; i++)
            {
                if (r0[k][i] >= L(i)/2)
                    r0[k][i] -= L(i);
                else if (r0[k][i] < -L(i)/2)
                    r0[k][i] += L(i);
            }
        }
    }
    
    VerletSolve(t1, t2, niter, r0, v0,
                L,  den_acc, den_spl,
                embed_acc, embed_spl, pair_acc, pair_spl,
                "./res.xyz");
    
}
            
