#pragma once

static inline constexpr double sqr(double x)  { return x * x; }

constexpr bool isPeriodicBoundary = true;

constexpr int n = 4; // Supercell size
constexpr int N = 256; // Number of atoms
constexpr int Nf = 3*N-3; // degrees of freedom

constexpr double a = 4.04263158; //(Angstrom)
constexpr double m = 26.98; // mass, atomic mass unit

// Граничные аргументы потенцилов
constexpr double a_den = 1.0017625; constexpr double b_den = 6.28720;
constexpr double a_pair = a_den; constexpr double b_pair = b_den;
constexpr double a_embed = 0.001; constexpr double b_embed = 2.;

constexpr double rcut = b_den; // potential cutoff distance
constexpr double rnei = 1.2*rcut; // neighbor-list distance
constexpr int Nm = 1; // neighbor list update frequency

constexpr double FromPicoSec = 98.226947433914546; // convert ps to inner time format

constexpr double kb = 8.617333*1e-5; // Bolzmann constant

constexpr double Text = 270; // K
constexpr int M = 2; // length of thermostat chain

constexpr double tau1_ps = 100*0.002;  // thermostat relaxation time in ps
constexpr double Q1 = 9*kb*Text*sqr(FromPicoSec*tau1_ps);
constexpr double Q2 = kb*Text*sqr(FromPicoSec*tau1_ps);

constexpr double Q[M] = {Q1, Q2};//, Q2, Q2, Q2};

// A small value of Q will lead to rapid fluctuations of the kinetic energy
// too large value of Q will greatly inhibit the energy exchange between
// the initial system and the thermostat, thereby leading to poor temperature control.
// Reccomended: tau_Q = 100*tau

constexpr double tau2_ps = 7*tau1_ps;
constexpr double Pext = 0;
constexpr double W = Nf*kb*Text*sqr(FromPicoSec*tau2_ps);
//constexpr double W = 3*365*kb*T*100*100*5*5/10;













constexpr double epslj = 4e-3;
constexpr double sigmalj = 0.2236;

constexpr double q = 2;
constexpr double kc = 5;
