#include "xlife++.h"
#include <fstream>
#include <iostream>
using namespace xlifepp;
using namespace std;

// Function to get the radius of the a circle
Real get_radius(Domain dom)
{
  Real x_max = 0;
  for (const auto &P : dom.meshDomain()->nodes())
  {
    x_max = std::max(x_max, P(1));
  }
  return x_max;
}

// Function d tilde PML
Complex dt(Real r, Real r0, Real r1, Real s0)
{
  return 1 + s0 * i_ * (r - r0) / (r1 - r0);
  // return 1.+s0*ic*log((r1-r0)/(r1-r))/r;
}

// Function d PML
Complex d_(Real r, Real r0, Real r1, Real s0)
{
  return 1 + s0 * i_ * (2 * r - r0) / (r1 - r0);
  // return 1.+s0*ic/(r1-r);
}

// Matrix for PML in stiffness matrix
Matrix<Complex> MatPML(const Point &P, Parameters &pa = defaultParameters)
{
  Real r0 = pa("r0"), r1 = pa("r1"), s0 = pa("s0");
  Matrix<Complex> M(2, _idMatrix);
  Real r2 = P(1) * P(1) + P(2) * P(2);
  Real r = sqrt(r2);
  Complex d0 = dt(r, r0, r1, s0);
  Complex d1 = d_(r, r0, r1, s0);
  M(1, 1) = (P(1) * P(1) * d0 / d1 + P(2) * P(2) * d1 / d0) / r2;
  M(2, 2) = (P(2) * P(2) * d0 / d1 + P(1) * P(1) * d1 / d0) / r2;
  M(1, 2) = P(1) * P(2) * (d0 / d1 - d1 / d0) / r2;
  M(2, 1) = M(1, 2);
  return M;
}

// Function for PML in mass matrix
Complex FctPML(const Point &P, Parameters &pa = defaultParameters)
{
  Real r0 = pa("r0"), r1 = pa("r1"), s0 = pa("s0");
  Real r = norm2(P);
  return dt(r, r0, r1, s0) * d_(r, r0, r1, s0);
}

int main(int argc, char **argv)
{
  init(_lang = en); // mandatory initialization of xlifepp
  numberOfThreads(1);

  Real eps = -1.1;

  Int nbvp = 2;

  //// Mesh and domains definition
  Mesh mesh = Mesh("../../GMSH/disk_T_PML.msh", "disk_T_PML", msh);
  Domain Omega = mesh.domain("Omega");
  Domain Cavity = mesh.domain("Cavity");
  Domain Vacuum = mesh.domain("Vacuum");
  Domain PML = mesh.domain("PML");
  Domain Gamma = mesh.domain("Gamma");
  Domain Sigma = mesh.domain("Sigma");
  Domain SigPML = mesh.domain("SigPML");

  //// PML Parameters
  Real r0 = get_radius(Sigma);
  Real r1 = get_radius(SigPML);

  Parameters pars;
  pars << Parameter(r0, "r0") << Parameter(r1, "r1");

  //// Space and unknown
  Space V(Omega, interpolation(Lagrange, _GaussLobattoPoints, 8, H1), "V");
  Unknown u(V, "u");
  TestFunction v(u, "v");

  vector<Real> sig0_fact = {1, 1.05};
  vector<String> sig0_name = {"s0", "s1"};
  vector<Complex> k0_vec = {
      0.17035918970056643 - i_ * 0.06871976789571771,
      0.501170235519639 - i_ * 0.02987688354121789,
      0.827284936320656 - i_ * 0.009463252032254791,
      1.1460876926752444 - i_ * 0.0024449822275129105,
      1.4583612406212996 - i_ * 0.0005494081963914154,
      1.7666093821137823 - i_ * 0.00011330489747763084,
      2.072620745403709 - i_ * 2.2173316592827555e-05,
      2.377350096298269 - i_ * 4.1937434161936945e-06,
      2.681287566664264 - i_ * 7.744879421209387e-07,
      2.984699668096411 - i_ * 1.405234448913448e-07,
      3.2877432588219477 - i_ * 2.51503198132076e-08,
      3.5905173384492284 - i_ * 4.452436062958641e-09,
      3.8930878895961456 - i_ * 7.812273373381207e-10,
      4.1955007721083115 - i_ * 1.3605740018399644e-10,
      4.497788931598282 - i_ * 2.3548818960193295e-11,
      4.799976681986642 - i_ * 4.051834862944728e-12};

  for (int is0 = 0; is0 < 2; is0++)
  {
    std::ofstream myfile;
    myfile.precision(17);
    String name("data/eigs_" + sig0_name[is0]);
    myfile.open(name.c_str());

    for (int ik0 = 0; ik0 < k0_vec.size(); ik0++)
    {
      Complex k0 = k0_vec[ik0];
      Real s0 = min(sig0_fact[is0] * 10 / (r1 * k0.real()), 12.5);
      Complex sigma = k0 * k0 + 0.1 * i_;

      std::cout << "m = " << ik0 << std::endl
                << "k0 = " << k0 << std::endl
                << "s0 = " << s0 << std::endl
                << "sigma = " << sigma << std::endl
                << std::endl;

      pars << Parameter(s0, "s0");

      Function Fstif(MatPML, pars);
      Function Fmass(FctPML, pars);

      //// Bilinear forme
      BilinearForm auv = (1 / eps) * intg(Cavity, grad(u) | grad(v)) + intg(Vacuum, grad(u) | grad(v)) + intg(PML, Fstif * grad(u) | grad(v));
      BilinearForm buv = intg(Cavity, u * v) + intg(Vacuum, u * v) + intg(PML, Fmass * u * v);
      EssentialConditions ecs = (u | SigPML = 0);

      //// TermMatrix and Solver
      cpuTime();
      TermMatrix A(auv, ecs, ReductionMethod(_pseudoReduction, 100.), "A");
      TermMatrix B(buv, ecs, ReductionMethod(_pseudoReduction, 0.01), "B");
      cpuTime("Time Terms");

      EigenElements eigs = eigenSolve(A, B, _nev = nbvp, _sigma = sigma);
      cpuTime("Time Arpack");

      Complex k;
      for (int j = 1; j <= nbvp; j++)
      {
        k = std::sqrt(eigs.value(j));
        myfile << ik0 << " " << j << " " << k.real() << " " << k.imag() << endl;
        saveToFile("data/eig_" + sig0_name[is0] + "_m" + std::to_string(ik0 + 1) + "_j" + std::to_string(j), eigs.vector(j), _vtu);
      }
    }
    myfile.close();
  }

  return 0;
}
