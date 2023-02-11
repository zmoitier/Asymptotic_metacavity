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

// hankelH1^prime(m,z) / hankelH1(m,z)
Complex DtN_m(Int m, Real k, Real r)
{
  Real z = k * r;
  return k * (hankelH1(z, m - 1) - hankelH1(z, m + 1)) / (2 * hankelH1(z, m));
}

// Spectral basis for DtN
Real csmt(const Point &P, Parameters &pa = defaultParameters)
{
  Real T = pa("T");
  Real c0 = 1 / std::sqrt(2 * pi_ * T);
  Real cm = 1 / std::sqrt(pi_ * T);
  Real theta = std::atan2(P(2), P(1));
  int i = pa("basis index") - 1;
  int m = (i + 1) / 2;
  if (i == 0)
  {
    return c0;
  }
  else if (i % 2)
  {
    return cm * cos(m * theta);
  }
  else
  {
    return cm * sin(m * theta);
  }
}

Vector<Complex> coef_DtN(Int M, Real k, Real T)
{
  Vector<Complex> lambda(2 * M + 1);
  lambda[0] = DtN_m(0, k, T);
  for (int m = 1; m <= M; m++)
  {
    lambda[2 * m - 1] = DtN_m(m, k, T);
    lambda[2 * m] = lambda[2 * m - 1];
  }
  return lambda;
}

//// Main
int main(int argc, char **argv)
{
  init(argc, argv, _lang = en); // mandatory initialization of xlifepp

  // Parameters
  Real eps = atof(argv[1]);
  Real k = atof(argv[2]);
  Int dGL = atoi(argv[3]);

  Int M = 24; // ceil(1.5 * k * T);
  cout << "+--------------------" << endl;
  cout << "| M = " << M << " (N = " << 2 * M + 1 << ")" << endl;
  cout << "+--------------------" << endl;

  //// Mesh and domains definition
  Mesh mesh = Mesh("../../GMSH/disk_T_DtN.msh", "disk_T_DtN", msh);
  Domain Omega = mesh.domain("Omega");
  Domain Cavity = mesh.domain("Cavity");
  Domain Vacuum = mesh.domain("Vacuum");
  Domain Gamma = mesh.domain("Gamma");
  Domain Sigma = mesh.domain("Sigma");

  Real T = get_radius(Sigma);

  Parameters params;
  params << Parameter(M, "M") << Parameter(k, "k") << Parameter(T, "T");

  //// Space and unknown
  Space V(Omega, interpolation(Lagrange, _GaussLobattoPoints, dGL, H1), "V");
  Unknown u(V, "u");
  TestFunction v(u, "v");

  //// Create spectral space and DtN kernel
  Space Sp(Sigma, Function(csmt, params), 2 * M + 1, "cos/sin(m*theta)");
  Unknown phi(Sp, "phi");
  TensorKernel tk(phi, coef_DtN(M, k, T));

  cout << "k = " << k << endl;

  //// Create bilinear form and linear form
  BilinearForm auv = (1 / eps) * intg(Cavity, grad(u) | grad(v)) +
                     intg(Vacuum, grad(u) | grad(v)) -
                     (k * k) * intg(Omega, u * v) -
                     intg(Sigma, Sigma, u * tk * v);

  //// TermMatrix and Solver
  cpuTime();
  TermMatrix A(auv, "A");
  cpuTime("Time Terms");
  A.saveToFile("tmp_matrix", _coo);
  cpuTime("Time Write");

  return 0;
}
