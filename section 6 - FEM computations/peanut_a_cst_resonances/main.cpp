#include "xlife++.h"
#include <fstream>
#include <iostream>
#include <tuple>
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

auto read_data(string filename)
{
  std::vector<Int> m_vec, j_vec;
  std::vector<Complex> k_vec;
  string line;
  ifstream myfile(filename);
  if (myfile.is_open())
  {
    while (getline(myfile, line))
    {
      std::istringstream iss(line);
      Int m, j;
      Real kr, ki;
      if (!(iss >> m >> j >> kr >> ki))
      {
        cout << "Can not parse the line: " << line << '\n';
        break;
      } // error
      m_vec.push_back(m);
      j_vec.push_back(j);
      k_vec.push_back(kr + i_ * ki);
    }
    myfile.close();
  }
  else
  {
    cout << "Unable to open file" << '\n';
  }
  return std::make_tuple(m_vec, j_vec, k_vec);
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
  String filename = "peanut_eCst_-1.1_pla";
  Int nev = 4; // > 2

  //// Mesh and domains definition
  Mesh mesh = Mesh("../../GMSH/peanut_T_PML_3o10.msh", "peanut_T_PML_3o10", msh);
  Domain Omega = mesh.domain("Omega");
  Domain Cavity = mesh.domain("Cavity");
  Domain Vacuum = mesh.domain("Vacuum");
  Domain PML = mesh.domain("PML");
  Domain Gamma = mesh.domain("Gamma");
  Domain Sigma = mesh.domain("Sigma");
  Domain SigPML = mesh.domain("SigPML");
  Domain OmegaV = merge(Cavity, Vacuum, "OmegaV");

  //// PML Parameters
  Real r0 = get_radius(Sigma);
  Real r1 = get_radius(SigPML);

  Parameters pars;
  pars << Parameter(r0, "r0") << Parameter(r1, "r1");

  //// Space and unknown
  Space V(Omega, interpolation(Lagrange, _GaussLobattoPoints, 8, H1), "V");
  Unknown u(V, "u");
  TestFunction v(u, "v");

  Int dofs = V.nbDofs();
  Int ncv = max((Int)20, min(dofs, 2 * nev + 1));
  Int max_iter = 10 * dofs;

  vector<Real> sig0_fact = {1, 1.05};
  vector<String> sig0_name = {"s0", "s1"};
  std::vector<Int> m_vec, j_vec;
  std::vector<Complex> k0_vec;
  tie(m_vec, j_vec, k0_vec) = read_data(filename);

  for (int is0 = 0; is0 < 2; is0++)
  {
    std::ofstream myfile;
    myfile.precision(17);
    String name("data/eigs_" + sig0_name[is0]);
    myfile.open(name.c_str());

    for (int ik0 = 0; ik0 < k0_vec.size(); ik0++)
    {
      Complex k0 = k0_vec[ik0];
      Real s0 = sig0_fact[is0] * (10 + r1 * std::abs(k0.imag())) / (r1 * k0.real());
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

      EigenElements eigs = eigenSolve(A, B, _nev = nev, _ncv = ncv, _maxIt = max_iter, _sigma = sigma, _sort = _incr_realpart);
      cpuTime("Time Arpack");

      Int i1;
      Complex k1;
      Real min1 = 1000;
      for (int j = 1; j <= nev; j++)
      {
        Complex k = std::sqrt(eigs.value(j));
        Real m = std::abs(k - k0);
        if (m < min1)
        {
          k1 = k;
          i1 = j;
          min1 = m;
        }
      }
      myfile << m_vec[ik0] << " " << j_vec[ik0] << " " << k1.real() << " " << k1.imag() << '\n';
      saveToFile("data/eig_" + sig0_name[is0] + "_m" + std::to_string(m_vec[ik0]) + "_j" + std::to_string(j_vec[ik0]), eigs.vector(i1).onDomain(OmegaV), _vtu);
    }
    myfile.close();
  }

  return 0;
}
