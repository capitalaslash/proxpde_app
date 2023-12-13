#include <proxpde/def.hpp>

#include <proxpde/assembly.hpp>
#include <proxpde/bc.hpp>
#include <proxpde/builder.hpp>
#include <proxpde/fe.hpp>
#include <proxpde/fespace.hpp>
#include <proxpde/iomanager.hpp>
#include <proxpde/mesh.hpp>

int main()
{
  using namespace proxpde;

  using Elem_T = Quad;
  using Mesh_T = Mesh<Elem_T>;
  using FESpace_T = FESpace<
      Mesh_T,
      LagrangeFE<Elem_T, 1>::RefFE_T,
      LagrangeFE<Elem_T, 1>::RecommendedQR>;

  auto const rhs = [](Vec3 const & p) { return M_PI * std::sin(M_PI * p(0)); };

  auto const exactSol = [](Vec3 const & p)
  { return std::sin(M_PI * p(0)) / M_PI + p(0); };

  uint const numElemsX = 20U;
  uint const numElemsY = 20U;

  Vec3 const origin{0.0, 0.0, 0.0};
  Vec3 const length{1.0, 1.0, 0.0};

  std::unique_ptr<Mesh_T> mesh{new Mesh_T};

  buildHyperCube(*mesh, origin, length, {numElemsX, numElemsY, 0U});

  FESpace_T feSpace{*mesh};

  auto bc = BCEss{feSpace, side::LEFT};
  bc << [](Vec3 const &) { return 0.0; };

  AssemblyStiffness stiffness{1.0, feSpace};
  AssemblyAnalyticRhs f{rhs, feSpace};

  Builder builder{feSpace.dof.size};
  builder.buildLhs(std::tuple{stiffness}, std::vector{bc});
  builder.buildRhs(std::tuple{f}, std::vector{bc});
  builder.closeMatrix();

  Var sol{"u"};
  LUSolver solver;
  solver.analyzePattern(builder.A);
  solver.factorize(builder.A);
  sol.data = solver.solve(builder.b);

  Var exact{"exact"};
  interpolateAnalyticFunction(exactSol, feSpace, exact.data);

  Var error{"error"};
  error.data = sol.data - exact.data;
  fmt::print("error: {}\n", error.data.norm());

  IOManager io{feSpace, "output/sol"};
  io.print({sol, exact, error});

  return 0;
}
