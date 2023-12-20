#include <proxpde/def.hpp>

#include <proxpde/assembly.hpp>
#include <proxpde/bc.hpp>
#include <proxpde/builder.hpp>
#include <proxpde/fe.hpp>
#include <proxpde/fespace.hpp>
#include <proxpde/iomanager.hpp>
#include <proxpde/mesh.hpp>

using namespace proxpde;

using Elem_T = Quad;
using Mesh_T = Mesh<Elem_T>;
using FESpace_T = FESpace<
    Mesh_T,
    LagrangeFE<Elem_T, 1>::RefFE_T,
    LagrangeFE<Elem_T, 1>::RecommendedQR>;

using FESpaceRT0_T = FESpace<
    Mesh_T,
    RaviartThomasFE<Elem_T, 0>::RefFE_T,
    RaviartThomasFE<Elem_T, 0>::RecommendedQR>;

struct MySolver
{
  MySolver() = default;

  void init();

  void test();

  std::unique_ptr<Mesh_T> mesh;
  FESpace_T feSpace;
  FEVar<FESpace_T> u;

  FESpaceRT0_T feSpaceRT0;
  FEVar<FESpaceRT0_T> uRT0;
  FEVar<FESpaceRT0_T> uRT0Old;
};

void MySolver::init()
{
  uint const numElemsX = 20U;
  uint const numElemsY = 20U;

  Vec3 const origin{-0.5, -0.5, 0.0};
  Vec3 const length{1.0, 1.0, 0.0};

  mesh.reset(new Mesh_T);

  buildHyperCube(
      *mesh,
      origin,
      length,
      {numElemsX, numElemsY, 0U},
      MeshFlags::INTERNAL_FACETS | MeshFlags::FACET_PTRS | MeshFlags::NORMALS);

  feSpace.init(*mesh);
  u.init("u", feSpace);

  feSpaceRT0.init(*mesh);
  uRT0.init("uRT0", feSpaceRT0);
  uRT0Old.init("uRT0Old", feSpaceRT0);

  auto const fun = [](Vec3 const & p)
  {
    double const x = p[0] + 0.5;
    double const y = p[1] + 0.5;
    return Vec3{
        -2.0 * std::sin(M_PI * x) * std::sin(M_PI * x) * std::sin(M_PI * y) *
            std::cos(M_PI * y),
        2.0 * std::sin(M_PI * x) * std::cos(M_PI * x) * std::sin(M_PI * y) *
            std::sin(M_PI * y),
        0.0};
  };

  interpolateAnalyticFunction(fun, *uRT0.feSpace, uRT0.data);
  uRT0Old.data = uRT0.data;
}

void MySolver::test()
{
  id_T const elemId = 200U;
  Elem_T const & elem = mesh->elementList[elemId];

  auto const ptHat = Vec2{0.0, 0.0};

  uRT0.reinit(elem);
  Vec3 uRT0Pt = uRT0.evaluateOnRef(ptHat);
  fmt::print("uRT0Pt: {}\n", uRT0Pt.transpose());

  uRT0Old.reinit(elem);
  Vec3 uRT0OldPt = uRT0Old.evaluateOnRef(ptHat);
  fmt::print("uRT0OldPt: {}\n", uRT0OldPt.transpose());
}

int main()
{
  MySolver s;

  s.init();

  s.test();

  return 0;
}
