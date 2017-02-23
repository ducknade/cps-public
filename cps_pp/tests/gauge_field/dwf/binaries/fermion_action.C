#include <util/gauge_field.h>
#include "fermion_action.h"

#include <alg/fermion_vector.h>

CPS_START_NAMESPACE

void conjugateGradient(GaugeField& gf, FermionAction& fa,
    FermionVectorTp& source, FermionVectorTp& sol, CgArg& cg_arg,
		int& iter, Float& true_res) {

  const char *fname = "conjugateGradient(gf&, fa&, source&, sol&, int&, Float&)";
  VRB.Func("Global", fname);

  // Set the node size of the full (non-checkerboarded) fermion field
  //----------------------------------------------------------------
  // Do inversion
  //----------------------------------------------------------------
  switch (fa.Fclass()) {
    case F_CLASS_DWF :
      {
        FermionVectorTp midsol;
        conjugateGradientDomainWall(gf, fa.getFermionActionDomainWall(),
            source, sol, midsol, cg_arg, iter, true_res);
      }
      break;
    default :
      ERR.NotImplemented("", fname, "Fermion class\n");
      break;
  }
}

void conjugateGradientDomainWall(GaugeField& gf, FermionActionDomainWall& fa,
    FermionVectorTp& source, FermionVectorTp& sol, FermionVectorTp& midsol,
    CgArg& cg_arg, int& iter, Float& true_res) {

  const char *fname = "conjugateGradientDomainWall(gf&, fa&, source&, sol&, midsol&, int&, Float&)";
  VRB.Func("Global", fname);

  GnoneFdwf lat;
  gf.setLattice(lat);

  int ls = fa.getS();
  int ls_glb = ls;
  int f_size = 4 * sizeof(Vector);
  int f_size_5d = f_size * ls;

  Vector *src_4d    = (Vector*)source.data();
  Vector *sol_4d    = (Vector*)sol.data();
  Vector *midsol_4d = (Vector*)midsol.data();
  Vector *src_5d    = (Vector*)smalloc("Global", fname, "src_5d", f_size_5d * sizeof(IFloat));
  Vector *sol_5d    = (Vector*)smalloc("Global", fname, "sol_5d", f_size_5d * sizeof(IFloat));

  lat.Ffour2five(src_5d, src_4d, 0, ls_glb-1);
  lat.Ffour2five(sol_5d, sol_4d, ls_glb-1, 0);

  iter = lat.FmatInv(sol_5d, src_5d, &cg_arg, &true_res,
      CNV_FRM_YES, PRESERVE_NO);

  // prop on walls
  lat.Ffive2four(sol_4d, sol_5d, ls_glb-1, 0);
  lat.Ffive2four(midsol_4d, sol_5d, ls_glb/2-1, ls_glb/2);

  sfree("Global",fname, "sol_5d", sol_5d);
  sfree("Global",fname, "src_5d", src_5d);
}

CPS_END_NAMESPACE
