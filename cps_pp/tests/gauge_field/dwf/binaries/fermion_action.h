#ifndef INCLUDED_FERMION_ACTION_H
#define INCLUDED_FERMION_ACTION_H

#include <util/gauge_field.h>
#include <alg/fermion_vector.h>

CPS_START_NAMESPACE

class FermionActionDomainWall;

class FermionAction {

  public :

    FermionAction(FclassType _fclass) {
      cname = "FermionAction";
      fclass = _fclass;
    }

    FermionActionDomainWall& getFermionActionDomainWall() {
      return *(FermionActionDomainWall *)this;
    }

    FclassType Fclass() {
      return fclass;
    }

  private :

    const char * cname;

    FclassType fclass;
};

class FermionActionDomainWall : public FermionAction {

  public :

    FermionActionDomainWall(int _s) : FermionAction(F_CLASS_DWF) {
      s = _s;
    }

    int getS() {
      return s;
    }

  private :

    const char * cname;

    int s;
};

void conjugateGradient(GaugeField& gf, FermionAction& fa,
    FermionVectorTp& source, FermionVectorTp& sol, CgArg& cg_arg,
		int& iter, Float& true_res);

void conjugateGradientDomainWall(GaugeField& gf, FermionActionDomainWall& fa,
    FermionVectorTp& source, FermionVectorTp& sol, FermionVectorTp& midsol,
    CgArg& cg_arg, int& iter, Float& true_res);

CPS_END_NAMESPACE
#endif

