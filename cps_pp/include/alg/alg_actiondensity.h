#ifndef __ALG__ACTIONDENSITY__
#define __ALG__ACTIONDENSITY__
#include <config.h>

#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pmalloc.h> 
#include <util/gauge_field.h>
#include "alg_base.h" 
#include "common_arg.h"
#include "no_arg.h"

CPS_START_NAMESPACE

class AlgActionDensity : public Alg
{
private:

  const char *cname;

  void ZeroReal(Matrix& m);
  Complex ActionDensity(Matrix clovers[]);
  void CloverLeaf(GaugeField& gf, Matrix& pl,  int* pos, int mu, int nu);
  
public:

  AlgActionDensity(Lattice&    latt, 
             CommonArg *c_arg ):
    Alg(latt,c_arg),
    cname("AlgActionDensity")
  {;}
  
  virtual ~AlgActionDensity() {;}
  
  void run();
  void smartrun();

};

CPS_END_NAMESPACE
#endif 
