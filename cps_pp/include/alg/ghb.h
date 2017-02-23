#ifndef INCLUDED_GHB_H
#define INCLUDED_GHB_H

#include <config.h>
#include <util/gauge_field.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/ghb_arg.h>
#include <vector>
CPS_START_NAMESPACE

//-----------------------------------------------------------------------------
//! Class implementing the gauge field global heatbath algorithm.
/*!
  The algorithm used is the Cabbibo-Marinari SU(N) heatbath update with
  the Kennedy-Pendleton method for updating SU(2) subgroups

  \ingroup alg 
 */
//------------------------------------------------------------------
class GheatBath {

 private:

    const char *cname;

    GhbArg *alg_ghb_arg;
        // The argument structure for the GheatBath algorithm

    void UpdateLink(Matrix* mp, const Matrix& stap, double beta, int rgen_pos_4d);

    bool offsets_list_initialized;
    std::vector<Offsets> offsets_list;

 public:

    GheatBath(GhbArg& arg) {
      cname = "GheatBath";
      alg_ghb_arg = &arg;
      offsets_list_initialized = false;
    }

    ~GheatBath() {
    }

    void run(GaugeField& gf, GaugeAction& ga);

};

CPS_END_NAMESPACE
#endif
