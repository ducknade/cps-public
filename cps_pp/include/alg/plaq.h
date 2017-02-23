#ifndef INCLUDED_PLAQ_H
#define INCLUDED_PLAQ_H

#include <config.h>
#include <util/gauge_field.h>
#include <util/smalloc.h>
#include <util/pmalloc.h>
#include <alg/alg_base.h>
#include <alg/common_arg.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//! A class implementing calculation of the average plaquette.
/*!
  This class computes the real trace of the plaquette averaged over
  the total number of plaquettes and the number of colours (which is three).
  Also computed is the variance of this mean, and one third of the real trace
  of the plaquette at the origin in the X-Y plane.

  \ingroup alg
*/
//------------------------------------------------------------------
class Plaq {

 private:

    const char * cname;

 public:

    Plaq() {
      cname = "Plaq";
    }

    ~Plaq() {
    }

    // fn is the output filename
    void run(GaugeField& gf, const char * fn);

    Float getPlaqSum(GaugeField& gf);
    // return sum of (1.0 - plaq.ReTr() / 3.0)

};

CPS_END_NAMESPACE

#endif
