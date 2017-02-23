#include <config.h>
#include <alg/plaq.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/time_cps.h>
#include <util/error.h>
#include <comms/glb.h>
#include <util/qcdio.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
//! Performs the computation.
/*!
  The real trace of the plaquette is averaged over the lattice and normalised
  to unity.
  \post The following data are written to the file specified in the
  common_arg structure:
  -# mean plaquette
  -# variance of mean plaquette
  -# maximum plaquette
  -# minimum plaquette
  -# mean temporal plaquette
  -# mean spatial plaquette
  */
//------------------------------------------------------------------
void Plaq::
run(GaugeField& gf, const char * fn) {
  const char * fname = "run()";

  Float dtime = -dclock();

  gf.require(1);
  gf.refresh();

  int total_sites = gf.localVolume() * GJP.Xnodes() * GJP.Ynodes() * GJP.Znodes() * GJP.Tnodes();
  Float norm_fac = 1.0 / (6.0 * total_sites) ;

  Float p_temporal, p_spatial;
  p_temporal = p_spatial = 0.0;

  Float p_sum    =  0.0 ;
  Float p_sq_sum =  0.0 ;
  Float p_min    =  3.0 ;
  Float p_max    = -3.0 ;

  int x[4];
  int dirs[4];
  Matrix plaq;
  for (int index = 0; index < gf.localVolume(); index++) {
    gf.coordinatesFromIndex(x, index);
    for (int mu = 0; mu < 3; mu++) {
      for (int nu = mu+1; nu < 4; nu++) {
        dirs[0] = mu;
        dirs[1] = nu;
        dirs[2] = mu + 4;
        dirs[3] = nu + 4;
        gfPathOrdProd(plaq, gf, x, dirs, 4);
        Float tmp_flt = plaq.ReTr();
        p_sum    += tmp_flt ;
        p_sq_sum += tmp_flt * tmp_flt ;
        p_min    =  p_min<tmp_flt ? p_min : tmp_flt;
        p_max    =  p_max>tmp_flt ? p_max : tmp_flt;
        if(nu==DIR_T||mu==DIR_T) p_temporal += tmp_flt;
        else p_spatial += tmp_flt;
      }
    }
  }

  glb_sum(&p_sum) ;
  glb_sum(&p_sq_sum) ;
  glb_min(&p_min) ;
  glb_max(&p_max) ;
  glb_sum(&p_spatial);
  glb_sum(&p_temporal);

  Float one_third = 1.0 / 3.0 ;

  p_sum    *= one_third ;
  p_sq_sum *= one_third * one_third ;
  p_min    *= one_third ;
  p_max    *= one_third ;

  Float p_var = norm_fac * (1.0 + norm_fac) \
                * (p_sq_sum - norm_fac*p_sum*p_sum) ;
  p_sum *= norm_fac ;
  p_spatial *= 2.0*norm_fac*one_third;
  p_temporal *= 2.0*norm_fac*one_third;

  // Print out results
  //----------------------------------------------------------------
  FILE *fp;

  if( (fp = Fopen(fn, "a")) == NULL )
    ERR.FileA(cname,fname, fn);

  Fprintf(fp, "%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n",
      p_sum, p_var,
      p_max, p_min,
      p_temporal, p_spatial);
  Fclose(fp);

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);
}

Float Plaq::
getPlaqSum(GaugeField& gf) {
  const char * fname = "getPlaqSum()";

  Float dtime = -dclock();

  gf.require(1);
  gf.refresh();

  Float p_sum = 0.0 ;

  for (int index = 0; index < gf.localVolume(); index++) {
    int x[4];
    int dirs[4];
    Matrix plaq;
    gf.coordinatesFromIndex(x, index);
    for (int mu = 0; mu < 3; mu++) {
      for (int nu = mu+1; nu < 4; nu++) {
        dirs[0] = mu;
        dirs[1] = nu;
        dirs[2] = mu + 4;
        dirs[3] = nu + 4;
        gfPathOrdProd(plaq, gf, x, dirs, 4);
        Float tmp_flt = 1.0 - plaq.ReTr() / 3.0;
        p_sum += tmp_flt ;
      }
    }
  }

  glb_sum(&p_sum) ;

  dtime += dclock();
  print_flops(cname, fname, 0, dtime);

  return p_sum;
}


CPS_END_NAMESPACE
