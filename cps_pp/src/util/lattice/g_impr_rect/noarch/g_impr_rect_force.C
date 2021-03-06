#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Implementation of GimprRect class.

  $Id: g_impr_rect_force.C,v 1.3 2008/09/18 15:23:17 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/09/18 15:23:17 $
//  $Header: /space/cvs/cps/cps++/src/util/lattice/g_impr_rect/noarch/g_impr_rect_force.C,v 1.3 2008/09/18 15:23:17 chulwoo Exp $
//  $Id: g_impr_rect_force.C,v 1.3 2008/09/18 15:23:17 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.3 $
//  $Source: /space/cvs/cps/cps++/src/util/lattice/g_impr_rect/noarch/g_impr_rect_force.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <stdlib.h>
#include <math.h>
#include <util/lattice.h>
#include <util/verbose.h>
#include <util/vector.h>
#include <util/gjp.h>
#include <util/gw_hb.h>
#include <util/time_cps.h>
#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
#include <util/gauge_field.h>
CPS_START_NAMESPACE

ForceArg GimprRect::EvolveMomGforce(Matrix *mom, Float dt) {
  const char * fname = "EvolveMomGforce(Matrix *mom, Float dt)";
  static bool initilized = false;
  static Offsets offsets;
  cps::GaugeField gf(*this);
  GaugeAction& ga = *(getGaugeAction(*this));
  gf.resize(gaInteractionRange(ga));
  if (!initilized) {
    gf.refresh();
    gf.recordOffsetStart(false, true);
  } else {
    gf.refresh(offsets);
  }
  ForceArg farg = gfEvolveMomGforce(mom, gf, ga, dt);
  if (!initilized) {
    offsets = gf.getRecordedOffsets();
    gf.recordOffsetEnd();
    initilized = true;
  }
  delete &ga;
  return farg;
}

// //------------------------------------------------------------------------------
// enum { MATRIX_SIZE = 18 };
//
//
// #define PROFILE
// //------------------------------------------------------------------------------
// // EvolveMomGforce(Matrix *mom, Float dt):
// // It evolves the canonical momentum mom by dt
// // using the pure gauge force.
// //------------------------------------------------------------------------------
//
//
// ForceArg GimprRect::EvolveMomGforce(Matrix *mom, Float dt){
//   char *fname = "EvolveMomGforce(M*,F)";
//   VRB.Func(cname,fname);
//
// #ifdef PROFILE
//   Float time = -dclock();
//   ForceFlops = 0;
// #endif
//
//   Float L1=0.0;
//   Float L2=0.0;
//   Float Linf=0.0;
//
// // Not doing it in parallel because of the
// // lat.Staple and lat.RectStaple are not thread safe
// //#pragma omp parallel
//   {
//     Float pL1=0.0;
//     Float pL2=0.0;
//     Float pLinf=0.0;
// //#pragma omp for nowait
//     for (int offset = 0; offset < 4 * GJP.VolNodeSites(); offset++) {
//       //  for(x[0] = 0; x[0] < GJP.XnodeSites(); ++x[0])
//       //  for(x[1] = 0; x[1] < GJP.YnodeSites(); ++x[1])
//       //  for(x[2] = 0; x[2] < GJP.ZnodeSites(); ++x[2])
//       //  for(x[3] = 0; x[3] < GJP.TnodeSites(); ++x[3]) {
//
//       //    for (int mu = 0; mu < 4; ++mu) {
//       int mu = offset % 4;
//       int x[4];
//       CoordsFromOffset(x, offset / 4);
//       //int uoff = GsiteOffset(x);
//
//       Matrix mt0;
//       Matrix mt1;
//       Matrix mt2;
//       Matrix *mp0 = &mt0;		// ihdot
//       Matrix *mp1 = &mt1;
//       Matrix *mp2 = &mt2;
//
// //#pragma omp critical
//       GforceSite(*mp0, x, mu);
//
//       //IFloat *ihp = (IFloat *)(mom+uoff+mu);
//       IFloat *ihp = (IFloat *)(mom+offset);
//       IFloat *dotp = (IFloat *)mp0;
//       fTimesV1PlusV2(ihp, dt, dotp, ihp, 18);
//       Float norm = ((Matrix*)dotp)->norm();
//       Float tmp = sqrt(norm);
//       pL1 += tmp;
//       pL2 += norm;
//       pLinf = (tmp>pLinf ? tmp : pLinf);
//       //}}
//     }
// //#pragma omp critical
//     {
//       L1 += pL1;
//       L2 += pL2;
//       Linf = pLinf > Linf ? pLinf : Linf;
//     }
//   }
//
//   ForceFlops +=GJP.VolNodeSites()*4*18*2;
// #ifdef PROFILE
//   time += dclock();
//   print_flops(cname,fname,ForceFlops,time);
// #endif
//
//   glb_sum(&L1);
//   glb_sum(&L2);
//   glb_max(&Linf);
//
//   L1 /= 4.0*GJP.VolSites();
//   L2 /= 4.0*GJP.VolSites();
//
//   VRB.FuncEnd(cname,fname);
//   return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);
// }

CPS_END_NAMESPACE
