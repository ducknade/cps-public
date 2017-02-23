// vim: set ts=2 sw=2 expandtab:
#include<config.h>
/*!\file
  \brief Definitions of the GheatBath class methods.
  
  $Id: alg_ghb.C,v 1.14 2006/07/03 04:48:06 chulwoo Exp $
*/

#include <stdlib.h>	// exit()
#include <util/qcdio.h>
#include <util/time_cps.h>
#include <math.h>
#include <time.h>
#include <alg/alg_ghb.h>
#include <alg/ghb.h>
#include <alg/common_arg.h>
#include <alg/ghb_arg.h>
#include <util/lattice.h>
#include <util/gauge_field.h>
#include <util/gjp.h>
#include <util/random.h>
#include <util/smalloc.h>
#include <util/vector.h>
#include <util/verbose.h>
#include <util/error.h>
#include <mem/p2v.h>
#include <vector>
#include <omp.h>
#if TARGET != NOARCH
#include <qmp.h>
#endif
CPS_START_NAMESPACE


void GheatBath::
run(GaugeField& gf, GaugeAction& ga) {
  using namespace std;

  char *fname = "run()";
  VRB.Func(cname,fname);

  LRG.SetInterval(1,-1);
  gf.require(gaInteractionRange(ga));

  // Run the gauge heat bath
  //----------------------------------------------------------------
  for(int i=0; i< alg_ghb_arg->num_iter; i++) {
    // Heat bath
    // Checkerboard everything, do even or odd sites
    // Scan over local subvolume doing all even(odd) links
    // index for "x" is consistant with GsiteOffset: "x,y,z,t" order

    // The y[4] are the local coordinates for the 2^4 cube. 
    // We traverse the points on this cube, and do the corresponding
    // points on every other hypercube before moving on.

    if (ga.Gclass() == G_CLASS_WILSON) {
      offsets_list.resize(32);
      Float dtime = 0.0;
      for (int checker = 0; checker < 2; checker++) {
        // checker board on 2^4 cube
        // assuming every dimension can be divided by 4
        for (int lindex = 0; lindex < 16; lindex++) {

          if (offsets_list_initialized) {
            gf.refresh(offsets_list[checker*16+lindex]);
          } else {
            gf.refresh();
            gf.recordOffsetStart(false, true);
          }

          dtime -= dclock();
#pragma omp parallel for
          for (int cindex = 0; cindex < gf.localVolume() / 16; cindex++) {

            //          if (QMP_get_node_number() == 0) {
            //#pragma omp critical
            //            {
            //              printf("calc %2d : %1d %2d %2d\n", omp_get_thread_num(), checker, lindex, cindex);
            //            }
            //          }

            int crem = cindex;
            int cx[4]; // coordinates of the 2^4 cube
            cx[0] = crem % (gf.nodeSites(0) / 2);
            crem /= gf.nodeSites(0) / 2;
            cx[1] = crem % (gf.nodeSites(1) / 2);
            crem /= gf.nodeSites(1) / 2;
            cx[2] = crem % (gf.nodeSites(2) / 2);
            crem /= gf.nodeSites(2) / 2;
            cx[3] = crem % (gf.nodeSites(3) / 2);
            if ((cx[0] + cx[1] + cx[2] + cx[3]) % 2 != checker) {
              continue;
            }
            int lrem = lindex;
            int x[4]; // final site coordinates
            x[0] = cx[0] * 2 + lrem % 2;
            lrem /= 2;
            x[1] = cx[1] * 2 + lrem % 2;
            lrem /= 2;
            x[2] = cx[2] * 2 + lrem % 2;
            lrem /= 2;
            x[3] = cx[3] * 2 + lrem % 2;

            // VRB.Result(cname, fname, "ghb %d %d %d %d\n", x[0], x[1], x[2], x[3]);

            int rgen_pos_4d = LRG.getGeneratorPos4d(x);
            for( int mu=0; mu<4; mu++ ) {
              Matrix mStaple;

              gfAllStaple(mStaple, gf, ga, x, mu);

              Matrix * pmLink = gf.unsafeGetLink(x, mu);

              UpdateLink(pmLink, mStaple, ga.getBeta(), rgen_pos_4d);
            }
          }
          if (!offsets_list_initialized) {
            offsets_list[checker*16+lindex] = gf.getRecordedOffsets();
            gf.recordOffsetEnd();
          }
          dtime += dclock();
        }
      }
      offsets_list_initialized = true;
      print_flops(cname, fname, 0, dtime);
    } else if (ga.Gclass() == G_CLASS_IMPR_RECT) {
      offsets_list.resize(32);
      Float dtime = 0.0;
      for (int checker = 0; checker < 2; checker++) {
        // checker board on 2^4 cube
        // assuming every dimension can be divided by 4
        for (int lindex = 0; lindex < 16; lindex++) {

          if (offsets_list_initialized) {
            gf.refresh(offsets_list[checker*16+lindex]);
          } else {
            gf.refresh();
            gf.recordOffsetStart(false, true);
          }

          dtime -= dclock();
#pragma omp parallel for
          for (int cindex = 0; cindex < gf.localVolume() / 16; cindex++) {

            //          if (QMP_get_node_number() == 0) {
            //#pragma omp critical
            //            {
            //              printf("calc %2d : %1d %2d %2d\n", omp_get_thread_num(), checker, lindex, cindex);
            //            }
            //          }

            int crem = cindex;
            int cx[4]; // coordinates of the 2^4 cube
            cx[0] = crem % (gf.nodeSites(0) / 2);
            crem /= gf.nodeSites(0) / 2;
            cx[1] = crem % (gf.nodeSites(1) / 2);
            crem /= gf.nodeSites(1) / 2;
            cx[2] = crem % (gf.nodeSites(2) / 2);
            crem /= gf.nodeSites(2) / 2;
            cx[3] = crem % (gf.nodeSites(3) / 2);
            if ((cx[0] + cx[1] + cx[2] + cx[3]) % 2 != checker) {
              continue;
            }
            int lrem = lindex;
            int x[4]; // final site coordinates
            x[0] = cx[0] * 2 + lrem % 2;
            lrem /= 2;
            x[1] = cx[1] * 2 + lrem % 2;
            lrem /= 2;
            x[2] = cx[2] * 2 + lrem % 2;
            lrem /= 2;
            x[3] = cx[3] * 2 + lrem % 2;

            // VRB.Result(cname, fname, "ghb %d %d %d %d\n", x[0], x[1], x[2], x[3]);

            int rgen_pos_4d = LRG.getGeneratorPos4d(x);
            for( int mu=0; mu<4; mu++ ) { // Local link must not be bufferred.
              Matrix mStaple;

              gfAllStaple(mStaple, gf, ga, x, mu);

              Matrix * pmLink = gf.unsafeGetLink(x, mu);

              UpdateLink(pmLink, mStaple, ga.getBeta(), rgen_pos_4d);
            }
          }
          if (!offsets_list_initialized) {
            offsets_list[checker*16+lindex] = gf.getRecordedOffsets();
            gf.recordOffsetEnd();
          }
          dtime += dclock();
        }
      }
      offsets_list_initialized = true;
      print_flops(cname, fname, 0, dtime);
    } else {
      ERR.NotImplemented(cname,fname," Do not support this type of action yet\n");
    }
  }
  gfUnitarize(gf);
}

void m_conjugate( Float* );
void cmhb_kernel( Float*, Float* );
void metropolis_kernel( Float*, Float*, double beta, int rgen_pos_4d);

void GheatBath::
UpdateLink(Matrix * pmLink, const Matrix & mStaple, double beta, int rgen_pos_4d) {

          Float* pfStapleCMHBorder;

/*  # ifdef _TARTAN
	  extern Float* fast_sigma; 
	  pfStapleCMHBorder =  fast_sigma; 
# else
*/
	  Float staple_buf[18];
	  pfStapleCMHBorder = (Float*)staple_buf;

// # endif

          *(pfStapleCMHBorder +  0) = *( (Float*)&mStaple +  0 );
          *(pfStapleCMHBorder +  9) = *( (Float*)&mStaple +  1 );
          *(pfStapleCMHBorder +  1) = *( (Float*)&mStaple +  2 );
          *(pfStapleCMHBorder + 10) = *( (Float*)&mStaple +  3 );
          *(pfStapleCMHBorder +  2) = *( (Float*)&mStaple +  4 );
          *(pfStapleCMHBorder + 11) = *( (Float*)&mStaple +  5 );
          *(pfStapleCMHBorder +  3) = *( (Float*)&mStaple +  6 );
          *(pfStapleCMHBorder + 12) = *( (Float*)&mStaple +  7 );
          *(pfStapleCMHBorder +  4) = *( (Float*)&mStaple +  8 );
          *(pfStapleCMHBorder + 13) = *( (Float*)&mStaple +  9 );
          *(pfStapleCMHBorder +  5) = *( (Float*)&mStaple + 10 );
          *(pfStapleCMHBorder + 14) = *( (Float*)&mStaple + 11 );
          *(pfStapleCMHBorder +  6) = *( (Float*)&mStaple + 12 );
          *(pfStapleCMHBorder + 15) = *( (Float*)&mStaple + 13 );
          *(pfStapleCMHBorder +  7) = *( (Float*)&mStaple + 14 );
          *(pfStapleCMHBorder + 16) = *( (Float*)&mStaple + 15 );
          *(pfStapleCMHBorder +  8) = *( (Float*)&mStaple + 16 );
          *(pfStapleCMHBorder + 17) = *( (Float*)&mStaple + 17 );

          Float* pfLinkCMHBorder;

/*  # ifdef _TARTAN
	  extern Float* fast_link; 
	  pfLinkCMHBorder =  fast_link; 
# else
*/
	  Float link_buf[18];
	  pfLinkCMHBorder = (Float*)link_buf;

//	  printf("link_buf = %p  pfLinkCMHBorder = %p pmLink = %p\n",link_buf,pfLinkCMHBorder,pmLink);
// # endif

          *(pfLinkCMHBorder +  0) = *((Float*)pmLink +  0 );
          *(pfLinkCMHBorder +  9) = *((Float*)pmLink +  1 );
          *(pfLinkCMHBorder +  1) = *((Float*)pmLink +  2 );
          *(pfLinkCMHBorder + 10) = *((Float*)pmLink +  3 );
          *(pfLinkCMHBorder +  2) = *((Float*)pmLink +  4 );
          *(pfLinkCMHBorder + 11) = *((Float*)pmLink +  5 );
          *(pfLinkCMHBorder +  3) = *((Float*)pmLink +  6 );
          *(pfLinkCMHBorder + 12) = *((Float*)pmLink +  7 );
          *(pfLinkCMHBorder +  4) = *((Float*)pmLink +  8 );
          *(pfLinkCMHBorder + 13) = *((Float*)pmLink +  9 );
          *(pfLinkCMHBorder +  5) = *((Float*)pmLink + 10 );
          *(pfLinkCMHBorder + 14) = *((Float*)pmLink + 11 );
          *(pfLinkCMHBorder +  6) = *((Float*)pmLink + 12 );
          *(pfLinkCMHBorder + 15) = *((Float*)pmLink + 13 );
          *(pfLinkCMHBorder +  7) = *((Float*)pmLink + 14 );
          *(pfLinkCMHBorder + 16) = *((Float*)pmLink + 15 );
          *(pfLinkCMHBorder +  8) = *((Float*)pmLink + 16 );
          *(pfLinkCMHBorder + 17) = *((Float*)pmLink + 17 );

          // Call the heat bath
          {
	    metropolis_kernel( pfStapleCMHBorder, pfLinkCMHBorder, beta, rgen_pos_4d);
//            cmhb_kernel( pfStapleCMHBorder, pfLinkCMHBorder );
          }

          // Copy the link back into the lattice, and
          // arrange the internal storage order back to
          // the "canonical" order.

          *((Float*)pmLink +  0 ) = *(pfLinkCMHBorder +  0);
          *((Float*)pmLink +  1 ) = *(pfLinkCMHBorder +  9);
          *((Float*)pmLink +  2 ) = *(pfLinkCMHBorder +  1);
          *((Float*)pmLink +  3 ) = *(pfLinkCMHBorder + 10);
          *((Float*)pmLink +  4 ) = *(pfLinkCMHBorder +  2);
          *((Float*)pmLink +  5 ) = *(pfLinkCMHBorder + 11);
          *((Float*)pmLink +  6 ) = *(pfLinkCMHBorder +  3);
          *((Float*)pmLink +  7 ) = *(pfLinkCMHBorder + 12);
          *((Float*)pmLink +  8 ) = *(pfLinkCMHBorder +  4);
          *((Float*)pmLink +  9 ) = *(pfLinkCMHBorder + 13);
          *((Float*)pmLink + 10 ) = *(pfLinkCMHBorder +  5);
          *((Float*)pmLink + 11 ) = *(pfLinkCMHBorder + 14);
          *((Float*)pmLink + 12 ) = *(pfLinkCMHBorder +  6);
          *((Float*)pmLink + 13 ) = *(pfLinkCMHBorder + 15);
          *((Float*)pmLink + 14 ) = *(pfLinkCMHBorder +  7);
          *((Float*)pmLink + 15 ) = *(pfLinkCMHBorder + 16);
          *((Float*)pmLink + 16 ) = *(pfLinkCMHBorder +  8);
          *((Float*)pmLink + 17 ) = *(pfLinkCMHBorder + 17);

}

CPS_END_NAMESPACE
