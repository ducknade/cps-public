//Hybrid Monte Carlo force term for an improved gauge action with plaquette
//and rectangle terms.  Takes advantage of lattice-wide operations to speed
//evolution of HMC trajectory.

#include <config.h>
#include <math.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/smalloc.h>
#include <util/pt_mat.h>
#include <util/time_cps.h>
#include <cassert>
CPS_START_NAMESPACE

#define USE_LINK_BUFFER 0
#if USE_LINK_BUFFER

ForceArg GimprRect::EvolveMomGforce(Matrix *mom, Float dt){
  char *fname = "EvolveMomGforce(M*,F)";
  VRB.Func(cname,fname);

  EnableLinkBuffer(0); VRB.Result(cname, fname, "Warning! Enabling link buffer!\n");
  assert(LinkBufferIsEnabled()); //needed for parallel for loop
  const int buffer_thickness = 2;
  printf("creating buffer\n");
  CreateBuffer(buffer_thickness, BUFFER_SHAPE_LOCAL_AND_FACES_AND_1ST_CORNERS);
  printf("finished creating buffer\n");

  Float L1=0.0;
  Float L2=0.0;
  Float Linf=0.0;

#ifdef PROFILE
  Float time = -dclock();
  ForceFlops = 0;
#endif
  
  const int vol_node_sites = GJP.VolNodeSites();
#pragma omp parallel num_threads(64) reduction(+:L1,L2)
  {
    Float Linf_priv = 0.0;
    Matrix mt0;
#pragma omp for 
    for(int uoff = 0; uoff < vol_node_sites; uoff++) {
      int x[4];
      CoordsFromOffset(x, uoff);

      for (int mu = 0; mu < 4; ++mu) {
        GforceSite(mt0, x, mu);

        IFloat *ihp = (IFloat *)(mom+uoff+mu);
        IFloat *dotp = (IFloat *)&mt0;
        fTimesV1PlusV2(ihp, dt, dotp, ihp, 18);

        Float norm = ((Matrix*)dotp)->norm();
        Float tmp = sqrt(norm);
        L1 += tmp;
        L2 += norm;
        if(tmp > Linf_priv) Linf_priv = tmp;
      }
    }

#pragma omp critical
    {
      if(Linf_priv > Linf) Linf = Linf_priv;
    }
  }

  ClearAllBufferedLink();

  ForceFlops +=GJP.VolNodeSites()*4*18*2;
#ifdef PROFILE
  time += dclock();
  print_flops(cname,fname,ForceFlops,time);
#endif

  glb_sum(&L1);
  glb_sum(&L2);
  glb_max(&Linf);

  L1 /= 4.0*GJP.VolSites();
  L2 /= 4.0*GJP.VolSites();

  VRB.FuncEnd(cname,fname);
  return ForceArg(dt*L1, dt*sqrt(L2), dt*Linf);
}

#else

//------------------------------------------------------------------
// EvolveMomGforce(Matrix *mom, Float dt):
// It evolves the canonical momentum mom by dt using the pure gauge
// force.
// ------------------------------------------------------------------
ForceArg GimprRect::EvolveMomGforce(Matrix *mom, Float dt)
{
    const char *fname = "EvolveMomGforce()";

    VRB.Result(cname, fname, "Hantao's code\n");

    const int vol = GJP.VolNodeSites();
    const int N = 4; //Num of dimensions

    //Pointer to block of unit matrices for each lattice site
    //Set all these matrices to the identity
    Matrix *Unit = (Matrix *) fmalloc(vol*sizeof(Matrix));
    for(int i = 0; i < vol; i++) Unit[i] = 1.;

    //Temporary matrices for use in calculation of the staples
    Matrix *tmp1[N];
    Matrix *tmp2[N];

    for(int i = 0;i<N;i++) {
        tmp1[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
        tmp2[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
        //Set all bytes in tmp2 to zero
        bzero((char*)tmp2[i],vol*sizeof(Matrix));
    }

    //Holds the sum of staples associated with each lattice site
    Matrix *result[N];
    for(int i = 0;i<N;i++) {
        result[i] = (Matrix *) fmalloc(vol*sizeof(Matrix));
    }

    //Array of four pointers to a block of unit matrices
    Matrix *Units[N] = {Unit, Unit, Unit, Unit};
    Matrix *gauge[N] = {GaugeField(), GaugeField(), GaugeField(), GaugeField()};

    {
        const int dirs_p[] = {0,2,4,6,0,2,4};   //Positive directions
        const int dirs_m[] = {1,3,5,7,1,3,5};   //Negative directions

        //Take a vector field V(x).  Suppose we wish to parallel transport
        //this vector field from the point x to the point x-mu.  That is,
        //we wish to find a new field V'(x-mu) in terms of V(x) such that
        //the new field V'(x-mu) changes appropriately under gauge transformation.
        //
        //The parallel transport that satisfies this property is
        //V'(x-mu) = U_mu(x-mu)V(x).  This combination transforms like
        //a vector field.
        //
        //When using the parallel transport class, one must specify an intial
        //field, a direction, and a final field.  However, the direction
        //that one must specify is the direction of the link matrix, not the
        //direction of the parallel transport.  Thus, to calculate V' as above,
        //use:
        //
        //pt.run(1,V',V,dir_Plus_mu)
        //
        //The new vector field will be indexed according to its new position.
        for(int nu = 1; nu < 4; ++nu) {

            //First calculate the staple in the positive nu direction
            pt_generic::pt(N, tmp1, gauge, Units, dirs_m + nu);
            pt_generic::pt(N, result, gauge, tmp1, dirs_m);
            pt_generic::pt(N, tmp1, gauge, result, dirs_p + nu);

            //tmp2 contains the sum of the staples for a given link
            for(int i = 0; i<N;i++) {
                tmp2[i]->FTimesV1PlusV2(plaq_coeff,tmp1[i],tmp2[i],vol);
            }

            //Calculating one rectangular staple
            pt_generic::pt(N,tmp1,gauge,Units,dirs_p);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m+nu);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m);
            pt_generic::pt(N,tmp1,gauge,result,dirs_p+nu);

            for(int i = 0; i<N;i++) {
                tmp2[i]->FTimesV1PlusV2(rect_coeff,tmp1[i],tmp2[i],vol);
            }

            //Calculating another rectangular staple;
            pt_generic::pt(N,tmp1,gauge,Units,dirs_m+nu);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m);
            pt_generic::pt(N,result,gauge,tmp1,dirs_p+nu);
            pt_generic::pt(N,tmp1,gauge,result,dirs_p);

            for(int i = 0; i<N;i++) {
                tmp2[i]->FTimesV1PlusV2(rect_coeff,tmp1[i],tmp2[i],vol);
            }

            //Calculating another rectangular staple;
            pt_generic::pt(N,tmp1,gauge,Units,dirs_m+nu);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m+nu);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m);
            pt_generic::pt(N,result,gauge,tmp1,dirs_p+nu);
            pt_generic::pt(N,tmp1,gauge,result,dirs_p+nu);

            for(int i = 0; i<N;i++) {
                tmp2[i]->FTimesV1PlusV2(rect_coeff,tmp1[i],tmp2[i],vol);
            }

            //Calculating the staple in the negative nu direction
            pt_generic::pt(N,tmp1,gauge,Units,dirs_p+nu);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m+nu);

            //Add this result into tmp2
            for(int i = 0; i<N;i++) {
                tmp2[i]->FTimesV1PlusV2(plaq_coeff,tmp1[i],tmp2[i],vol);
            }

            //Calculating one rectangular staple
            pt_generic::pt(N,tmp1,gauge,Units,dirs_p);
            pt_generic::pt(N,result,gauge,tmp1,dirs_p+nu);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m+nu);

            for(int i = 0; i<N;i++) {
                tmp2[i]->FTimesV1PlusV2(rect_coeff,tmp1[i],tmp2[i],vol);
            }

            //Calculating another rectangular staple;
            pt_generic::pt(N,tmp1,gauge,Units,dirs_p+nu);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m+nu);
            pt_generic::pt(N,tmp1,gauge,result,dirs_p);

            for(int i = 0; i<N;i++) {
                tmp2[i]->FTimesV1PlusV2(rect_coeff,tmp1[i],tmp2[i],vol);
            }

            //Calculating another rectangular staple;
            pt_generic::pt(N,tmp1,gauge,Units,dirs_p+nu);
            pt_generic::pt(N,result,gauge,tmp1,dirs_p+nu);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m);
            pt_generic::pt(N,result,gauge,tmp1,dirs_m+nu);
            pt_generic::pt(N,tmp1,gauge,result,dirs_m+nu);

            for(int i = 0; i<N;i++) {
                tmp2[i]->FTimesV1PlusV2(rect_coeff,tmp1[i],tmp2[i],vol);
            }
        }
        //Multiply on the left by our original link matrix to get force term
        pt_generic::pt(N,result,gauge,tmp2,dirs_p);
    }

    ForceArg ret;

#pragma omp parallel
    {
        ForceArg f_threaded;

#pragma omp for nowait
        for(int index = 0; index < 4 * vol; ++index) {
            int i  = index % vol;
            int mu = index / vol;
            
            Matrix mp1(result[mu][i]);
            mp1.TrLessAntiHermMatrix();
            mp1 *= dt;
            mom[mu + i * 4] += mp1;

            updateForce(f_threaded, mp1);
        }

#pragma omp critical
        {
            ret.combine(f_threaded);
        }
    }

    //Free some memory
    ffree(Unit);
    for(int i = 0;i<N;i++){
        ffree(tmp1[i]);
        ffree(tmp2[i]);
    }
    for(int i = 0;i<4;i++) 
        ffree(result[i]);

    ret.glb_reduce();
    return ret;
}

#endif

CPS_END_NAMESPACE
