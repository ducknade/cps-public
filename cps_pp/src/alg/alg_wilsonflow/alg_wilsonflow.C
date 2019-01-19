#include <alg/alg_wilsonflow.h>

#include <config.h>
#include <math.h>
#include <util/qcdio.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/gauge_field.h>
#include <util/lat_cont.h>

#include <mpi.h>

CPS_START_NAMESPACE

extern void generateSU3(Matrix &mat, Float Q[8]);
extern const Float SU3_lambda[8][18];

AlgWilsonFlow::AlgWilsonFlow(Lattice &lat, CommonArg *ca, Float dtime, bool proj, Float tol):
	Alg(lat,ca),su3_proj(proj),tolerance(tol),dt(dtime)
{
  cname = "AlgWilsonFlow";

  Z_Lie = new Float[GJP.VolNodeSites()*4*8];
  if(Z_Lie==NULL){ERR.Pointer(cname,cname, "Z_Lie");}
}

AlgWilsonFlow::~AlgWilsonFlow()
{
  delete [] Z_Lie;
}

void AlgWilsonFlow::logRun()
{
  const char* fname = "logRun";

  if(common_arg->filename != 0) {
    FILE* f = Fopen(common_arg->filename, "a");
    if(!f) ERR.FileA(cname, fname, common_arg->filename);
    Fprintf(f, "AlgWilsonFlow: dt = %e\n", dt);
    Fclose(f);
  } 
}

void AlgWilsonFlow::calculateZ(GaugeField &gf, const int pos[4], int mu, Float Z[8])
{
  Matrix stap;
  gfPlaqStaple(stap, gf, pos, mu);

  const Matrix *link = gf.localGetLink(pos,mu);
  Matrix loop;
  mDotMEqual((Float *)(&loop),(Float*)link, (Float*)(&stap));

  Float tmp[18];
  for(int i=0;i<8;i++)
  {
    mDotMEqual(tmp, SU3_lambda[i],(Float*)(&loop));
    Z[i]=-(tmp[1]+tmp[9]+tmp[17]); //-ImTr(tmp)
  }
}

void AlgWilsonFlow::DoRKStep(Lattice &lat, GaugeField& gf, int rk_step, int site) 
{
  const char* fname = "DoRKStep";
  //Follow Luscher's paper: arXiv:1006.4518v1 (June 2010)
  
  int x[4];
  gf.coordinatesFromIndex(x, site);
	
  Float QZ[8];
  Matrix expQ;
      
  for(int mu = 0; mu < 4; mu++) 
  {
    //QZ[a] = derivative of S with respect to changes in U_mu in the direction of the generator T^a
    calculateZ(gf, x, mu, QZ);      

    switch(rk_step) 
    {
      case 0:
        for(int i = 0; i < 8; i++) Z_Lie[8*(4*site+mu)+i] = QZ[i]/4.0 * dt;
        break;

      case 1:
        for(int i = 0; i < 8; i++) Z_Lie[8*(4*site+mu)+i] = QZ[i]*8.0/9.0*dt - 17.0/9.0*Z_Lie[8*(4*site+mu)+i];
        break;

      case 2:
	for(int i = 0; i < 8; i++) Z_Lie[8*(4*site+mu)+i] = QZ[i]*3.0/4.0*dt - Z_Lie[8*(4*site+mu)+i];
        break;
    }	

    generateSU3(expQ, &(Z_Lie[8*(4*site+mu)]));
      
    //new U_mu = expQ * U_mu
    mDotMEqual((Float *)(lat.GaugeField() + lat.GsiteOffset(x) + mu), 
               (Float *)&expQ, 
               (Float *)gf.localGetLink(x, mu));
  }
    
}

void AlgWilsonFlow::run()
{
  const char fname[] = "run()";

  logRun();

  Lattice& lat(AlgLattice());
    
  static bool initialized = false;
  static Offsets offsets;

  //Do each of the 3 steps of the RK integrator in turn:
  for(int rk_step = 0; rk_step < 3; rk_step++) 
  {
    cps::GaugeField gf(lat);
    const int buffer_thickness = 1; //Required expansion in each direction
    gf.resize(buffer_thickness);
    if (!initialized) {
      gf.refresh();
      gf.recordOffsetStart(false, true);
    } else {
      gf.refresh(offsets);
    }

#pragma omp parallel for
    for(int site = 0; site < GJP.VolNodeSites(); ++site) {
      DoRKStep(lat, gf, rk_step, site);
    }

    if (!initialized) {
      offsets = gf.getRecordedOffsets();
      gf.recordOffsetEnd();
      initialized = true;
    }
  }
  
  //Check if smearing messes up open boundary conditions (it shouldn't!)
  if(GJP.TopenBc() && !lat.TestZeroTboundary()) ERR.General(cname, fname, "Nonzero boundary link!\n");
}

double Linf_distance(Lattice& lat, const GaugeField& gfp){
  double local_max = 0.;
  for(int site = 0; site < GJP.VolNodeSites(); ++site){ 
    int x[4];
    gfp.coordinatesFromIndex(x, site);
    for(int mu = 0; mu < 4; mu++){
      double* pA = reinterpret_cast<double*>( lat.GaugeField()+lat.GsiteOffset(x)+mu );
      // double* pA = reinterpret_cast<double*>(gf.localGetLink(x, mu));
      double* pB = reinterpret_cast<double*>( gfp.localGetLink(x, mu) );
      double sum = 0.;
      for(int idx = 0; idx < 18; idx++){
        // if(site == 4 && cps::UniqueID() == 0){
        //   printf("pA = %+.8f pB = %+.8f\n", pA[idx], pB[idx]);
        // }
        sum += (pA[idx]-pB[idx])*(pA[idx]-pB[idx]);
      }
      if(sum > local_max){
        local_max = sum;
      }
    }
  }
  local_max = sqrt(local_max/18.);
  printf("node #%03d: local_max =  %.12f\n", cps::UniqueID(), local_max);
  double global_max;
  MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(cps::UniqueID() == 0) printf("Linf distance: global_max = %.12f\n", global_max);
  return global_max;
} 

void AlgWilsonFlow::do_rk_step_adaptive(Lattice &lat, GaugeField& gf, GaugeField& gfp, int rk_step, int site) 
{
  const char* fname = "DoRKStep";
  //Follow Luscher's paper: arXiv:1006.4518v1 (June 2010)
  
  int x[4];
  gf.coordinatesFromIndex(x, site);
	
  Float QZ[8];
  Matrix expQ;
      
  for(int mu = 0; mu < 4; mu++) 
  {
    //QZ[a] = derivative of S with respect to changes in U_mu in the direction of the generator T^a
    calculateZ(gf, x, mu, QZ);      

    Float z_lie_p[8]; 
    Matrix expQp, cp;

    switch(rk_step) 
    {
      case 0:
        for(int i = 0; i < 8; i++) Z_Lie[8*(4*site+mu)+i] = QZ[i]/4.0 * dt;
        break;

      case 1:
        for(int i = 0; i < 8; i++){
          z_lie_p[i] = 2.*QZ[i]*dt-4.*Z_Lie[8*(4*site+mu)+i];
          Z_Lie[8*(4*site+mu)+i] = QZ[i]*8.0/9.0*dt - 17.0/9.0*Z_Lie[8*(4*site+mu)+i];
        }
        generateSU3(expQp, z_lie_p);
        cp = expQp * (*gfp.localGetLink(x, mu));
        (*gfp.localGetLink(x, mu)) = cp;
        break;
      case 2:
	      for(int i = 0; i < 8; i++) Z_Lie[8*(4*site+mu)+i] = QZ[i]*3.0/4.0*dt - Z_Lie[8*(4*site+mu)+i];
        break;
    }	

    generateSU3(expQ, &(Z_Lie[8*(4*site+mu)]));
      
    //new U_mu = expQ * U_mu
    mDotMEqual((Float *)(lat.GaugeField() + lat.GsiteOffset(x) + mu), 
               (Float *)&expQ, 
               (Float *)gf.localGetLink(x, mu));
  }
    
}

double AlgWilsonFlow::run_adaptive(double input_dt)
{
  const char fname[] = "run()";

  static const double tolerance = 1e-3;
  static const double multiplier = 0.95;

  Lattice& lat(AlgLattice());
    
  static bool initialized = false;
  static Offsets offsets;

  const int buffer_thickness = 1; //Required expansion in each direction

  LatticeContainer stash;
  stash.Get(lat);

  dt = input_dt;

  while(true){
    stash.Set(lat);
    cps::GaugeField gfp(lat);
    gfp.resize(buffer_thickness);

    //Do each of the 3 steps of the RK integrator in turn:
    for(int rk_step = 0; rk_step < 3; rk_step++) 
    {
      cps::GaugeField gf(lat);
      gf.resize(buffer_thickness);
      if (!initialized) {
        gf.refresh();
        gf.recordOffsetStart(false, true);
      } else {
        gf.refresh(offsets);
      }

#pragma omp parallel for
      for(int site = 0; site < GJP.VolNodeSites(); ++site) {
        do_rk_step_adaptive(lat, gf, gfp, rk_step, site);
      }

      if (!initialized) {
        offsets = gf.getRecordedOffsets();
        gf.recordOffsetEnd();
        initialized = true;
      }
    }

    double d = Linf_distance(lat, gfp);
    double previous_dt = dt;
    if(d < tolerance){
      if(cps::UniqueID() == 0) printf("Adaptive Wilson flow ACCEPTED  with d = %12.8e, tolerance = %12.8e.\n", d, tolerance);
      logRun();
      dt *= multiplier*std::pow(tolerance/d, 1./3.);
      if(cps::UniqueID() == 0) printf("dt: %12.8e --> %12.8e.\n", previous_dt, dt);
      break;
    }else{
      if(cps::UniqueID() == 0) printf("Adaptive Wilson flow DISCARDED with d = %12.8e, tolerance = %12.8e.\n", d, tolerance);
      dt *= multiplier*std::pow(tolerance/d, 1./3.);
      if(cps::UniqueID() == 0) printf("dt: %12.8e --> %12.8e.\n", previous_dt, dt);
    }
  }

  //Check if smearing messes up open boundary conditions (it shouldn't!)
  if(GJP.TopenBc() && !lat.TestZeroTboundary()) ERR.General(cname, fname, "Nonzero boundary link!\n");

  return dt;
}


CPS_END_NAMESPACE

