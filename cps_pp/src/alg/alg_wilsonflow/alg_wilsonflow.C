#include <alg/alg_wilsonflow.h>

#include <config.h>
#include <math.h>
#include <util/qcdio.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/error.h>
#include <util/gauge_field.h>


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


CPS_END_NAMESPACE

