#ifndef ALG_WILSONFLOW
#define	ALG_WILSONFLOW
#include <config.h>
#include <alg/alg_base.h>
#include<util/site.h>

//using namespace cps;
CPS_START_NAMESPACE
class AlgWilsonFlow:public Alg
{
private:
	char *cname;
	bool su3_proj;//generally speaking it is not needed
	Float tolerance;
	Float dt;

	Float *Z_Lie;

	void logRun();

  void calculateZ(GaugeField &gf, const int pos[4], int mu, Float Z[8]);
  void DoRKStep(Lattice &lat, GaugeField& gf, int rk_step, int site); 
  
  void do_rk_step_adaptive(Lattice &lat, GaugeField& gf, GaugeField& gfp, int rk_step, int site); 

public:
	AlgWilsonFlow(Lattice& lat, CommonArg *ca, Float dtime=0.01, bool proj=true, Float tol=1e-8);
	virtual ~AlgWilsonFlow();

	void run();
	
  double run_adaptive(double input_dt);

	void su3projon(){su3_proj=true;}
	void su3projoff(){su3_proj=false;}
	void set_tol(Float x){tolerance=x;}
	Float get_tol() const {return tolerance;}
	void set_dt(Float dt){this->dt = dt;}
	Float get_dt() const {return dt;}
};

CPS_END_NAMESPACE

#endif

