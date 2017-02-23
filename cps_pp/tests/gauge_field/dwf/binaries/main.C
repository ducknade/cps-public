#include <util/vector.h>
#include <util/gauge_field.h>
#include "fermion_action.h"

#include <util/time_cps.h>

#include <util/gjp.h>
#include <util/random.h>

#include <alg/do_arg.h>
#include <alg/common_arg.h>

#include <alg/plaq.h>
#include <alg/ghb.h>

#include <util/lattice.h>
#include <alg/alg_hmc.h>
#include <alg/alg_fix_gauge.h>
#include <alg/qpropw.h>

#include <iostream>
#include <sstream>
#include <string>

#include <sys/stat.h>
#include <stdio.h>
#include <omp.h>

using namespace std;
using namespace cps;

const char * cname = "Global";

string addIntExtension(const string& prefix, int ext) {
  ostringstream out;
  out << prefix << ext;
  return out.str();
}

void setDoArg(DoArg& do_arg, int *total_sites) {
  do_arg.x_sites = total_sites[0];
  do_arg.y_sites = total_sites[1];
  do_arg.z_sites = total_sites[2];
  do_arg.t_sites = total_sites[3];
  do_arg.s_sites = 8;
  do_arg.dwf_height = 1.8;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_APRD;
  do_arg.start_conf_kind = START_CONF_ORD;
  do_arg.start_seed_kind = START_SEED_INPUT;
  do_arg.start_seed_value = 123121;
  do_arg.x_nodes = 0;
  do_arg.y_nodes = 0;
  do_arg.z_nodes = 0;
  do_arg.t_nodes = 0;
  do_arg.s_nodes = 0;
  do_arg.x_node_sites = 0;
  do_arg.y_node_sites = 0;
  do_arg.z_node_sites = 0;
  do_arg.t_node_sites = 0;
  do_arg.s_node_sites = 0;
  do_arg.gfix_chkb = 1;
}

void runHmc(GaugeField& gf, GaugeAction& ga) {
  GJP.GetDoArg()->beta = ga.getBeta();
  if (G_CLASS_IMPR_RECT == ga.Gclass()) {
    GJP.GetDoArg()->c_1 = ga.getGaugeActionImprRect().getC1();
  }

  Lattice& lat1 = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  gf.setLattice(lat1);
  LatticeFactory::Destroy();

  AlgMomentum mom;
  ActionGaugeArg gauge_arg;
  gauge_arg.gluon = ga.Gclass();
  gauge_arg.action_arg.force_measure = FORCE_MEASURE_NO;
  gauge_arg.action_arg.force_label= "Gauge";
  AlgActionGauge gauge(mom, gauge_arg);

  IntABArg ab1_arg;
  ab1_arg.type = INT_FORCE_GRAD_QPQPQ;
  ab1_arg.A_steps = 1;
  ab1_arg.B_steps = 1;
  ab1_arg.level = TOP_LEVEL_INTEGRATOR;
  ab1_arg.lambda = 0.22;
  AlgIntAB &ab1 = AlgIntAB::Create(mom, gauge, ab1_arg);

  CommonArg c_arg;
  HmcArg hmc_arg;
  hmc_arg.steps_per_traj = 6;
  hmc_arg.step_size = 1.0 / hmc_arg.steps_per_traj;
  //hmc_arg.metropolis = METROPOLIS_YES;
  hmc_arg.metropolis = METROPOLIS_NO;
  hmc_arg.reunitarize = REUNITARIZE_YES;
  hmc_arg.reverse = REVERSE_NO;
  hmc_arg.reproduce = REPRODUCE_NO;
  hmc_arg.reproduce_attempt_limit = 1;
  hmc_arg.wfm_md_sloppy = 0;
  AlgHmc hmc(ab1, c_arg, hmc_arg);
  hmc.run();

  AlgIntAB::Destroy(ab1);

  Lattice& lat2 = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  gf.getLattice(lat2);
  LatticeFactory::Destroy();
}

void runFixGauge(Lattice& lat) {
  const char * fname = "runFixGauge(lat)";
  FixGaugeArg fix_gauge_arg;
  fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_LANDAU;
  //fix_gauge_arg.fix_gauge_kind = FIX_GAUGE_COULOMB_T;
  //fix_gauge_arg.hyperplane_start = 0;
  //fix_gauge_arg.hyperplane_step = 1;
  //fix_gauge_arg.hyperplane_num = 1;
  fix_gauge_arg.stop_cond = 1.0e-12;
  fix_gauge_arg.max_iter_num = 1000;

  if (NULL != lat.FixGaugePtr()) {
    VRB.Result(cname, fname, "Unexpected FixGaugeFree.\n");
    lat.FixGaugeFree();
  }
  CommonArg c_arg;
  AlgFixGauge fg(lat, &c_arg, &fix_gauge_arg);
  fg.run();
}

double runQPropW(GaugeField& gf) {
  const char * fname = "runQPropW(gf)";
  Lattice& lat = LatticeFactory::Create(F_CLASS_DWF, G_CLASS_NONE);
  gf.setLattice(lat);
  runFixGauge(lat);

  CgArg cg_arg;
  cg_arg.mass = 0.1;
  cg_arg.epsilon = 0;
  cg_arg.max_num_iter = 1;
  cg_arg.stop_rsd = 1.0e-8;
  cg_arg.true_rsd = 1.0e-8;
  cg_arg.RitzMatOper = MAT_HERM;
  cg_arg.Inverter = CG;
  cg_arg.bicgstab_n = 0;

  QPropWArg qpropw_arg;
  qpropw_arg.cg = cg_arg;
  qpropw_arg.file = "";
  qpropw_arg.x = 0;
  qpropw_arg.y = 0;
  qpropw_arg.z = 0;
  qpropw_arg.t = 0;
  qpropw_arg.gauge_fix_src = 1;
  qpropw_arg.gauge_fix_snk = 1;
  qpropw_arg.store_midprop = 0;
  qpropw_arg.save_prop = 0;
  qpropw_arg.save_ls_prop = 0;
  qpropw_arg.do_half_fermion = 0;
  qpropw_arg.SeqSmearSink = WALL;
  qpropw_arg.ensemble_label = "lquark";
  qpropw_arg.ensemble_id = "0001";
  qpropw_arg.seqNum = 0;
  qpropw_arg.StartSrcSpin = 0;
  qpropw_arg.EndSrcSpin = 4;
  qpropw_arg.StartSrcColor = 0;
  qpropw_arg.EndSrcColor = 3;

  CommonArg c_arg;
  int mom[] = { 0, 0, 0 };
  QPropWMomSrc qpropw(lat, &qpropw_arg, mom, &c_arg);
//  QPropWPointSrc qpropw(lat, &qpropw_arg, &c_arg);
  VRB.Result(cname, fname, "qpropw result = %.16e %.16e\n",
      qpropw.MomSinkProp(1,mom)(0,0,0,0).real(), qpropw.MomSinkProp(1,mom)(0,0,0,0).imag());
  if (NULL != lat.FixGaugePtr()) {
    lat.FixGaugeFree();
  }
  LatticeFactory::Destroy();
}

void init_bfm(int *argc, char **argv[]);

int main(int argc, char *argv[]) {
  const char * fname = "main()";
  mkdir("../results", 0755);
  mkdir("../results/plaq", 0755);
  setlinebuf(stdout);

  Start(&argc, &argv);

  int total_sites[] = { 8, 8, 8, 16 };
  DoArg do_arg;
  setDoArg(do_arg, total_sites);

  GJP.Initialize(do_arg);
  LRG.Initialize();

  GaugeField gf(total_sites);
  gfUnitMatrix(gf);

  GaugeActionImprRect ga(8.3, 0.05);
  gf.resize(gaInteractionRange(ga));

  Plaq plaq;

  GhbArg ghb_arg;
  ghb_arg.num_iter = 10;
  GheatBath ghb(ghb_arg);

  int traj = 0;

  plaq.run(gf, addIntExtension("../results/plaq/plaq.", traj).c_str());

  VRB.Result(cname, fname, "plaq.getPlaqSum = %.16e ; Traj = %d\n", plaq.getPlaqSum(gf), traj);
  runQPropW(gf);

  traj++;
  runHmc(gf, ga);
  plaq.run(gf, addIntExtension("../results/plaq/plaq.", traj).c_str());

  VRB.Result(cname, fname, "plaq.getPlaqSum = %.16e ; Traj = %d\n", plaq.getPlaqSum(gf), traj);
  runQPropW(gf);

  traj++;
  runHmc(gf, ga);
  plaq.run(gf, addIntExtension("../results/plaq/plaq.", traj).c_str());

  VRB.Result(cname, fname, "plaq.getPlaqSum = %.16e ; Traj = %d\n", plaq.getPlaqSum(gf), traj);
  runQPropW(gf);

  End();
  return 0;
}
