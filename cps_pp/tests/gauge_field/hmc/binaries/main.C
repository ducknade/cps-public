#include <util/gauge_field.h>
#include <util/vector.h>

#include <util/time_cps.h>

#include <util/gjp.h>
#include <util/random.h>

#include <alg/do_arg.h>
#include <alg/common_arg.h>

#include <alg/plaq.h>
#include <alg/ghb.h>

#include <util/lattice.h>
#include <alg/alg_hmc.h>

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

void setDoArg(DoArg& do_arg, int * total_sites) {
  do_arg.x_sites = total_sites[0];
  do_arg.y_sites = total_sites[1];
  do_arg.z_sites = total_sites[2];
  do_arg.t_sites = total_sites[3];
  do_arg.s_sites = 32;
  do_arg.x_bc = BND_CND_PRD;
  do_arg.y_bc = BND_CND_PRD;
  do_arg.z_bc = BND_CND_PRD;
  do_arg.t_bc = BND_CND_PRD;
  do_arg.start_conf_kind = START_CONF_DISORD;
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
}

void runHmc(GaugeField& gf, GaugeAction& ga, int traj) {
  GJP.GetDoArg()->beta = ga.getBeta();
  if (G_CLASS_IMPR_RECT == ga.Gclass()) {
    GJP.GetDoArg()->c_1 = ga.getGaugeActionImprRect().getC1();
  }

  Lattice &lat1 = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
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

  CommonArg common_arg;
  HmcArg hmc_arg;
  hmc_arg.steps_per_traj = 6;
  hmc_arg.step_size = 1.0 / hmc_arg.steps_per_traj;
  hmc_arg.metropolis = METROPOLIS_YES;
  hmc_arg.reunitarize = REUNITARIZE_YES;
  hmc_arg.reverse = REVERSE_NO;
  hmc_arg.reproduce = REPRODUCE_NO;
  hmc_arg.reproduce_attempt_limit = 1;
  hmc_arg.wfm_md_sloppy = 0;
  AlgHmc hmc(ab1, common_arg, hmc_arg);
  hmc.run();

  AlgIntAB::Destroy(ab1);

  Lattice &lat2 = LatticeFactory::Create(F_CLASS_NONE, G_CLASS_NONE);
  gf.getLattice(lat2);
  LatticeFactory::Destroy();
}

int main(int argc, char *argv[]) {
  const char * fname = "main()";
  mkdir("../results", 0755);
  mkdir("../results/plaq", 0755);
  setlinebuf(stdout);

  Start(&argc, &argv);

//  int total_sites[] = { 8, 8, 8, 8 };
  int total_sites[] = { 16, 16, 16, 16 };
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

  for (int i = 0; i < 5; i++) {
    Float time = - dclock();
    gfConstruct3rdRow(gf);
    ghb.run(gf, ga);
    traj++;
    time += dclock();
    print_time(cname, fname, time);

    plaq.run(gf, addIntExtension("../results/plaq/plaq.", traj).c_str());
    VRB.Result(cname, fname, "plaq.getPlaqSum = %.16e ; Traj = %d\n", plaq.getPlaqSum(gf), traj);
  }

  for (int i = 0; i < 5; i++) {
    Float time = - dclock();
    gfConstruct3rdRow(gf);
    runHmc(gf, ga, traj);
    traj++;
    time += dclock();
    print_time(cname, fname, time);

    plaq.run(gf, addIntExtension("../results/plaq/plaq.", traj).c_str());
    VRB.Result(cname, fname, "plaq.getPlaqSum = %.16e ; Traj = %d\n", plaq.getPlaqSum(gf), traj);
  }

  End();
  return 0;
}
