#include <util/vector.h>

#include <util/time_cps.h>

#include <util/gjp.h>
#include <util/random.h>

#include <alg/do_arg.h>
#include <alg/common_arg.h>

#include <alg/alg_wilsonflow.h>
#include <alg/alg_actiondensity.h>

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

int main(int argc, char *argv[]) {
  const char * fname = "main()";
  mkdir("../results", 0755);
  setlinebuf(stdout);

  Start(&argc, &argv);

  int total_sites[] = {16, 16, 16, 32};
  DoArg do_arg;
  setDoArg(do_arg, total_sites);

  GJP.Initialize(do_arg);
  LRG.Initialize();


  GwilsonFnone lat;

  CommonArg common_arg;
  common_arg.set_filename("../results/wflow.out");

  Float dt = 0.05;
  AlgWilsonFlow wflow(lat, &common_arg, dt);
  AlgActionDensity action(lat, &common_arg);

  for(int i = 0; i < 60; i++) {
    VRB.Result(cname, fname, "Running round %d\n", i);

    Float timer = -dclock();
    action.run();
    timer += dclock();
    VRB.Result(cname, fname, "AlgActionDensity round %d: %e seconds\n", i, timer);

    timer = -dclock();
    wflow.run();
    timer += dclock();
    VRB.Result(cname, fname, "AlgWilsonFlow round %d: %e seconds\n", i, timer);
  }
  action.run();

  End();
  return 0;
}
