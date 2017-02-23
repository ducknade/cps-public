#include <config.h>
#include <util/gjp.h>
#include <util/lattice.h>
#include <util/time_cps.h>
#include <util/gauge_field.h>
#include <comms/nga_reg.h>
#include <comms/glb.h>
#include <comms/cbuf.h>
CPS_START_NAMESPACE

ForceArg Gwilson::EvolveMomGforce(Matrix *mom, Float dt) {
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

CPS_END_NAMESPACE
