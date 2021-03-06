#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTransStagTypes class constructor and destructor.

  $Id: pt_gauge.C,v 1.4 2011/02/26 00:19:27 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2011/02/26 00:19:27 $
//  $Header: /space/cvs/cps/cps++/src/util/parallel_transport/pt_gauge/pt_gauge.C,v 1.4 2011/02/26 00:19:27 chulwoo Exp $
//  $Id: pt_gauge.C,v 1.4 2011/02/26 00:19:27 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: pt_gauge.C,v $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/src/util/parallel_transport/pt_gauge/pt_gauge.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/pt.h>
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
/*!
  \param latt The lattice on which this operation is defined
 */
//------------------------------------------------------------------

static StrOrdType old_str_ord;
ParTransGauge::ParTransGauge(Lattice & latt) :
                                   ParTrans(latt)
{
  cname = "ParTransGauge";
  char *fname = "ParTransGauge(Lattice&)";
  VRB.Func(cname,fname);
//  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL ){
  if (lat.StrOrd() != CANONICAL ){
    old_str_ord = lat.StrOrd();
    lat.Convert(CANONICAL);
  }
  pt_init(lat);
  pt_init_g();
}


//------------------------------------------------------------------
ParTransGauge::~ParTransGauge() {
  char *fname = "~ParTransGauge()";
  VRB.Func(cname,fname);
  pt_delete_g();
  pt_delete();
//  if (old_str_ord != WILSON && old_str_ord != CANONICAL ){
  if (old_str_ord != CANONICAL ){
    lat.Convert(old_str_ord);
  }
}

CPS_END_NAMESPACE
