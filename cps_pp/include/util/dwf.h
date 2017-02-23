#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008/04/21 14:19:17 $
//  $Header: /space/cvs/cps/cps++/include/util/dwf.h,v 1.6 2008/04/21 14:19:17 chulwoo Exp $
//  $Id: dwf.h,v 1.6 2008/04/21 14:19:17 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: dwf.h,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/include/util/dwf.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf.h
//
// C header file for the dwf fermion library
//
//------------------------------------------------------------------

#ifndef INCLUDED_DWF_H
#define INCLUDED_DWF_H

CPS_END_NAMESPACE
#include <util/wilson.h>
#include <util/vector.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// External declarations
//------------------------------------------------------------------

extern int dwfso_wire_map[];
// Set in dwf_int. For a given index 0-1 corresponding to
// S+, S-  it gives the corresponding wire.


//------------------------------------------------------------------
// Type definitions
//------------------------------------------------------------------
// The Dwf structure typedef
typedef struct{
  IFloat dwf_kappa;    // 1/[2*(5-dwf_height)]
  int vol_4d;         // The 4-dimensional volume   
  int ls;             // The extent of the 5th direction
  IFloat *frm_tmp1;    // Pointer to temorary fermion field 1
  IFloat *frm_tmp2;    // Pointer to temorary fermion field 2
  IFloat *comm_buf;    // Communications buffer (12 words for dwfso)
  Wilson *wilson_p;   // Pointer to the wilson structure
#if TARGET == QCDOC
  SCUDirArgIR *PlusArg[2];
  SCUDirArgIR *MinusArg[2];
  SCUDirArgMulti *Plus;
  SCUDirArgMulti *Minus;
#endif
} Dwf;


//------------------------------------------------------------------
// Function declarations
//------------------------------------------------------------------


//------------------------------------------------------------------
// This routine performs all initializations needed before dwf 
// library funcs are called. It sets the addressing related arrays 
// and reserves memory for the needed temporary buffers. It only 
// needs to be called only once at the begining of the program 
// (or after a dwf_end call)before any number of calls to dwf 
// functions are made.
//------------------------------------------------------------------
void dwf_init(Dwf *dwf_p);             // pointer to Dwf struct.


//------------------------------------------------------------------
// This routine frees any memory reserved by dwf_init
//------------------------------------------------------------------
void dwf_end(Dwf *dwf_p);              // pointer to Dwf struct.


//------------------------------------------------------------------
// dwf_mdagm M^dag M where M is the fermion matrix.
// The in, out fields are defined on the checkerboard lattice.
// <out, in> = <dwf_mdagm*in, in> is returned in dot_prd.
//------------------------------------------------------------------
void  dwf_mdagm(Vector *out, 
		Matrix *gauge_field, 
		Vector *in, 
		Float *dot_prd,
		Float mass,
		Dwf *dwf_lib_arg);


//------------------------------------------------------------------
// dwf_dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash(Vector *out, 
		Matrix *gauge_field, 
		Vector *in,
		Float mass,
		int cb,
		int dag,
		Dwf *dwf_lib_arg);


//------------------------------------------------------------------
// dwf_m is the fermion matrix.  
// The in, out fields are defined on the checkerboard lattice.
//------------------------------------------------------------------
void  dwf_m(Vector *out, 
	    Matrix *gauge_field, 
	    Vector *in,
	    Float mass,
	    Dwf *dwf_lib_arg);


//------------------------------------------------------------------
// dwf_mdag is the dagger of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
//------------------------------------------------------------------
void dwf_mdag(Vector *out, 
	      Matrix *gauge_field, 
	      Vector *in,
	      Float mass,
	      Dwf *dwf_lib_arg);



//------------------------------------------------------------------
// dwf_dslash_4 is the derivative part of the space-time part of
// the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
// cb = 0/1 <--> even/odd checkerboard of in field.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash_4(Vector *out, 
		  Matrix *gauge_field, 
		  Vector *in,
		  int cb,
		  int dag,
		  Dwf *dwf_lib_arg);



//------------------------------------------------------------------
// dwf_dslash_5_plus is the derivative part of the 5th direction
// part of the fermion matrix. This routine accumulates the result
// on the out field 
// The in, out fields are defined on the checkerboard lattice.
// The action of this operator is the same for even/odd
// checkerboard fields because there is no gauge field along
// the 5th direction.
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//------------------------------------------------------------------
void dwf_dslash_5_plus(Vector *out, 
		       Vector *in,
		       Float mass,
		       int dag,
		       Dwf *dwf_lib_arg);
// subroutine for possible better cache management  (CJ)
void dwf_dslash_5_plus_slice(Vector *out, 
		       Vector *in,
		       Float mass,
		       int dag,
		       Dwf *dwf_lib_arg,
     		       int s_slice);
void dwf_dslash_5_plus_start(Vector *out, 
		       Vector *in,
		       Float mass,
		       int dag,
		       Dwf *dwf_lib_arg);
void dwf_dslash_all(Vector *out, 
		Matrix *gauge_field, 
		Vector *in,
		Float mass,
		int cb,
		int dag,
		Dwf *dwf_lib_arg);

#endif

CPS_END_NAMESPACE