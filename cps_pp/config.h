/*!\file
  \brief   Global options for compiling the Colombia code:

  Generated automatically from config.h.in by configure procedure.

  $Id: config.h.in,v 1.29 2012/03/26 13:50:11 chulwoo Exp $
*/
/* Global options for compiling the Columbia code:  
 * config.h.  Generated from config.h.in by configure.
 * 
 *--------------------------------------------------------------------
 *  CVS keywords
 *
 *  $Author: chulwoo $
 *  $Date: 2012/03/26 13:50:11 $
 *  $Header: /space/cvs/cps/cps++/config.h.in,v 1.29 2012/03/26 13:50:11 chulwoo Exp $
 *  $Id: config.h.in,v 1.29 2012/03/26 13:50:11 chulwoo Exp $
 *  $Name: v5_0_16_hantao_io_test_v7 $
 *  $Locker:  $
 *  $RCSfile: config.h.in,v $
 *  $Revision: 1.29 $
 *  $Source: /space/cvs/cps/cps++/config.h.in,v $
 *  $State: Exp $
 */
/* ------------------------------------------------------------------*/

#ifndef INCLUDED_CONFIG_H_
#define INCLUDED_CONFIG_H_                  //!< Prevent multiple inclusion 

#include <conf.h>

#define NOARCH 0
#define QCDOC  1
#define QCDSP  2
#define BGL    3
#define BGP    4
#define BGQ    5


#define VERSION_MAJOR 5
#define VERSION_MINOR 0
#define VERSION_SUB 16
#define VERSION_STR "CPS_V5.0.16"

#define TARGET NOARCH
#define PARALLEL 1

// The configure procedure should make this unnecessary, but just in case...
#ifndef TARGET
#define TARGET NOARCH
#endif

#if TARGET == QCDOC
#include<qalloc.h>
#include<qcdocos.h>
#endif

#if TARGET == BGL
#define CPS_FLOAT_ALIGN __attribute__((aligned(16)))
#else
#define CPS_FLOAT_ALIGN
#endif

#if TARGET == QCDOC
// temporary hack until qos is more mature
#define CWDPREFIX(A) "/"A
#else
#define CWDPREFIX(A) A
#endif



/*! Explicit casting away of the const-ness  */
#define CAST_AWAY_CONST(x) ( const_cast<char*>(x) )

/*!  Precision in the global sum (undefined gives QCDSP behaviour). */
#define GLOBALSUM_TYPE double

#define CPS_END_NAMESPACE    }  
#define CPS_START_NAMESPACE  namespace cps {
#define USING_NAMESPACE_CPS  using namespace cps;
#define CPS_NAMESPACE	     cps


#if TARGET == cpsMPI
/*! Data size for the MPI comms layer: */
#define COMMS_DATASIZE (sizeof(float))
/* Override printf to only print from only one processor */
#include<util/qcdio_qprintf.h>
#elif TARGET == BGL
/*! Data size for the MPI comms layer: */
#define COMMS_DATASIZE (sizeof(double))
/* Override printf to only print from only one processor */
#include<util/qcdio_qprintf.h>
#else
#define COMMS_DATASIZE (sizeof(double))
#endif

#undef UNIFORM_SEED_TESTING
#undef UNIFORM_SEED_NO_COMMS

/* ------------------------------------------------------------------*/

#endif /* INCLUDED_CONFIG_H_ */





