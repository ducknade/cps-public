#include<config.h>
CPS_START_NAMESPACE
/*----------------------------------------------------------*/
/*!\file
  \brief  Definitions for the MPI implementation of the QCDSP/QCDOC communications layer.
  
  $Id: sysfunc.C,v 1.10 2008/02/08 18:35:06 chulwoo Exp $
*/
/*----------------------------------------------------------------------
/* The Sysfunc Comms Interface: sysfunc.C

  The MPI implementation of the QCDSP SCU comms-layer.

  A.N.Jackson: ajackson@epcc.ed.ac.uk                      
  -----------------------------------------------------------
  CVS keywords
 
  $Author: chulwoo $
  $Date: 2008/02/08 18:35:06 $
  $Header: /space/cvs/cps/cps++/src/comms/mpi/scu/sysfunc.C,v 1.10 2008/02/08 18:35:06 chulwoo Exp $
  $Id: sysfunc.C,v 1.10 2008/02/08 18:35:06 chulwoo Exp $
  $Name: v5_0_16_hantao_io_test_v7 $
  $Locker:  $
  $RCSfile: sysfunc.C,v $
  $Revision: 1.10 $
  $Source: /space/cvs/cps/cps++/src/comms/mpi/scu/sysfunc.C,v $
  $State: Exp $  */
/*----------------------------------------------------------*/

CPS_END_NAMESPACE
#include<comms/sysfunc_cps.h>
#include <util/qcdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
CPS_START_NAMESPACE

/*!\namespace cps
  \brief Main namespace for all CPS classes, functions <em>etc</em>.
*/
/*!\namespace cps::MPISCU
  \brief Namespace for the MPI emulations of the SCU.
*/


// File-scoped data used by MPISCU functions.

namespace MPISCU{

// Emulation layer things

    bool Is_Initialised = false;
//!< Whether the MPI-SCU layer has been initialised.
/*!<
  The initial value of the global comms-layer
  initialization flag is FALSE, indicating that the communications system has
  not been set up.
*/

//! Name of the environment variable with the parallel execution parameters.
/*! This environment variable might define the parameters directly or
  name a file which does.
*/
    static char *ENVVAR = "COMMS_DEF";

//! Default name for the file containing the parallel execution parameters.
    static char *default_filename = "commsMPI.def";

    static const int STRING_MAX_LEN = 10000;
    
    static char  logFileName[STRING_MAX_LEN];
/*!< String containing the name of the file to be used for verbose output. */

    static FILE *logFile;      /*!< Pointer to the logfile output stream. */
    static unsigned int RNGseed;       /*!< Seed for the RNG */
    static char         seedFileName[STRING_MAX_LEN];
/*!< Filename of the file that holds the list of RNG seeds. */


// MPI things

/*!
  Definition of the default size (in bytes) for the basic elements (ints and
  IFloats) to be communicated .
*/ 
//    static int Datasize = COMMS_DATASIZE;

/*! The total number of MPI int types */
    static const int N_INT_TYPES =   4;                  

/*! The total number of MPI IFloat types */
    static const int N_FLOAT_TYPES = 3;                  

    static MPI_Comm Cart_Comm;          /*!< MPI communicator for the  
				     cartesian PE topology. */
    static MPIRequestManager   *ReqMan;   /*!< Pointer to the MPI
					    request manager */
    static MPI_Datatype *mpi_dt;        /*!< Current MPI datatype */
    

// Grid geometry

//! Number of grid dimensions.
    static const int NDIM = 5;

    static int peGrid[NDIM]; // initialise to invalid value.
//!< Number of processors in each direction.

    static int peRank;        /*!< Rank/identify of this  process */
    static int peNum;          /*!< Total number of processors */
    static int pePos[NDIM];  /*!< Position of this process in the grid.*/ 
    static int root_pe;        /*!< Specify the root processor by rank */ 

    static int          nnList[2*NDIM];/*!< Look-up table of NN PEs, 
						indexed by SCUDIR. */

} //namespace


/*-------------------------------------------------------------------------*/
/* Definitions of the actual emulated SCU system functions                 */
/*-------------------------------------------------------------------------*/

int UniqueID() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::peRank+1; }

int CoorT() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::pePos[SCU_T]; }
int CoorX() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::pePos[SCU_X]; }
int CoorY() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::pePos[SCU_Y]; }
int CoorZ() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::pePos[SCU_Z]; }
int CoorS() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::pePos[SCU_S]; }

int SizeT() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::peGrid[SCU_T]; }
int SizeX() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::peGrid[SCU_X]; }
int SizeY() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::peGrid[SCU_Y]; }
int SizeZ() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::peGrid[SCU_Z]; }
int SizeS() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::peGrid[SCU_S]; }

int NumNodes() { if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::peNum; }

//----------------------------------------------------------------
/*
  The seed can be different for each node and can
  change every time the machine is reset.

  \note The behaviour of the MPI and serial implementations
  may differ from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int Seed(){ return MPISCU::ReadSeedFile(); }
//----------------------------------------------------------------
/*
  The seed is the same for each node (spatially fixed, hence the S), but
  can change every time the machine is reset.

  \note The behaviour of the MPI and serial implementations may differ
  from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedS(){ if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); return MPISCU::RNGseed; }
//----------------------------------------------------------------
/*
  SeedT can be different for each node, but is fixed in time (the T), so it is
  unchanged by a reset.

  \note The behaviour of the MPI and serial implementations may
  differ from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedT(){ return Seed(); }
//----------------------------------------------------------------
/*
  SeedST is the same for each node (spatially fixed, hence the S), and the
  same after every reset (fixed time, hence T).

  \note The behaviour of the MPI and serial implementations may differ
  from that of the QCDSP/QCDOC system function.
*/
//----------------------------------------------------------------
unsigned int SeedST(){ return SeedS(); }

//----------------------------------------------------------------
/*
  This function blocks further code execution until all
  nodes in the machine have begun executing the code in the sync()
  routine.
  \return 0
*/
//----------------------------------------------------------------
unsigned int sync() {
    if( !MPISCU::Is_Initialised ) MPISCU::CommsInit(); 
    MPI_Barrier( MPISCU::Cart_Comm );
    return 0;
}

//----------------------------------------------------------------
/*
  On QCDSP this function returns the explicit wire
  number (0 - 7) of the physics direction given by \a dir. In the MPI
  version this returns the internal direction from the cartesian
  communicator which corresponds to the given physics direction.
  \param dir The physics (lattice) direction.
  \return The number used by the comms layer to represents that direction.

  Possibly.
*/
/* In this implementation, this just returns the integer value
  associated with the direction from the SCUDir enum */
int SCURemap( SCUDir dir ) {
    return (int)dir;
}


/*! SCUTrans (multiple, overloaded):
  N.B. For each direction, send and recieve requests must be in the
  same order!  */
/*!
  Performs the communication specified in \a arg.
  \param arg The object fully specifying what data to send to (or receive from)
  where.
*/
//----------------------------------------------------------------
void SCUTrans( SCUDirArg * arg ) {
    SCUTrans( &arg, 1 );
}

//----------------------------------------------------------------
/*!
  Performs the multiple communications specified in \a arg.
  \param arg A pointer to an array of objects, each fully specifying what
  data to send to (or receive from) where.
*/
//----------------------------------------------------------------
void SCUTrans( SCUDirArg ** arg, int n ) {

    for(int i=0; i<n; i++ ) 
	MPISCU::Trans( arg[i]->Addr(), 
		      arg[i]->Datatype(), 
		      arg[i]->CommDir(), 
		      arg[i]->CommType() );

}

//----------------------------------------------------------------
/*!
  This function does a number of transfers in the same direction of data
  with the same block length, stride and number
  of blocks but with different addresses. These addresses are specified as
  offsets to the base address.
  \param arg The object containing the information about the structure of the
  data, its base address, \e etc.
  \param offset The array of offsets from the base address, defining addresses
  to/from where data will be sent/received.
  \param n The number of data transfers.
*/
//----------------------------------------------------------------
void SCUTrans( SCUDirArg * arg, unsigned int * offset, int n ) {

    for(int i=0; i<n; i++ ) 
	MPISCU::Trans( (void*)((unsigned int)arg->Addr() + offset[i]), 
		      arg->Datatype(), 
		      arg->CommDir(), 
		      arg->CommType() );

}

//----------------------------------------------------------------
/*!
  The block length, stride and number of blocks involved in the data
  transfer are given to to the underlying communications layer,
  but no transfers are done.
  \param arg The object containing the information about the structure of the
  data.
*/
//----------------------------------------------------------------
void SCUSetDMA( SCUDirArg * arg ) { SCUSetDMA( &arg, 1 ); }

//----------------------------------------------------------------
/*!
  The block length, stride and number of blocks involved in a number of data
  transfer are given to to the underlying communications layer,
  but no transfers are done.
  \param arg A pointer ot an array of objects containing the information about
  the structure of the data.
  \param n The number of sets of data to be transferred.
//----------------------------------------------------------------*/
void SCUSetDMA( SCUDirArg ** arg, int n ) {

    // remove old settings if necessary:
    if( MPISCU::mpi_dt != NULL ) delete[] MPISCU::mpi_dt;

    // Remember the set of datatypes:
    MPISCU::mpi_dt = new MPI_Datatype[n];
    for(int  i=0; i<n; i++ ) MPISCU::mpi_dt[i] = arg[i]->Datatype();

}

//----------------------------------------------------------------
/*!
  Performs the communication specified by its arguments.

  \pre The transfer must have been set up using ::SCUTransAddr

  \param arg The object specifiying the base address of the data, the
  direction of the transfer and whether to send of receive the data. 
*/
//----------------------------------------------------------------
void SCUTransAddr( SCUDirArg * arg ) { SCUTransAddr( &arg, 1 ); }

//----------------------------------------------------------------
/*!
  Performs the communications specified by its arguments.

  \pre The transfers must have been set up using ::SCUTransAddr

  \param arg A pointer to an array of objects specifiying, for each transfer,
  the base address of the data, the direction of the transfer and whether to
  send of receive the data.
  \param n The number of sets of data to be transferred.
*/  
//----------------------------------------------------------------
void SCUTransAddr( SCUDirArg ** arg, int n ) {

    for(int i=0; i<n; i++ ) 
	MPISCU::Trans( arg[i]->Addr(), 
		      MPISCU::mpi_dt[i], 
		      arg[i]->CommDir(), 
		      arg[i]->CommType() );
    
}


//----------------------------------------------------------------
/*!
  This function returns only when all pending communications operations
  have completed.
*/
//----------------------------------------------------------------
void SCUTransComplete( void ) {

    

    // Check how many requests are pending:
    int numreq = MPISCU::ReqMan->NumReq();
    
    if( numreq > 0 ) {
#ifdef MPISCU_DEBUG
	MPISCU::fprintf_all(MPISCU::logFile,
			    "SCUTransComplete: numreq = %i\n", numreq ); 
#endif
	/* Grab the requires status array */
	MPI_Status *stat = new MPI_Status[numreq];
	/* WAIT for all the previously initiated communications to finish */
	MPI_Waitall( numreq, MPISCU::ReqMan->ReqArray(), stat );
	/* Throw away the status */
	delete[] stat;
#ifdef MPISCU_DEBUG
	MPISCU::fprintf_all(MPISCU::logFile,
			    "SCUTransComplete: finished.\n", numreq );      
#endif
    }

    // Clear the request manager:
    MPISCU::ReqMan->Clear();

}


// Functions used by the SCU syscall emulations.

namespace MPISCU{

// To record whether the following function has been called.
    static bool grid_is_set = false;

/*!
  If MPI is started outwith the MPI-SCU layer then this function should
  be used to tell the MPI-SCU layer what the grid dimensions are.
  This mechanism need and should not be used if the MPI parameters are read
  from a file.
  This function can only be called once; if it is called again it does nothing.
  \param x The grid dimension in the X direction.
  \param y The grid dimension in the Y direction.
  \param z The grid dimension in the Z direction.
  \param t The grid dimension in the T direction.
*/
    void set_pe_grid(int x, int y, int z, int t, int s){

	if(grid_is_set) return; // issue a warning?
    
	peGrid[0] = t;
	peGrid[1] = x;
	peGrid[2] = y;
	peGrid[3] = z;
	peGrid[4] = s;
	
	grid_is_set = true;

    }
    

//----------------------------------------------------------------
/*!
  This function finds the parameters relevant for the parallel
  decomposition of the lattice, sets up the communications layer
  and defines the grid topology.

  It also  defines a root node (which is useful for IO)
  and opens logfiles (if required).

  \ingroup mpicomms
*/
//----------------------------------------------------------------
    void CommsInit(  ) {
    
	int  grid_periodicity[NDIM] = {1,1,1,1,1};  /* Array used to specify periodic BCs */
	int  pe_reorder = 0;      /* Flag to disallow PE reordering for the cart-comm */

	// If we have already been initialized, don't try to do it twice:
	if( Is_Initialised ) return;

   
    
	// Set-up the default values for the comms parameters:
	RNGseed = 1;
	strcpy(seedFileName,"rng.dat");
	strcpy(logFileName, "stderr");

	// If MPI has not already been started somewhere else, start it here.
    
	int initialised_remotely;
	MPI_Initialized(&initialised_remotely);
	if(initialised_remotely){
	    if(!grid_is_set)
		RaiseError("MPISCU_CommsInit: Processor grid must be set using MPISCU::set_pe_grid.");
	}else{
	    int dummy_argc;
	    char  **dummy_argv;
	    MPI_Init(&dummy_argc, &dummy_argv);
	    grid_is_set = true; // just in case.
	}


	// Check the environment variable MPISCU::ENVVAR. If it is defined and
	// is a non-zero length string, then that is either the name of a file
	// from which to read the read the parameters, or it is the parameters
	// themselves. If it is defined but zero length, then we use the 
	// default parameter filename.
	// If it is not defined, then MPISCU::peGrid had better be
	// initialised already.
    
	char *envvar = getenv(ENVVAR);
	if(envvar){ 
	    if(strlen(envvar) == 0) envvar = default_filename;
	    ParseCommsParam(envvar);
	}


	/* Define the (cartesian, periodic) topology of the MPI communicator */

	MPI_Cart_create( MPI_COMM_WORLD,   /* Original communicator */
			 NDIM,             /* No. dimensions */
			 peGrid,           /* No. PEs in each direction */
			 grid_periodicity, /* Periodicity of PE grid in each direction */
			 pe_reorder,       /* True/false flag to allow PE re-ordering */
			 &Cart_Comm      /* The new, cartesian, communicator */
	    );

	/* Look up process number */
	MPI_Comm_rank( Cart_Comm, &peRank );

	/* Look-up processor position */
	MPI_Cart_coords( Cart_Comm, peRank, NDIM, pePos );

	/* Look up number of processors */
	MPI_Comm_size( Cart_Comm, &peNum );
#define MPISCU_DEBUG    
#ifdef MPISCU_DEBUG
	/* Initialise the log-file, which may actually be stdout or stderr */

	if( strcmp(logFileName,"stderr") == 0 )
	    logFile = stderr;
	else if( strcmp(logFileName,"stdout") == 0 )
	    logFile = stdout;
	else {
	    /* Create a logfile name with the PE number as a suffix,
	       such that we have logfile.01, logfile.02 ... logfile.15 etc */

	    /* ordPEnum is no. of orders of magnitude (base10) of number of PEs */ 
	    int ordPEnum = 1 + (int)log10(peNum);
	    char *PEnum = (char*)malloc( ( ordPEnum + 1 ) * sizeof(char) );
	    strcat(logFileName,".");

	    /* Create filenumber based on PE number + leading zeros */
	    sprintf(PEnum,"%i", ((int)(exp(((double)ordPEnum)*log(10.0))))+peRank);
	    strcat(logFileName,&PEnum[1]); /* Skip the leading 1 */
	    free( PEnum );
	    /* Open the logfile associated with this PE */
	    logFile = Fopen(logFileName,"w");
	}
#endif

	/* Inform user that initialization has started */
	printf_all("MPISCU::CommsInit:  Initializing.\n");

	/* identify the root processor as that which lies at pos[i]=0 forall i */
	/* calculate the identifier on every process */
	int root_check = 0;
	for(int idirn = 0; idirn < NDIM; idirn++ ) root_check+=pePos[idirn];
	/* Gather the values of root_check from every PE onto every PE */
	int *root_array = (int*)malloc(peNum*sizeof(int));
	if( root_array == NULL ) 
	    RaiseError("MPISCU::CommsInit: malloc failed for root_array.");
	MPI_Allgather( &root_check, /* Pointer to number to be gathered */
		       1,           /* i.e. gathering a single item */
		       MPI_INT,     /* which is a standard C integer */
		       root_array,  /* Pointer to the array which will recieve the data */
		       1,           /* One thing from each PE */
		       MPI_INT,     /* and that thing is an int. */
		       Cart_Comm /* Using the cartesian communicator */
	    );

	/* Every PE goes through the list and identifies the root PE */
	for(int ir = 0; ir < peNum; ir++ ) 
	    if( root_array[ir] == 0 ) root_pe = ir;



	/* Free the memory associated with the root-checking array */
	free(root_array);

	/* Log that the initialization has completed and give this PEs rank */
#ifdef MPISCU_DEBUG
	fprintf_all(logFile,"MPISCU::CommsInit:  Initialization complete [PE=%i of %i, ROOT_PE=%i].\n",peRank, peNum, root_pe );
#endif


	/* Initialise the MPI Request handler */
	ReqMan = new MPIRequestManager();

	/* Initialise the table of NNs, indexed by SCUDir */
	int dir_index, dummy;
	for( int idim = 0; idim < NDIM; idim++ ) {
	    for( int idir = -1; idir <=1 ; idir+=2 ) {
  		if( idim == 0 && idir == +1 ) dir_index = SCU_TP;
  		if( idim == 0 && idir == -1 ) dir_index = SCU_TM;
		if( idim == 1 && idir == +1 ) dir_index = SCU_XP;
		if( idim == 1 && idir == -1 ) dir_index = SCU_XM;
		if( idim == 2 && idir == +1 ) dir_index = SCU_YP;
		if( idim == 2 && idir == -1 ) dir_index = SCU_YM;
		if( idim == 3 && idir == +1 ) dir_index = SCU_ZP;
		if( idim == 3 && idir == -1 ) dir_index = SCU_ZM;
		if( idim == 4 && idir == +1 ) dir_index = SCU_SP;
		if( idim == 4 && idir == -1 ) dir_index = SCU_SM;
		MPI_Cart_shift( Cart_Comm,  // Using the cartesian communicator
				idim,       // Do this dimension 
				idir,       // Look up nearest neighbour 
				&dummy,     // Rank of this PE
				&(nnList[dir_index]) // Rank of neighbour PE 
		    );
	    }
	}

	Is_Initialised = true;
    
    }


/*!   \ingroup mpicomms */

    void SCUCommsFinalize( void ) {
	if( Is_Initialised ) MPI_Finalize();
    }


/* The global summation */
/*! \ingroup mpicomms collectivecomms */

    void SCUGlobalSum(Type_tag t,   /*!< In: Type of data being summed */
		      size_t tsize, /*!< In: Size of the data type */
		      int n,        /*!< In: Number of values to sum */
		      void *ivec,   /*!< In: Vector of input values */
		      void *ovec    /*!< Out: Vector of output values */
	) {
	MPI_Datatype mpitype; /* This will hold the MPI_Datatype for type (t + size) */

	if( !Is_Initialised ) CommsInit(); 

#ifdef MPISCU_DEBUG
	MPISCU::fprintf(logFile,"SCUGlobalSum: Performing a global summation.\n");
#endif

	/* Check args make sense */
	if( n <= 0 )
	    RaiseError("SCUGlobalSum: no. of values to sum is <= 0!");
	if( ivec == NULL )
	    RaiseError("SCUGlobalSum: input vector points to NULL!");
	if( ovec == NULL )
	    RaiseError("SCUGlobalSum: output vector points to NULL!");

	/* Map the requested type onto an MPI_Datatype */
	mpitype = MPITypeConv( t, tsize );

	/* Invoke the relevent MPI call, so that all processors get the global sum*/
    
	MPI_Allreduce(ivec,         /* Array containing data to be summed */
		      ovec,         /* Array to receive the summations */
		      n,            /* Number of items in the array */
		      mpitype,      /* MPI datatype corresponding to Type_tag */
		      MPI_SUM,      /* Do a global sum operation */
		      Cart_Comm  /* Use the cartesian communicator */
	    );

    }

/* SCU-layer error handler:
   Should map onto the ERR class for the QCDSP code. */
/*!
  Prints an error message to \c stdout and causes the program to exit
  immediately with the value \a EXIT_FAILURE.
  \param errstr The messsage.

  \ingroup mpicomms  
*/

    void RaiseError( char* errstr ) {

	/* Report the error: */
	::fprintf(stderr, "Error: %s\n", errstr);  

	/* Finish with MPI if it has been initialised: */
	if( Is_Initialised ) MPI_Finalize();

	exit(EXIT_FAILURE);
    }


// Extra error wrapper to deal with string literals. */
/*!
  Prints an error message to \c stdout and causes the program to exit
  immediately with the value \a EXIT_FAILURE.
  \param errstr The messsage.

  \ingroup mpicomms
*/
    void RaiseError( const char* errstring ) { 
	RaiseError( const_cast<char*>(errstring) ); 
    }


/*-------------------------------------------------------------------------*/
/*                   Implementation-specific subroutines:                  */
/*              If this were a class, these would be private.              */
/*-------------------------------------------------------------------------*/

//----------------------------------------------------------------
/*!
  The lowest level MPI comms subroutine, on which all other comms calls
  are based.  
*/
//----------------------------------------------------------------
    void Trans(void* addr, MPI_Datatype mpi_dt, SCUDir dir, SCUXR sendrx){
    
	MPI_Request request;

	// Determine the NN in the given direction:
	int nnPE = nnList[dir];

	// Initiate the send or recieve:
	if( sendrx == SCU_SEND ) 
	    MPI_Issend( addr,            // base-address of the data 
			1,               // Number of items to send, one datatype 
			mpi_dt,          // MPI datatype to send 
			nnPE,            // ID of destination PE 
			dir,             // Message-tag based on dirn 
			Cart_Comm,     // The communicator 
			&request         // RETURNS, the request handle 
		);
	else 
	    MPI_Irecv( addr,              // base-address of the data 
		       1,                 // Number of items to recieve, one struct
		       mpi_dt,            // MPI datatype to recv 
		       nnPE,              // ID of source PE 
		       dir-((dir%2)*2-1), // Tag based on dirn 
		       Cart_Comm,       // The communicator 
		       &request           // RETURNS, the request handle 
		);
    
	// Add the new request to the req. handler:
	ReqMan->AddRequest(request);

	return;
    }

//----------------------------------------------------------------
/*!
  Looks up and parses the run-time user parameters specified via ENVVAR.
  i.e. Lots of messy string handling et cetera.
 
  \todo Currently, all PEs open the file and look up the required 
  information.  It would perhaps be quicker to get one PE to look 
  in the file and then distribute the information.
  In fact, if only one node is capable of I/O, this would be 
  neccessary, so it should be done.
*/                                                                     
/* ----------------------------------------------------------------- */
    void ParseCommsParam(char *envvar) {
    
	enum { NULL_READ, GRID_READ, LOGF_READ, SEED_READ, SEEDFILE_READ};

	char *parameters;
	char  f_line[STRING_MAX_LEN], *def_token, *tok_pos;
	int  idirn, read_state, io_state;

	/* NULL the file pointer in case this routine fails */
	logFile = NULL;

	// Determine if the string holds the parameters or a file name
    
	bool read_from_file = true;
	for(int i=0; i<strlen(envvar) && read_from_file==true; i++ )
	    if( envvar[i] == '{' ) read_from_file = false;
	
	if(read_from_file) {	

	    FILE *fp = Fopen(envvar,"r");
	    if( fp == NULL ) RaiseError("MPIParseCommsParam: Could not open comms parameter file!");
	    /* Lookup the file size and define a suitably sized buffer */
	    fseek(fp, 0, SEEK_END);

	    parameters = (char*)malloc( ftell(fp) * sizeof(char) );
	    fseek(fp, 0, SEEK_SET);

	    // Read the file line-by-line and put the whole thing into parameters.
	
	    if(!parameters)
		RaiseError("MPIParseCommsParam: malloc failed for file buffer.");
	    strcpy(parameters,"");

	    while( fscanf(fp,"%[^\n]\n",f_line) != EOF ) strcat(parameters,f_line);
	    Fclose(fp);

	}else
	    parameters = envvar;

	/* Set initial (null/invalid) values for the user parameters: */
	for(int i=0; i<NDIM; i++) peGrid[i] = -1;
	strcpy(logFileName,"!!empty!!");

	/* Now attempt to decipher the definition string held in comm_def */
	/* This is done by breaking the string down into a stream of tokens */
	read_state = NULL_READ;
	tok_pos = parameters;
	while( def_token = (char*)CommsStringTokenizer(parameters, 
							      "{} =,\n\t", 
							      &tok_pos) ) {
		
	    /* Look up the number of processors in each direction. */
	    /* If we find the `grid' token, change into GRID_READ mode */
	    if( strcmp(def_token,"GRID") == 0 ) {
		read_state = GRID_READ;
		idirn = 0;
	    } else if( read_state == GRID_READ ) { 
		/* After `grid', read NDIM*ints into the peGrid array */
		if( idirn < NDIM ) {
		    io_state = sscanf(def_token,"%i",&peGrid[idirn]);
		    idirn++;
		    if( idirn == NDIM ) read_state = NULL_READ; /* Have we finished? */
		}
	    }

	    /* Get the name of the logfile for verbose output */
	    if( strcmp(def_token,"LOGFILE") == 0 ) {
		read_state = LOGF_READ;
	    } else if( read_state == LOGF_READ ) {
		/* Grab the filename token */
		strcpy(logFileName,def_token);
		read_state = NULL_READ; /* and finish */
	    }

	    /* Get the specified RNG seed */
	    if( strcmp(def_token,"SEED") == 0 ) {
		read_state = SEED_READ;
	    } else if( read_state == SEED_READ ) {
		io_state = sscanf(def_token,"%i", &RNGseed );
		read_state = NULL_READ;
	    }

	    /* Get the specified RNG seeds filename*/
	    if( strcmp(def_token,"SEEDFILE") == 0 ) {
		read_state = SEEDFILE_READ;
	    } else if( read_state == SEEDFILE_READ ) {
		io_state = sscanf(def_token,"%s", &seedFileName );
		read_state = NULL_READ;
	    }

	}

	/* Free the memory associated with the file buffer if required */
	if(read_from_file) free(parameters);

	/* If any necessary parameters have not been specified properly, exit */
	/* Checking the processor-element grid specification: */
	for( idirn = 0; idirn < NDIM; idirn++ ) 
	    if( peGrid[idirn] < 0 ) 
		RaiseError("MPIParseCommsParam: Processor array dimensions have not been specified correctly.");
	

    }

//----------------------------------------------------------------
/*!
  String tokenizer, coded here to ensure portability:
*/
//----------------------------------------------------------------

    char *CommsStringTokenizer(char* str, const char* delim, char** tok_pos ) {
	char *tokenstr, *substr;
	int i, tokenstate, toki, isgap, idel, tok_find;
    
	substr = *tok_pos;

	if( substr[0] != '0' ) {
	    // Not at the end of the string, so find the next token:
	    tokenstr = (char*)malloc( strlen(substr) );
	    tokenstate = 0; toki = 0;
	    for( i=0; i<=strlen(substr); i++ ) {
		// Determine if the current character is one of the delimiters;
		idel = 0; tok_find = 0;
		while( idel < strlen(delim)+1 ) { //<The end-of-string 0 is a delimiter.
		    if( substr[i] == delim[idel] ) {
			tok_find = idel+1;
			idel = strlen(delim)+1;
		    }
		    idel++;
		}
		if( tok_find == 0 ) {
		    isgap = 0;
		} else {
		    isgap = 1;
		}

		if( i == 0 && isgap == 0 ) tokenstate = 1;
		if( tokenstate == 0 && isgap == 0 ) {
		    tokenstate = 1;
		    // A token has begun:
		    tokenstr[toki] = substr[i];
		    toki++;
		} else if( tokenstate == 1 ) {
		    if( isgap == 0 ) {
			// A token continues:
			tokenstr[toki] = substr[i];
			toki++;
		    } else {
			// We have found a token, so return it:
			*tok_pos = &(substr[i]);
			tokenstr[toki] = 0;
			return( tokenstr );
		    }
		}
		//printf("tokenizer: %i %i %i\n", tokenstate, isgap, toki );
	    }
	}

	// If the code gets to here we are at the end of the string:
	return(NULL);
    }

//----------------------------------------------------------------
/*!
  On-the-fly type+size -> MPI_Datatype conversion.  There are
  probably somewhat quicker ways of doing this via look-up tables, but
  for now, the look-up has been left explicit.
*/
//----------------------------------------------------------------

    MPI_Datatype MPITypeConv( Type_tag t, size_t tsize ) {

	/* Check arguments make sense */
	if( t != TYPE_IFloat && t != TYPE_int )
	    RaiseError("SCUMPITypeConv: unknown data type!");
	if( tsize <= 0 )
	    RaiseError("SCUMPITypeConv: size of data type is <= 0!");


	MPI_Datatype int_types[N_INT_TYPES] = { MPI_CHAR, MPI_SHORT, MPI_INT, MPI_LONG };
	MPI_Datatype IFloat_types[N_FLOAT_TYPES] = { MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE };
	char err_str[STRING_MAX_LEN];

	int i, mpisize;

	/* Go through int OR IFloat types */

	if( t == TYPE_int ) {
	    for( i=0; i < N_INT_TYPES; i++ ) {
		MPI_Type_size( int_types[i], &mpisize );  /* Get size of this MPI type */
		if( tsize == mpisize )               // If we have a match 
		    return( int_types[i] );          // return the matched type. 
	    
	    }
	    /* if this executes, no suitable type has not been found, so raise an error */
	    sprintf(err_str,"SCUMPITypeConv: no suitable %i-byte int type among MPI primitive types",tsize);
	    RaiseError(err_str);

	    /* IFloat types */
	} else if( t == TYPE_IFloat ) {
	    for( i=0; i < N_FLOAT_TYPES; i++ ) {
		MPI_Type_size( IFloat_types[i], &mpisize );  /* Get size of this MPI type */
		if( tsize == mpisize )                // If we have a match
		    return( IFloat_types[i] );        // return the matched type.
	    
	    }
	    /* if this executes, no suitable type has not been found, so raise an error */
	    sprintf(err_str,"no suitable %i-byte IFloat type among MPI primitive types",tsize);
	    RaiseError(err_str);

	}

	/* This statement should never execute, however just to check, and keep 
	   the compiler happy, we shall say: */
	RaiseError("SCUMPITypeConv: running unrunnable section of SCUMPITypeConv!  Possible memory access problem?");
	return(( MPI_Datatype)0 );

    }

//----------------------------------------------------------------
/*!
  Reads a seed for every PE from a file specified during intialisation.
*/
//----------------------------------------------------------------

    unsigned int ReadSeedFile( void ) {
	FILE *seedfp = NULL;
	int i, io_state, n = peNum;
	unsigned int *iseed, seed; 

	// Check we have actually been initialised:
	if( !Is_Initialised ) CommsInit();

#ifdef MPISCU_DEBUG
	fprintf_all(logFile,"MPISCU::ReadSeedFile: Opening seed file %s.\n", seedFileName);
#endif

	// Create the seeds buffer:
	iseed = new unsigned int[n+1];
    
	// Open the file:
	seedfp = Fopen(seedFileName, "r" );
	if( seedfp == NULL ) 
	    RaiseError("SCUReadSeedFile: could not open seed file!\n");
    
	// Read in n seeds:
	i = 0; while( i < n && fscanf(seedfp,"%u",&(iseed[i])) != EOF ) i++;
    
	// Close the file:
	Fclose(seedfp); seedfp = NULL;
    
	// Die if the file ended before all seeds had been read in:
	if( i < n )
	    RaiseError("SCUReadSeedFile: not enough seeds have been supplied in the seed file!");
	// XXX, EXTREME warning.  This killed one thread and then hung.
    
	// Get the seed which belongs to this PE:
	seed = iseed[peRank];

	// Delete the seeds buffer:
	delete [] iseed;

	// Return this PE's seed:
	return seed;
    }
    
} //namespace MPISCU






CPS_END_NAMESPACE
