/*
  $Id: main.C,v 1.4 2008/02/08 18:35:08 chulwoo Exp $
*/

#include<config.h>
#include <util/qcdio.h>
#include <math.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/dirac_op.h>
#include<util/error.h>
#include<util/time_cps.h>
#include<util/lat_data.h>
#include<comms/scu.h>
#include<alg/alg_hmd.h>
#include<alg/do_arg.h>
#if TARGET == QCDOC
#include <qalloc.h>
extern "C"{
void _mcleanup(void);
}
#endif




USING_NAMESPACE_CPS



const char *f_wilson_test_filename = CWDPREFIX("f_wilson_test");
const char *psi_filename = CWDPREFIX("psi");


int main(int argc,char *argv[]){

#if TARGET == QCDOC
    DefaultSetup();
    printf("Sizes = %d %d %d %d %d %d\n",SizeX(),SizeY(),SizeZ(),SizeT(),SizeS(),SizeW());
    printf("Coors = %d %d %d %d %d %d\n",CoorX(),CoorY(),CoorZ(),CoorT(),CoorS(),CoorW());
#endif
    FILE *fp;
    double dtime;

    //----------------------------------------------------------------
    // Initializes all Global Job Parameters
    //----------------------------------------------------------------
    DoArg do_arg;
    int nx,ny,nz,nt;

    if (argc < 5){
        ERR.General("f_wilson_test","main()","usage: %s nx ny nz nt\n",argv[0]);
    }
    sscanf(argv[1],"%d",&nx);
    sscanf(argv[2],"%d",&ny);
    sscanf(argv[3],"%d",&nz);
    sscanf(argv[4],"%d",&nt);
    printf("total sites = %d %d %d %d\n",nx,ny,nz,nt);
#if TARGET == QCDOC
    do_arg.x_node_sites = nx/SizeX();
    do_arg.y_node_sites = ny/SizeY();
    do_arg.z_node_sites = nz/SizeZ();
    do_arg.t_node_sites = nt/SizeT();
    do_arg.s_node_sites = 0;
    do_arg.x_nodes = SizeX();
    do_arg.y_nodes = SizeY();
    do_arg.z_nodes = SizeZ();
    do_arg.t_nodes = SizeT();
    do_arg.s_nodes = 1;
#else
    do_arg.x_node_sites = nx;
    do_arg.y_node_sites = ny;
    do_arg.z_node_sites = nz;
    do_arg.t_node_sites = nt;
    do_arg.s_node_sites = 0;
    do_arg.x_nodes = 1;
    do_arg.y_nodes = 1;
    do_arg.z_nodes = 1;
    do_arg.t_nodes = 1;
    do_arg.s_nodes = 1;
#endif 
    do_arg.x_bc = BND_CND_PRD;
    do_arg.y_bc = BND_CND_PRD;
    do_arg.z_bc = BND_CND_PRD;
    do_arg.t_bc = BND_CND_APRD;

    do_arg.start_conf_kind = START_CONF_DISORD;

    do_arg.start_seed_kind = START_SEED_FIXED;
    do_arg.beta = 5.5;
    do_arg.dwf_height = 0.9;
    do_arg.clover_coeff = 2.0171;
 //   do_arg.verbose_level = -1205;

    do_arg.asqtad_KS = (1.0/8.0)+(6.0/16.0)+(1.0/8.0);
    do_arg.asqtad_naik = -1.0/24.0;
    do_arg.asqtad_3staple = (-1.0/8.0)*0.5;
    do_arg.asqtad_5staple = ( 1.0/8.0)*0.25*0.5;
    do_arg.asqtad_7staple = (-1.0/8.0)*0.125*(1.0/6.0);
    do_arg.asqtad_lepage = -1.0/16;

	VRB.Level(0);
	VRB.ActivateLevel(VERBOSE_RNGSEED_LEVEL);
	VRB.ActivateLevel(VERBOSE_FUNC_LEVEL);
	VRB.ActivateLevel(VERBOSE_FLOW_LEVEL);
    
    CgArg cg_arg;

    cg_arg.mass = 0.1;
    cg_arg.stop_rsd = 1e-12;
    cg_arg.max_num_iter = 500;
    GJP.Initialize(do_arg);

//    VRB.Level(GJP.VerboseLevel());

    fp = Fopen(ADD_ID,"f_wilson_test.out","w");

    GwilsonFwilson lat;

#if 0
    Vector *result = 
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    Vector *X_out =
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    Vector *X_out2 =
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));

    if(!result) ERR.Pointer("","","result");
    if(!X_out) ERR.Pointer("","","X_out");
    if(!X_out2) ERR.Pointer("","","X_out2");
    Vector *X_in =
	(Vector*)smalloc(GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
	bzero(X_in,GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
    if(!X_in) ERR.Pointer("","","X_in");
#else
    LatVector result(lat.SpinComponents());
    LatVector X_out(lat.SpinComponents());
    LatVector X_out2(lat.SpinComponents());
    LatVector X_in(lat.SpinComponents());
#endif


    int s[4];
#if 1
	lat.RandGaussVector(X_in,1.0);
#else

    Matrix *gf = lat.GaugeField();
    IFloat *gf_p = (IFloat *)lat.GaugeField();

    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++)
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {

		    int n = lat.FsiteOffset(s);

#if TARGET==QCDOC
		  if(GJP.XnodeCoor()==0 && GJP.YnodeCoor()==0 && GJP.ZnodeCoor()==0 && GJP.TnodeCoor()==0 &&n==0) crd=1.0; else crd = 0.0;
#else
	if(n==0) crd = 1.0; else crd = 0.0;
#endif
					
		    for(int v=0; v<6; v+=2){ 
			if (v==0)
			*((IFloat*)&X_in[n]+v) = crd;
			else
			*((IFloat*)&X_in[n]+v) = 0;
			*((IFloat*)&X_in[n]+v+1) = 0.0;
		    }
		}
#endif

    Vector *out;
    DiracOpWilson dirac(lat,X_out,X_in,&cg_arg,CNV_FRM_NO);

	for(int k = 0; k< 1; k++){
    	double maxdiff=0.;
		printf("k=%d ",k);
		if (k ==0)
			out = result;
		else
			out = X_out;
		bzero((char *)out, GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
		lat.Fconvert(out,WILSON,CANONICAL);
		lat.Fconvert(X_in,WILSON,CANONICAL);
		int offset = GJP.VolNodeSites()*lat.FsiteSize()/ (2*6);
#if 1
#if 1
#if TARGET==QCDOC
		int vol = nx*ny*nz*nt/(SizeX()*SizeY()*SizeZ()*SizeT());
#else
		int vol = nx*ny*nz*nt;
#endif
		dtime = -dclock();
   		int iter = dirac.MatInv(out,X_in);
		dtime +=dclock();
#else
		dirac.Dslash(out,X_in+offset,CHKB_EVEN,DAG_NO);
		dirac.Dslash(out+offset,X_in,CHKB_ODD,DAG_NO);
#endif
#endif

		if (k == 0){
			bzero((char *)X_out2, GJP.VolNodeSites()*lat.FsiteSize()*sizeof(IFloat));
			dirac.Dslash(X_out2,out+offset,CHKB_EVEN,DAG_NO);
			dirac.Dslash(X_out2+offset,out,CHKB_ODD,DAG_NO);
			lat.Fconvert(X_out2,CANONICAL,WILSON);
		}
		lat.Fconvert(out,CANONICAL,WILSON);
		lat.Fconvert(X_in,CANONICAL,WILSON);
		X_out2->FTimesV1PlusV2(-0.5/(cg_arg.mass+4.0),X_out2,out,GJP.VolNodeSites()*lat.FsiteSize());
    
    Float dummy;
    Float dt = 2;

    for(s[3]=0; s[3]<GJP.NodeSites(3); s[3]++) 
	for(s[2]=0; s[2]<GJP.NodeSites(2); s[2]++)
	    for(s[1]=0; s[1]<GJP.NodeSites(1); s[1]++)
		for(s[0]=0; s[0]<GJP.NodeSites(0); s[0]++) {

		    int n = lat.FsiteOffset(s)*lat.SpinComponents();
			for(int i=0; i<(lat.FsiteSize()/2); i++){
#if TARGET == QCDOC
		    if ( k==0 )
				Fprintf(ADD_ID, fp," %d %d %d %d %d ", CoorX()*GJP.NodeSites(0)+s[0], CoorY()*GJP.NodeSites(1)+s[1], CoorZ()*GJP.NodeSites(2)+s[2], CoorT()*GJP.NodeSites(3)+s[3], i);
#else
		    if ( k==0 )
				Fprintf(ADD_ID, fp," %d %d %d %d %d ", s[0], s[1], s[2], s[3], i);
#endif
		    if ( k==0 )
				Fprintf(ADD_ID, fp," (%0.7e %0.7e) (%0.7e %0.7e)",
				*((IFloat*)&result[n]+i*2), *((IFloat*)&result[n]+i*2+1),
				*((IFloat*)&X_in[n]+i*2), *((IFloat*)&X_in[n]+i*2+1));
#if 0
				Fprintf(ADD_ID, fp," (%0.2e %0.2e)\n",
#if 0
		*((IFloat*)&X_out2[n]+i*2)-*((IFloat*)&X_in[n]+i*2), 
	*((IFloat*)&X_out2[n]+i*2+1)-*((IFloat*)&X_in[n]+i* 2+1));
#else
		*((IFloat*)&X_out2[n]+i*2),
	*((IFloat*)&X_out2[n]+i*2+1));
#endif
#else
				Fprintf(ADD_ID, fp,"\n");
#endif
	double diff =	*((IFloat*)&X_out2[n]+i*2)-*((IFloat*)&X_in[n]+i*2);
        if (fabs(diff)>maxdiff) maxdiff = fabs(diff);
 	diff = *((IFloat*)&X_out2[n]+i*2+1)-*((IFloat*)&X_in[n]+i* 2+1);
        if (fabs(diff)>maxdiff) maxdiff = fabs(diff);
			}
		}
    printf("Max diff between X_in and M*X_out = %0.2e\n", maxdiff);
}
    Fclose(fp);
    
    sfree(X_in);
    sfree(result);
    sfree(X_out);
    sfree(X_out2);
    return 0; 
}
