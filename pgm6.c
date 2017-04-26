/*
 *  ATMS 502 / CSE 566 -- Spring, 2017
 *  Demo for pgm1:  Linear and nonlinear advection
 *  =====>>>>> Linling Miao <<<<<=====
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ncarg/ncargC.h>
#include <ncarg/gks.h>
#define IWTYPE 1
#define WKID   1

#include "param.h"

main()
{

/*
 * Definitions
 */
    
char *name  = "Linling Miao";
char title[25];

/* Arrays and other variables */

    float t1[NXDIM][NYDIM][NZDIM],t2[NXDIM][NYDIM][NZDIM],
    p1[NXDIM][NYDIM][NZDIM],p2[NXDIM][NYDIM][NZDIM],p3[NXDIM][NYDIM][NZDIM],
    u1[NXDIM][NYDIM][NZDIM],u2[NXDIM][NYDIM][NZDIM],u3[NXDIM][NYDIM][NZDIM],
    v1[NXDIM][NYDIM][NZDIM],v2[NXDIM][NYDIM][NZDIM],v3[NXDIM][NYDIM][NZDIM],
    w1[NXDIM][NYDIM][NZDIM],w2[NXDIM][NYDIM][NZDIM],w3[NXDIM][NYDIM][NZDIM];
    float umax,umin,wmax,wmin,tmax,tmin,pmax,pmin;
    float tstep;
	int i,j,k,n,nstep,nplot;

/* Variables to reverse default black/white colors in NCAR Graphics */

	/*Gcolr_rep rgb1,rgb2;*/

/* Function prototype declarations */

    void ic(float t1[][NYDIM][NZDIM],float p1[][NYDIM][NZDIM],float p2[][NYDIM][NZDIM],float u1[][NYDIM][NZDIM],float u2[][NYDIM][NZDIM],float v1[][NYDIM][NZDIM],float v2[][NYDIM][NZDIM],float w1[][NYDIM][NZDIM],float w2[][NYDIM][NZDIM]);
    
	void stats(float s[][NYDIM][NZDIM],float *smax,float *smin);
    
    void bc(float t[][NYDIM][NZDIM],float p[][NYDIM][NZDIM],float u[][NYDIM][NZDIM],float v[][NYDIM][NZDIM],float w[][NYDIM][NZDIM]);
    
    void advection(float t2[][NYDIM][NZDIM],float u2[][NYDIM][NZDIM],float u3[][NYDIM][NZDIM],float v2[][NYDIM][NZDIM],float v3[][NYDIM][NZDIM],float w2[][NYDIM][NZDIM],float w3[][NYDIM][NZDIM],float tstep);
    
    void diffusion(float u1[][NYDIM][NZDIM],float u3[][NYDIM][NZDIM],float v1[][NYDIM][NZDIM],float v3[][NYDIM][NZDIM],float w1[][NYDIM][NZDIM],float w3[][NYDIM][NZDIM],float t1[][NYDIM][NZDIM],float t2[][NYDIM][NZDIM],float tstep);
    
    void pgf(float u3[][NYDIM][NZDIM],float v3[][NYDIM][NZDIM],float w3[][NYDIM][NZDIM],float p1[][NYDIM][NZDIM],float p3[][NYDIM][NZDIM],float t1[][NYDIM][NZDIM],float tstep);
    
    void update(float u1[][NYDIM][NZDIM],float u2[][NYDIM][NZDIM],float u3[][NYDIM][NZDIM],float v1[][NYDIM][NZDIM],float v2[][NYDIM][NZDIM],float v3[][NYDIM][NZDIM],float w1[][NYDIM][NZDIM],float w2[][NYDIM][NZDIM],float w3[][NYDIM][NZDIM],float p1[][NYDIM][NZDIM],float p2[][NYDIM][NZDIM],float p3[][NYDIM][NZDIM],float t1[][NYDIM][NZDIM],float t2[][NYDIM][NZDIM]);

/* Parameters and input .................................... */

	printf("Program #6       Numerical Fluid Dynamics\n\n");
	printf("NX=%d, NY=%d NZ=%d BC_WIDTH=%d",NX,NY,NZ,BC_WIDTH);
	
	printf("Reading input file....\n");
    FILE *input = fopen("input.txt","r");
    fscanf(input,"dx = %f\n",&dx);
    dy = dx; dz = dx;
    fscanf(input,"dt = %f\n",&dt);
    fscanf(input,"nstep = %d\n",&nstep);
    fscanf(input,"nplot = %d\n",&nplot);
    fscanf(input,"center = %f %f %f %f %f %f\n",&xc[0],&yc[0],&zc[0],&xc[1],&yc[1],&zc[1]);
    fscanf(input,"radius = %f %f %f %f %f %f\n",&xradius[0],&yradius[0],&zradius[0],&xradius[1],&yradius[1],&zradius[1]);
    fscanf(input,"magnitude = %f %f\n",&dtp[0],&dtp[1]);
    fscanf(input,"km_v = %f\n",&km);
    fscanf(input,"km_t = %f\n",&kmt);
    fclose(input);

 	printf("Time step dt = %6.3f, Grid spacing dx = dz = %6.3f\n",dt,dx);
    printf("Temperature perturbation Info:\n");
    printf("Center: (%.1f,%.1f,%.1f) (%.1f,%.1f,%.1f)\n Radius: (%.1f,%.1f) (%.1f,%.1f)\n Magnitude: %f,%f\n",xc[0],yc[0],zc[0],xc[1],yc[0],zc[1],xradius[0],zradius[0],xradius[1],zradius[1],dtp[0],dtp[1]);
    printf("Running for %d steps ...\n",nstep);
/*
 * Open the NCAR Graphics package and set colors.
 */
	/*gopen_gks("stdout",0);
	gopen_ws(WKID, NULL, IWTYPE);
	gactivate_ws(WKID);*/

	/* omit following four lines to invert black/white colors */
	/*rgb1.rgb.red = 1.; rgb1.rgb.green = 1.; rgb1.rgb.blue = 1.;
	rgb2.rgb.red = 0.; rgb2.rgb.green = 0.; rgb2.rgb.blue = 0.;
    gset_colr_rep(WKID,0,&rgb1);
    gset_colr_rep(WKID,1,&rgb2);*/
/*
 * Set and plot the initial condition
 */
	ic(t1,p1,p2,u1,u2,v1,v2,w1,w2);
    
    printf("density: %.2f %.2f\n",du[NZ/2],dw[NZ/2]);
    
/*  . . . Set boundary conditions				*/
    bc(t1,p1,u1,v1,w1);
    bc(t2,p2,u2,v2,w2);

    
/*  . . . Set up tstep				*/
    tstep = dt;
    
/*
 * .. Integrate .....................................................
 */

	for (n=1; n<=nstep; n++) {
/*  . . . Initialize				*/
        for(i = 0; i<NXDIM; i++){
            for(j = 0; j<NYDIM; j++){
                for(k = 0; k<NZDIM; k++){
                    u3[i][j][k] = u1[i][j][k];
                    w3[i][j][k] = w1[i][j][k];
                    p3[i][j][k] = p1[i][j][k];
                    t2[i][j][k] = t1[i][j][k];
                }
            }
        }
        
/*  . . . Compute values at next step				*/
		advection(t2,u2,u3,v2,v3,w2,w3,tstep);
        diffusion(u1,u3,v1,v3,w1,w3,t1,t2,tstep);
        pgf(u3,v3,w3,p1,p3,t1,tstep);

/*  . . . Do array update at end of time step			*/
		if(n==1) update(u2,u2,u3,v2,v2,v3,w2,w2,w3,p2,p2,p3,t1,t2);
        else update(u1,u2,u3,v1,v2,v3,w1,w2,w3,p1,p2,p3,t1,t2);
        
/*  . . . Update tstep				*/
        if(n==1) tstep = 2.0*dt;
        
/*  . . . Set boundary conditions				*/
        bc(t1,p1,u1,v1,w1);
        bc(t2,p2,u2,v2,w2);

/*  . . . Plot fields when needed				*/
		if (n == nstep || n%nplot == 0) {
			/* Plot Contours */
            
	    }
        
    }	/* end of time loop n = 1,...,nstep */

    /* Error Calculation */
	/*errors(NX, NY, NYDIM, s1, strue, I1, I2, J1, J2);*/

/*
 * Run complete - do final plots
 */

    

/*
 *  Close the graphics package and stop.
 */

	/*gdeactivate_ws(WKID);
	gclose_ws(WKID);
	gclose_gks();*/

	exit;
}
