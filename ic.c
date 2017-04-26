/*
 * ============================ ic =====================
 * IC sets the initial condition
 * Arguments:
 *
 *	s1	real array	IC data. Set I1..I2 here;
 *				  [I1-1],[I2+1] = ghost zones
 *				  if 1 ghost point on each side
 *	dx,dy	real		grid spacing
 *	I1,I2,J1,J2	integers	indices bounding array data
 *	u,v	real	array	Flow velocity grids
 *	nx,ny	real	number of grid points
 *	nydim	real	number of grid points including ghost
 */
#include <math.h>
#include <stdio.h>
#include "param.h"

void ic(float t1[][NYDIM][NZDIM],float p1[][NYDIM][NZDIM],float p2[][NYDIM][NZDIM],float u1[][NYDIM][NZDIM],float u2[][NYDIM][NZDIM],float v1[][NYDIM][NZDIM],float v2[][NYDIM][NZDIM],float w1[][NYDIM][NZDIM],float w2[][NYDIM][NZDIM])
{
    int i,j,k,m;
    /* Initial velocity field */
    for(i = I1; i<=I2+1; i++){
        for(j = J1; j<=J2+1; j++){
            for(k = K1; k<=K2+1; k++){
                u1[i][j][k] = 0.0;
                u2[i][j][k] = 0.0;
                v1[i][j][k] = 0.0;
                v2[i][j][k] = 0.0;
                w1[i][j][k] = 0.0;
                w2[i][j][k] = 0.0;
            }
        }
    }
    /* Perturbation Pressure */
    for(i = I1; i<=I2; i++){
        for(j = J1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                p1[i][j][k] = 0.0;
                p2[i][j][k] = 0.0;
            }
        }
    }
        
    /* base-state density */
    float g = 9.81,cp = 1004.0,rd = 287.0,p0 = 100000.0;
    for(k = 0; k<NZ; k++){
        float z = dz/2.0+dz*(float)k;
        float T = 300.0 - g/cp*z;
        float pbar = p0*pow(T/TBAR,cp/rd);
        
        du[k] = pbar/rd/T;
    }
    for(k = 0; k<NZ+1; k++){
        float z = dz*(float)k;
        float T = 300.0 - g/cp*z;
        float pbar = p0*pow(T/TBAR,cp/rd);
        
        dw[k] = pbar/rd/T;
    }
    
    /* Temperature */
    for(i = I1; i<=I2; i++){
        for(j = J1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                t1[i][j][k] = TBAR;
                float x = dx/2.0 + dx*(float)(i-I1);
                float y = dy/2.0 + dy*(float)(j-J1);
                float z = dz/2.0 + dz*(float)(k-K1);
                for(m = 0; m<2; m++){
                    float rm = sqrt(pow((x-xc[m])/xradius[m],2.0)+pow((y-yc[m])/yradius[m],2.0)+pow((z-zc[m])/zradius[m],2.0));
                    if(rm<=1.0){
                        t1[i][j][k] += dtp[m]*(cos(rm*M_PI)+1.0)/2.0;
                    }
                }
            }
        }
    }
    
    
    

	return;
}

