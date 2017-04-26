/* 
 * This routine calculate the diffusion terms
 */

#include "param.h"

void diffusion(float u1[][NYDIM][NZDIM],float u3[][NYDIM][NZDIM],float v1[][NYDIM][NZDIM],float v3[][NYDIM][NZDIM],float w1[][NYDIM][NZDIM],float w3[][NYDIM][NZDIM],float t1[][NYDIM][NZDIM],float t2[][NYDIM][NZDIM],float tstep)
{
    int i,j,k;
    /* u */
    for(i = I1+1; i<=I2; i++){
        for(j = J1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                float a1=u1[i+1][j][k]-2.0*u1[i][j][k]+u1[i-1][j][k];
                float a2=u1[i][j+1][k]-2.0*u1[i][j][k]+u1[i][j-1][k];
                float a3=u1[i][j][k+1]-2.0*u1[i][j][k]+u1[i][j][k-1];
                u3[i][j][k]+=tstep*km*(a1/dx/dx+a2/dy/dy+a3/dz/dz);
            }
        }
    }
    
    /* v */
    for(i = I1; i<=I2; i++){
        for(j = J1+1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                float a1=v1[i+1][j][k]-2.0*v1[i][j][k]+v1[i-1][j][k];
                float a2=v1[i][j+1][k]-2.0*v1[i][j][k]+v1[i][j-1][k];
                float a3=v1[i][j][k+1]-2.0*v1[i][j][k]+v1[i][j][k-1];
                v3[i][j][k]+=tstep*km*(a1/dx/dx+a2/dy/dy+a3/dz/dz);
            }
        }
    }
    
    /* w */
    for(i = I1; i<=I2; i++){
        for(j = J1+1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                float a1=w1[i+1][j][k]-2.0*w1[i][j][k]+w1[i-1][j][k];
                float a2=w1[i][j+1][k]-2.0*w1[i][j][k]+w1[i][j-1][k];
                float a3=w1[i][j][k+1]-2.0*w1[i][j][k]+w1[i][j][k-1];
                w3[i][j][k]+=tstep*km*(a1/dx/dx+a2/dy/dy+a3/dz/dz);
            }
        }
    }
    
    /* theta */
    for(i = I1; i<=I2; i++){
        for(j = J1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                float a1=t1[i+1][j][k]-2.0*t1[i][j][k]+t1[i-1][j][k];
                float a2=t1[i][j+1][k]-2.0*t1[i][j][k]+t1[i][j-1][k];
                float a3=t1[i][j][k+1]-2.0*t1[i][j][k]+t1[i][j][k-1];
            t2[i][j][k] += dt*kmt*(a1/dx/dx+a2/dy/dy+a3/dz/dz);
            }
        }
    }
    
}
