/*
 * PGF routine
 */

#include "param.h"

void pgf(float u3[][NYDIM][NZDIM],float v3[][NYDIM][NZDIM],float w3[][NYDIM][NZDIM],float p1[][NYDIM][NZDIM],float p3[][NYDIM][NZDIM],float t1[][NYDIM][NZDIM],float tstep)
{
    void bc(float t[][NYDIM][NZDIM],float p[][NYDIM][NZDIM],float u[][NYDIM][NZDIM],float v[][NYDIM][NZDIM],float w[][NYDIM][NZDIM]);
    
    int i,j,k;
    float g = 9.81;
    /* u */
    for(i = I1+1; i<=I2; i++)
        for(j = J1; j<=J2; j++)
            for(k = K1; k<=K2; k++)
                u3[i][j][k] -= tstep*(p1[i][j][k] - p1[i-1][j][k])/dx/du[k-K1];
    
    /* v */
    for(i = I1; i<=I2; i++)
        for(j = J1+1; j<=J2; j++)
            for(k = K1; k<=K2; k++)
                v3[i][j][k] -= tstep*(p1[i][j][k] - p1[i][j-1][k])/dy/du[k-K1];
    
    /* w */
    for(i = I1; i<=I2; i++){
        for(j = J1; j<=J2; j++){
            for(k = K1+1; k<=K2; k++){
                w3[i][j][k] -= tstep*(p1[i][j][k] - p1[i][j][k-1])/dz/dw[k-K1];
                w3[i][j][k] += tstep*g*(t1[i][j][k]+t1[i][j][k-1]-2.0*TBAR)/2.0/TBAR;
            }
        }
    }
    
    
    /* Apply bc to u3 and w3 */
    bc(t1,p1,u3,v3,w3);
    
    /* p' */
    for(i = I1; i<=I2; i++){
        for(j = J1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                float a1 = du[k-K1]*(u3[i+1][j][k] - u3[i][j][k]);
                float a2 = du[k-K1]*(v3[i][j+1][k] - v3[i][j][k]);
                float a3 = dw[k+1-K1]*w3[i][j][k+1] - dw[k-K1]*w3[i][j][k];
                p3[i][j][k] -= tstep*CS*CS*(a1/dx+a2/dy+a3/dz);
            }
        }
    }
    
}
