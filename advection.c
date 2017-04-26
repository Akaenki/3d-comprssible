/*
 * This routine calculate the advection terms
 * Box and spliting
 */
#include "param.h"

void advection(float t2[][NYDIM][NZDIM],float u2[][NYDIM][NZDIM],float u3[][NYDIM][NZDIM],float v2[][NYDIM][NZDIM],float v3[][NYDIM][NZDIM],float w2[][NYDIM][NZDIM],float w3[][NYDIM][NZDIM],float tstep)
{
	void advect1d(float s1[], float s2[], float c[], float dx, float dt,int i1, int i2, int nx);
    void bc(float t[][NYDIM][NZDIM],float p[][NYDIM][NZDIM],float u[][NYDIM][NZDIM],float v[][NYDIM][NZDIM],float w[][NYDIM][NZDIM]);
    float dummy[NXDIM][NYDIM][NZDIM];
    
    int i,j,k;
    /* u */
    for(i = I1+1; i<=I2; i++){
        for(j = J1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                float a1 = (u2[i+1][j][k]+u2[i][j][k]) * (u2[i+1][j][k]-u2[i][j][k])
                         + (u2[i][j][k]+u2[i-1][j][k]) * (u2[i][j][k]-u2[i-1][j][k]);
                float a2 = (v2[i][j+1][k]+v2[i-1][j+1][k]) * (u2[i][j+1][k]-u2[i][j][k])
                         + (v2[i][j][k]+v2[i-1][j][k]) * (u2[i][j][k]-u2[i][j-1][k]);
                float a3 = (w2[i][j][k+1]+w2[i-1][j][k+1]) * (u2[i][j][k+1]-u2[i][j][k])
                         + (w2[i][j][k]+w2[i-1][j][k]) * (u2[i][j][k]-u2[i][j][k-1]);
                u3[i][j][k] -= tstep*(1.0/(4.0*dx)*a1 + 1.0/(4.0*dy)*a2 + 1.0/(4.0*dz)*a3);
            }
        }
    }
    
    /* v */
    for(i = I1; i<=I2; i++){
        for(j = J1+1; j<=J2; j++){
            for(k = K1; k<=K2; k++){
                float a1 = (u2[i+1][j][k]+u2[i+1][j-1][k]) * (v2[i+1][j][k]-v2[i][j][k])
                         + (u2[i][j][k]+u2[i][j-1][k]) * (v2[i][j][k]-v2[i-1][j][k]);
                float a2 = (v2[i][j+1][k]+v2[i][j][k]) * (v2[i][j+1][k]-v2[i][j][k])
                         + (v2[i][j][k]+v2[i][j-1][k]) * (v2[i][j][k]-v2[i][j-1][k]);
                float a3 = (w2[i][j][k+1]+w2[i][j-1][k+1]) * (v2[i][j][k+1]-v2[i][j][k])
                         + (w2[i][j][k]+w2[i][j-1][k]) * (v2[i][j][k]-v2[i][j][k-1]);
                v3[i][j][k] -= tstep*(1.0/(4.0*dx)*a1 + 1.0/(4.0*dy)*a2 + 1.0/(4.0*dz)*a3);
            }
        }
    }

    /* w */
    for(i = I1; i<=I2; i++){
        for(j = J1; j<=J2; j++){
            for(k = K1+1; k<=K2; k++){
                float a1 = (u2[i+1][j][k]+u2[i+1][j][k-1]) * (w2[i+1][j][k]-w2[i][j][k])
                         + (u2[i][j][k]+u2[i][j][k-1]) * (w2[i][j][k]-w2[i-1][j][k]);
                float a2 = (v2[i][j+1][k]+v2[i][j+1][k-1]) * (w2[i][j+1][k]-w2[i][j][k])
                         + (v2[i][j][k]+v2[i][j][k-1]) * (w2[i][j][k]-w2[i][j-1][k]);
                float a3 = (w2[i][j][k+1]+w2[i][j][k]) * (w2[i][j][k+1]-w2[i][j][k])
                         + (w2[i][j][k]+w2[i][j][k-1]) * (w2[i][j][k]-w2[i][j][k-1]);
                v3[i][j][k] -= tstep*(1.0/(4.0*dx)*a1 + 1.0/(4.0*dy)*a2 + 1.0/(4.0*dz)*a3);
            }
        }
    }
    
    float t1x[NXDIM],t1y[NZDIM],t1z[NZDIM],t1x2[NXDIM],t1y2[NZDIM],t1z2[NZDIM],
    u1x[NXDIM],v1y[NYDIM],w1z[NZDIM];
	/* theta: 3D spliting Fx(dt/2) */
	for (k = K1; k <= K2; k++) {
        for (j = J1; j<=J2; j++) {
            for (i = 0; i < NXDIM; i++) {
			t1x[i] = t2[i][j][k];
            }
            for (i = 0; i < NXDIM; i++) {
			u1x[i] = u2[i][j][k];
            }
            advect1d(t1x, t1x2, u1x, dx, dt/2.0, I1, I2, NX);

            /* Update */
            for (i = 0; i < NXDIM; i++) {
			t2[i][j][k] += t1x2[i];
            }
        }
	}
    bc(t2,dummy,dummy,dummy,dummy);
    
    /* theta: 3D spliting Fy(dt/2) */
    for (i = I1; i <= I2; i++) {
        for (k = K1; k<=K2; k++) {
            for (j = 0; j < NYDIM; j++) {
                t1y[j] = t2[i][j][k];
            }
            for (j = 0; j < NYDIM; j++) {
                v1y[k] = v2[i][j][k];
            }
            advect1d(t1y, t1y2, v1y, dy, dt/2.0, J1, J2, NY);
            
            for (j = 0; j < NYDIM; j++) {
                t2[i][j][k] += t1y2[k];
            }
        }
    }
    bc(t2,dummy,dummy,dummy,dummy);
    
	/* theta: 3D spliting Fz(dt) */
	for (i = I1; i <= I2; i++) {
        for (j = J1; j<=J2; j++) {
            for (k = 0; k < NZDIM; k++) {
                t1z[k] = t2[i][j][k];
            }
            for (k = 0; k < NZDIM; k++) {
                w1z[k] = w2[i][j][k];
            }
            advect1d(t1z, t1z2, w1z, dx, dt, K1, K2, NZ);
		
            for (k = 0; k < NZDIM; k++) {
                t2[i][j][k] += t1z2[k];
            }
        }
	}
    bc(t2,dummy,dummy,dummy,dummy);
    
    /* theta: 3D spliting Fy(dt/2) */
    for (i = I1; i <= I2; i++) {
        for (k = K1; k<=K2; k++) {
            for (j = 0; j < NYDIM; j++) {
                t1y[j] = t2[i][j][k];
            }
            for (j = 0; j < NYDIM; j++) {
                v1y[k] = v2[i][j][k];
            }
            advect1d(t1y, t1y2, v1y, dy, dt/2.0, J1, J2, NY);
            
            for (j = 0; j < NYDIM; j++) {
                t2[i][j][k] += t1y2[k];
            }
        }
    }
    bc(t2,dummy,dummy,dummy,dummy);
    
    /* theta: 3D spliting Fx(dt/2) */
    for (k = K1; k <= K2; k++) {
        for (j = J1; j<=J2; j++) {
            for (i = 0; i < NXDIM; i++) {
                t1x[i] = t2[i][j][k];
            }
            for (i = 0; i < NXDIM; i++) {
                u1x[i] = u2[i][j][k];
            }
            advect1d(t1x, t1x2, u1x, dx, dt/2.0, I1, I2, NX);
            
            /* Update */
            for (i = 0; i < NXDIM; i++) {
                t2[i][j][k] += t1x2[i];
            }
        }
    }

	return;
}

