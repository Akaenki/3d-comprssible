/*
 * This routine calculate the 1D advection 
 * Piece-wise linear, non-monotonic
 */
#include <math.h>

void advect1d(float s1[],float s2[],float c[],float dx,float dt,
                   int i1,int i2,int nx)
{
	int i;
	float ds[nx+2*i1+1],f[nx+2*i1+1];
    
    /* compute ds */
    for(i = i1-1; i<=i2+1; i++)
        ds[i] = 0.5*(s1[i+1]-s1[i-1]);

    /* comput flux */
    for(i = i1; i<=i2+1; i++){
        float courant = fabs(dt/dx*c[i]);
        if(c[i]>=0)
            f[i] = courant*(s1[i-1]+(1.0-courant)/2.0*ds[i-1]);
        else
            f[i] = courant*(-s1[i]+(1.0-courant)/2.0*ds[i]);
    }
    
    for(i = i1; i<=i2; i++){
        s2[i] = -(f[i+1]-f[i]) + dt/dx*s1[i]*(c[i+1]-c[i]);
    }

	return;
}

