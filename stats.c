/*
 * ========================= stats =====================
 * Stats computes and prints out the max S values
 * Arguments:
 *
 *	s2	real array	Latest data. Check I1..I2;
 *				  [I1-1],[I2+1] = ghost zones
 *	I1,I2,J1,J2	integers	indices bounding array data
 *	nydim	integer		size of data array with ghost
 *	n	integer		time step counter
 *	smax	real		holds max value of s2
 *	smin	real		holds min value of s2
 */

#include "param.h"

void stats(float s[][NYDIM][NZDIM],float *smax,float *smin)
{
	int i,j,k;
    int imax = 0,jmax = 0,kmax = 0;

	float max = s[I1][J1][K1];
	float min = s[I1][J1][K1];
	for (i = I1; i <= I2; i++){
		for (j = J1; j <= J2; j++){
            for(k = K1; k <= K2; k++){
                if (s[i][j][k] > max){
                    max = s[i][j][k];
                    imax = i; jmax = j; kmax = k;
                }
                else if (s[i][j][k] < min) min = s[i][j][k];
            }
		}
	}
    
	*smax = max;
	*smin = min;

	return;
}

