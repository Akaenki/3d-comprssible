/*
 * ========================= update ====================
 * Update: replace old values with new ones
 * We are not copying ghost points here.
 * Arguments:
 *
 *	s1,s2	real arrays	old, new data arrays
 *	i1,i2	integers	indices bounding array data
 *	nydim	integer		size of arrays with ghost
 *
 */
#include "param.h"

void update(float u1[][NYDIM][NZDIM],float u2[][NYDIM][NZDIM],float u3[][NYDIM][NZDIM],float v1[][NYDIM][NZDIM],float v2[][NYDIM][NZDIM],float v3[][NYDIM][NZDIM],float w1[][NYDIM][NZDIM],float w2[][NYDIM][NZDIM],float w3[][NYDIM][NZDIM],float p1[][NYDIM][NZDIM],float p2[][NYDIM][NZDIM],float p3[][NYDIM][NZDIM],float t1[][NYDIM][NZDIM],float t2[][NYDIM][NZDIM])
{
	int i,j,k;

	for(i = I1; i <= I2; i++){
        for(j = J2; j <= J2; j++){
            for(k = K1; k <= K2; k++){
                p1[i][j][k] = p2[i][j][k];
                p2[i][j][k] = p3[i][j][k];
                t1[i][j][k] = t2[i][j][k];
            }
        }
	}
    
    for(i = I1; i <= I2+1; i++){
        for(j = J1; j <= J2; j++){
            for(k = K1; k <= K2; k++){
                u1[i][j][k] = u2[i][j][k];
                u2[i][j][k] = u3[i][j][k];
            }
        }
    }
    
    for(i = I1; i <= I2; i++){
        for(j = J1; j <= J2+1; j++){
            for(k = K1; k <= K2; k++){
                v1[i][j][k] = v2[i][j][k];
                v2[i][j][k] = v3[i][j][k];
            }
        }
    }
    
    for(i = I1; i <= I2; i++){
        for(j = J1; j <= J2; j++){
            for(k = K1; k <= K2+1; k++){
                w1[i][j][k] = w2[i][j][k];
                w2[i][j][k] = w3[i][j][k];
            }
        }
    }

	return;
}

