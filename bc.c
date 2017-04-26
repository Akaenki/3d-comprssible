/*
 * This routine applies boundary conditions
 */
#include "param.h"

void bc(float t[][NYDIM][NZDIM],float p[][NYDIM][NZDIM],float u[][NYDIM][NZDIM],float v[][NYDIM][NZDIM],float w[][NYDIM][NZDIM])
{
    int i,j,k;
    
    /* top-bottom */
    for(i = I1; i<=I2; i++){
        for(j = J2; j<=J2; j++){
            for(k = 0; k<=K1; k++){
                w[i][j][K1-k] = 0.0;
                w[i][j][K2+1+k] = 0.0;
            }
        }
    }
    for(i = I1; i<=I2; i++){
        for(j = J1; j<=J2+1; j++){
            for(k = 1; k<=K1; k++){
                v[i][j][K1-k] = v[i][j][K1];
                v[i][j][K2+k] = v[i][j][K2];
            }
        }
    }
    for(i = I1; i<=I2+1; i++){
        for(j = J1; j<=J2; j++){
            for(k = 1; k<=K1; k++){
                u[i][j][K1-k] = u[i][j][K1];
                u[i][j][K2+k] = u[i][j][K2];
            }
        }
    }
    
	for(i = I1; i<=I2; i++) {
        for(j = J1; j<=J2; j++){
            for(k = 1; k<=K1; k++) {
                t[i][j][K1-k] = t[i][j][K1];
                t[i][j][K2+k] = t[i][j][K2];
                p[i][j][K1-k] = p[i][j][K1];
                p[i][j][K2+k] = p[i][j][K2];
            }
        }
	}
    
    /* left-right */
    for(k = K1; k<=K2; k++){
        for(j = J1; j<=J2; j++){
            for(i = 0; i<=I1; i++){
                u[I1-i][j][k] = -u[I1+i+1][j][k];
                u[I2+i+1][j][k] =  -u[I2-i][j][k];
            }
        }
    }
    
    for(k = K1; k<=K2; k++){
        for(j = J1; j<=J2+1; j++){
            for(i = 1; i<=I1; i++){
                v[I1-i][j][k] = -v[I1+i+1][j][k];
                v[I2+i+1][j][k] =  -v[I2-i][j][k];
            }
        }
    }
    
    for(k = K1; k<=K2+1; k++){
        for(j = J1; j<=J2; j++){
            for(i = 1; i<=I1; i++){
                w[I1-i][j][k] = w[I1+i][j][k];
                w[I2+i][j][k] = w[I2-i][j][k];
            }
        }
    }
    
    for(k = K1; k<=K2; k++){
        for(j = J1; j<=J2; j++){
            for(i = 1; i<=I1; i++){
                t[I1-i][j][k] = t[I1+i][j][k];
                t[I2+i][j][k] = t[I2-i][j][k];
                p[I1-i][j][k] = p[I1+i][j][k];
                p[I2+i][j][k] = p[I2-i][j][k];
            }
        }
    }
    
    /* front-back */
    for(k = K1; k<=K2; k++){
        for(i = I1; i<=I2; i++){
            for(j = 0; j<=J1; j++){
                v[i][J1-j][k] = -v[i][J1+j+1][k];
                v[i][J1+j+1][k] =  -v[i][J2-j][k];
            }
        }
    }
    
    for(k = K1; k<=K2; k++){
        for(i = I1; i<=I2+1; i++){
            for(j = 1; j<=J1; j++){
                u[i][J1-j][k] = -u[i][J1+j+1][k];
                u[i][J2+j+1][k] =  -u[i][J2-j][k];
            }
        }
    }
    
    for(k = K1; k<=K2+1; k++){
        for(i = I1; i<=I2; i++){
            for(j = 1; j<=J1; j++){
                w[i][J1-j][k] = w[i][J1+j][k];
                w[i][J2+j][k] = w[i][J2-j][k];
            }
        }
    }
    
    for(k = K1; k<=K2; k++){
        for(i = I1; i<=I2; i++){
            for(j = 1; j<=J1; j++){
                t[i][J1-j][k] = t[i][J1+j][k];
                t[i][J2+j][k] = t[i][J2-j][k];
                p[i][J1-j][k] = p[i][J1+j][k];
                p[i][J2+j][k] = p[i][J2-j][k];
            }
        }
    }
        
	return;
}

