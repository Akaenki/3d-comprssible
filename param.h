//
//  param.h
//  pgm6
//
//  Created by Linling  Miao on 4/25/17.
//  Copyright © 2017 Linling Miao. All rights reserved.
//

#ifndef param_h
#define param_h

/* Global Constants */
#define NX 201
#define NY 201
#define NZ 41
#define BC_WIDTH 2
#define I1 BC_WIDTH
#define I2 I1+NX-1
#define J1 BC_WIDTH
#define J2 J1+NY-1
#define K1 BC_WIDTH
#define K2 K1+NZ-1
#define NXDIM NX+2*BC_WIDTH
#define NYDIM NY+2*BC_WIDTH
#define NZDIM NZ+2*BC_WIDTH
#define TBAR 300
#define CS 100

/* Global variables */
float dt,dx,dy,dz,km,kmt;
float xc[2],yc[2],zc[2];
float xradius[2],yradius[2],zradius[2],dtp[2];
float du[NZ],dw[NZ+1];


#endif /* param_h */
