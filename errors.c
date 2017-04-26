/*
* ======================== errors =====================
* Calculte the Takacs errors
* Arguments:
*
*	s1	real array	IC data. Set i1..i2 here;
*				  [i1-1],[i2+1] = ghost zones
*				  if 1 ghost point on each side
*	strue	real	array	True solution
*	i1,i2,j1,j2	integers	indices bounding array data
*	nydim	real	number of grid points including ghost
*/
#include <math.h>

void errors(int nx, int ny, int nydim, float s1[][nydim], float strue[][nydim], 
		int i1, int i2, int j1, int j2)
{
	int i, j;
	float s1_avg = 0, strue_avg = 0, s1_var = 0, strue_var = 0;

	/* Calculate Averages */
	for (i = i1; i <= i2; i++) {
		for (j = j1; j <= j2; j++) {
			s1_avg += s1[i][j] / (nx*ny);
			strue_avg += strue[i][j] / (nx*ny);
		}
	}

	/* Calculate Variances */
	for (i = i1; i <= i2; i++) {
		for (j = j1; j <= j2; j++) {
			s1_var += (s1[i][j] - s1_avg)*(s1[i][j] - s1_avg)/nx/ny;
			strue_var += (strue[i][j] - strue_avg)*(strue[i][j] - strue_avg)/nx/ny;
		}
	}
	/* Calculate the Correlation Coefficient */
	float numerator = 0;
	for (i = i1; i <= i2; i++) {
		for (j = j1; j <= j2; j++) {
			numerator += (s1[i][j] - s1_avg)*(strue[i][j] - strue_avg);
		}
	}
	float denominator = sqrt(s1_var*nx*ny*strue_var*nx*ny);
	float corr = numerator / denominator;
	printf("Correlation Coefficient = %.5f\n",corr);

	/* Calculate Takacs Errors */
	float e_total = 0;
	for (i = i1; i <= i2; i++) {
		for (j = j1; j <= j2; j++) {
			e_total += (strue[i][j] - s1[i][j])*(strue[i][j] - s1[i][j])/nx/ny;
		}
	}

	float e_diss = (sqrt(strue_var) - sqrt(s1_var))*(sqrt(strue_var) - sqrt(s1_var))
		+ (strue_avg - s1_avg)*(strue_avg - s1_avg);
	float e_disp = 2.0*(1.0 - corr)*sqrt(strue_var)*sqrt(s1_var);
	printf("e_total: %.5f, e_diss: %.5f, e_disp: %.5f\n", e_total, e_diss, e_disp);

	return;
}