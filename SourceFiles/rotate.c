/* Algorithm from appendix of paper by A. D. Mclachlan, J. Mol Biol (1979) */
/* vol. 128, pp 74-77.							   */
/* Given two vectors a[] and b[] of position vectors with weights w[i]     */
/* find a rotation matrix r which minimizes SUM (w[i] x (r X a[i] - b[i])) */

#include <stdio.h>
#include <math.h>
#ifdef CMSWATC
extern double fabs(); /* CMSWATC */
#endif

 static 
crossprod(z, x, y)
register double *x, *y, *z;
{
	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];
}


 int
rotate (a, b, w, m, r)
register double a[][3];		/* input:  m fixed 3d points (double) */
register double b[][3];		/* input:  m guide 3d points (double) */
double w[];			/* input:  m weights */
register int m;			/* input:  number of points */
register double r[3][3];	/* output: 3x3 matrix to rotate a onto b */
{
#define SMALL 1.0e-12
#define ERROR 0
#define OK 1

	double u[3][3];
	double evec[6][6];	/* columns are eigenvectors */
	double eval[6];		/* eigenvalues */
	double omega[6][6];	/* ( 0  U ) */
				/* ( U' 0 ) */
	register int i, j, p;
	double sum, tmp;
	double h[6][3];
	double k[6][3];
	double d[6];		/* eigenvalues */
	double det;
	int rank;
	double x[3];
	double sqrt2;

#ifdef DEBUG
printf("rotate( %d points:\n",m);fflush(stdout);
	for (i=0; i<m; i++) {
    printf("point %d weighted %g\n", i,w[i]);
    printf(" fixed %g %g %g\n",a[i][0],a[i][1],a[i][2]);
    printf(" guide %g %g %g\n",b[i][0],b[i][1],b[i][2]);
	}
	printf("\n");
#endif

	/* form u matrix where u[i][j] = SUM (w[p] x a[p][i] x b[p][j]) */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++) {
			u[i][j] = 0.0;
			for (p=0; p<m; p++)
				u[i][j] += w[p] * a[p][i] * b[p][j];
		}
#ifdef DEBUG
    dmatdump(u,3,3," (rotate) weighted product mat u =");
#endif

	/* form omega matrix */
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			omega[i][j] = 0.0;
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			omega[i][j+3] = u[i][j];
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			omega[i+3][j] = u[j][i];
	for (i=3; i<2*3; i++)
		for (j=3; j<2*3; j++)
			omega[i][j] = 0.0;

	eigen(omega, evec, eval, 2*3);

#ifdef DEBUG
dmatdump(evec,6,6,"(rotate) eigen returned evec=");
dmatdump(eval,1,6,"    eval=");
#endif


	/* extract h, k, d from eigenvectors and eigenvalues */
	sqrt2 = sqrt(2.0);
	for (i=0; i<6; i++)
		for (j=0; j<3; j++){
			h[i][j] = evec[j][5-i] * sqrt2;
			k[i][j] = evec[j+3][5-i] * sqrt2;
		}
	for(i=0; i<6; i++)
		d[i] = eval[5-i];

#ifdef DEBUG
    printf("h (extracted from evec) \n");
	for (i=0; i<6; i++) {
		for(j=0; j<3; j++)
            printf("%6.3g  ", h[i][j]);
		printf("\n");
	}
	printf("k\n");
	for (i=0; i<6; i++) {
		for(j=0; j<3; j++)
            printf("%6.3g  ", k[i][j]);
		printf("\n");
	}
    printf("d (eigenvalues returned by eigen)\n");
	for(j=0; j<6; j++)
        printf("%g  ", d[j]);
	printf("\n");
#endif

	crossprod(x, h[0], h[1]);

#ifdef DEBUG
	printf("h[0] X h[1]\n");
	for(j=0; j<3; j++)
        printf("%g  ", x[j]);
	printf("\n");
#endif

	/* set sign of h[2] so that h[2] = h[0] X h[1] */
	sum = 0.0;
	for (i=0; i<3; i++) {
		tmp = x[i] - h[2][i];
		sum += tmp * tmp;
	}
	if (sum > SMALL*SMALL) {
#ifdef DEBUG
        fprintf(stdout, "rotate changed sign sum=%g\n",sum);
#endif
		for (i=0; i<3; i++) {
			h[2][i] = - h[2][i];
			k[2][i] = - k[2][i];
		}
	}
	det =	u[0][0] * (u[1][1]*u[2][2] - u[1][2]*u[2][1]) +
		u[0][1] * (u[2][0]*u[1][2] - u[1][0]*u[2][2]) +
		u[0][2] * (u[1][0]*u[2][1] - u[2][0]*u[1][1]) ;
#ifdef DEBUG
	printf("det by brute force %g\n", det);
    dmatdump(u,3,3, " u mat was ");
#endif
	rank = 0;
	for (i=0; i<3; i++)
		if (fabs(d[i]) > SMALL) rank++;

	if (det < 0) {
		if (d[1] != d[2]) {
			/* reflection not allowed so  use k[2] = - k[2] */
#ifdef DEBUG
			fprintf(stdout, "reflection case\n");
#endif
			for (i=0; i<3; i++)
				k[2][i] = - k[2][i];
		}
        else {
#ifdef DEBUG
        printf("rotate returning err, less than 3 degrees freedom\n");
#endif
        return(ERROR); /* < 3 degrees of rotational freedom */
	}
    }
	else if (rank == 2) {
#ifdef DEBUG
		fprintf(stdout, "planar data\n");
#endif
		for (i=0; i<3; i++)
			h[2][i] = x[i];
		crossprod(k[2], k[0], k[1]); /* k[2] = k[0] X k[1] */
	}
    else if (rank < 2) {
#ifdef DEBUG
       printf("rotate returning error, rank %d less than 2\n",rank);
#endif
       return(ERROR); /* rank == 0 or 1 */
       }
	for (i=0; i<3; i++)
		for (j=0; j<3; j++) {
			r[i][j] = 0.0;
			for (p=0; p<3; p++)
				r[i][j] += k[p][i] * h[p][j];
			if (fabs(r[i][j]) < SMALL) r[i][j] = 0.0;
		}
	return(OK);

}

