#include "mymacs.h"
#include "graph.h"

/* Finds rotation matrix r and vector t such that	*/
/* SUM ( sqr (w* (b - r*a + t))) is minimized, where	*/
/* a and b are vectors of m points, weighted by w.	*/
/* These transformations are then applied to a to produce c */
/* By Tom Williams. Return of "vector" and "matrix" added by Mike Pique */

 int
rotlsqfit (a, b, c, w, m, matrix, vector)
register double a[][3];		/* input:  m fixed  3d points (double) */
register double b[][3];		/* input:  m guide  3d points (double) */
register double c[][3];		/* output: m result 3d points (double) */
double w[];			/* input:  m weights */
register int m;			/* input:  number of points */
register double matrix[3][3];	/* output: 3x3 rotation matrix 'r' */
register double vector[3];	/* output: 1x3 translation vector 't' (post-rotation) */
{
	double r[3][3];			/* 3x3 matrix to rotate a onto b */
	double ta[3];			/* translates origin to a */
	double tb[3];			/* translates origin to b */
	register int i, j;
	register double *newa;	/* a translated to origin */
	register double *newb;	/* b translated to origin */
	double wsum;
#define as(i,j) newa[i*3+j]
#define bs(i,j) newb[i*3+j]
	double *calloc();

#ifdef DEBUG
printf("rotlsqfit(%d points:\n",m);
for(i=0;i<m;i++) /* dump arrays   */ {
  printf(
  "pnt %d  weight %.2f fixed(%.1f %.1f %.1f) guide(%.1f %.1f %.1f)\n"
    ,i, w[i],
  a[i][0], a[i][1],  a[i][2],
  b[i][0], b[i][1],  b[i][2]);
	}
fflush(stdout);
#endif
	/* find center of mass of a and b */
	FOR_3(j)
		ta[j] = tb[j] = 0.0;
	FOR_END;
	wsum = 0.0;

	for (i=0; i<m; i++) {
		FOR_3(j)
			ta[j] += w[i] * a[i][j];
			tb[j] += w[i] * b[i][j];
		FOR_END;
		wsum += w[i];
	}
	if(wsum == 0) wsum=1;
	FOR_3(j)
		ta[j] = ta[j] / wsum;
		tb[j] = tb[j] / wsum;
	FOR_END;

	newa = calloc (m*3, sizeof (double));
	newb = calloc (m*3, sizeof (double));
	for (i=0; i<m; i++)
		FOR_3(j)
			as(i,j) = a[i][j] - ta[j];
			bs(i,j) = b[i][j] - tb[j];
		FOR_END;

    if (! rotate (newa, newb, w, m, r))   {
#ifdef DEBUG
	dmatdump(r,m,m,"rotlsqfit - rotate failed,  matrix =");
	     free(newa); free(newb);
#endif
     return(FALSE);
     }
#ifdef DEBUG
dmatdump(r,m,m,"rotlsqfit - rotate returned matrix =");
#endif
     free(newa); free(newb);
	if(c!=NULL) for (i=0; i<m; i++) {
		FOR_3(j)
			c[i][j] = a[i][j] - ta[j];
		FOR_END;
		matvec(r, c[i], c[i]);
		FOR_3(j)
			c[i][j] += tb[j];
		FOR_END;
	}
	if(matrix!=NULL) FOR_3(i) FOR_3(j) matrix[i][j] = r[i][j]; FOR_END; FOR_END;
	if(vector!=NULL) FOR_3(j) 
		vector[j] = tb[j] - (r[j][0]*ta[0] + r[j][1]*ta[1] + r[j][2]*ta[2]);
		FOR_END;
	return(TRUE);
}
