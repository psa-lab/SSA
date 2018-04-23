#include "mymacs.h"
#include <stdio.h>
#include <math.h>
#ifdef CMSWATC
extern double fabs(); /*CMSWATC */
#endif

#ifndef HUGE
#define HUGE 1e10    /* large positive number */
#endif
typedef double *matrix;
typedef double *vector;

 static
tred2 (a, z, d, e, tol, n)
matrix a;
#define as(i, j) (a[(i)*n + (j)])
register matrix z;
#define zs(i, j) (z[(i)*n + (j)])
vector d;
register vector e;
double tol;
int n;

/*	algorithm from NUMERISCHE MATHEMATIK 11, 181-195 (1968)	*/
/*	for Householder tridiagonalization of a, transformation		*/
/*	matrix is returned in z, diagonal in d, and subdiagonal		*/
/*	in e, tol is machine dependent.					*/

{ /**/
	register int i, j, k;
	int L;
	register double g, h, f, hh;

#ifdef DEBUG
printf("tred2( 0x%x, 0x%x, 0x%x, 0x%x, tol=%g, n=%d)\n",
a, z, d, e, tol, n);     fflush(stdout);
dmatdump(a, n, n, " tred2 ain= ");

#endif

	/* z = a */
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			zs(i, j) = as(i, j);

	for (i=n-1; i>=1; i--) {
		L = i - 2;
		f = zs(i, i-1);
		g = 0.0;
#ifndef CMSWATC
		for (k=0; k<=L; k++) g += (zs(i, k) * zs(i, k));
#else
        for (k=0; k<=L; k++) g=g+ (zs(i, k) * zs(i, k)); /*CMSWATC was +=*/
#endif
		h = g + f * f;

		/* skip this transformation if g is too small */
		if (g <= tol) {
			e[i] = f;
			h = 0.0;
		}
		else {
			L++;
			if (f >=0.0)	g = -sqrt(h);
			else		g = sqrt(h);
			e[i] = g;
#ifndef CMSWATC
			h -= (f*g);
#else
            h= h-(f*g);
#endif
			zs(i, i-1) = f - g;
			f = 0.0;
			for (j=0; j<=L; j++) {
				zs(j, i) = zs(i, j) / h;
				g = 0.0;
				for (k=0; k<=j; k++)
#ifndef CMSWATC
					g += (zs(j, k) * zs(i, k));
#else
                    g=g+ (zs(j, k) * zs(i, k)); /*CMSWATC was += */
#endif
				for (k=j+1; k<=L; k++)
#ifndef CMSWATC
					g += (zs(k, j) * zs(i, k));
#else
                    g =g+(zs(k, j) * zs(i, k));  /* CMSWATC was += */
#endif
				e[j] = g/h;
#ifndef CMSWATC
				f += (g * zs(j, i));
#else
                f=f+ (g * zs(j, i));
#endif
			}
			hh  = f / (h + h);
			for (j=0; j<=L; j++) {
				f = zs(i, j);
				g = e[j] - hh * f;
				e[j] = g;
				for (k=0; k<=j; k++)
#ifndef CMSWATC
					zs(j, k) -= (f*e[k] + g*zs(i, k));
#else
                 zs(j, k) =zs(j, k)-(f*e[k] + g*zs(i, k));
#endif
			}
		}
		d[i] = h;
	}
	d[0] = 0.0;
	e[0] = 0.0;
	for (i=0; i<n; i++) {
		L = i - 1;
		if (d[i] != 0.0)
			for (j=0; j<=L; j++) {
				g = 0.0;
#ifndef CMSWATC
				for (k=0; k<=L; k++) g += (zs(i, k) * zs(k, j));
				for (k=0; k<=L; k++) zs(k, j) -= (g*zs(k, i));
#else
                for (k=0; k<=L; k++) g=g+ (zs(i, k) * zs(k, j));
                for (k=0; k<=L; k++) zs(k, j) = zs(k, j)-(g*zs(k, i));
#endif
			}
		d[i] = zs(i, i);
		zs(i, i) = 1.0;
		for (j=0; j<=L; j++) {
			zs(i, j) = 0.0;
			zs(j, i) = 0.0;
		}
	}
}

 static
tql2 (z, d, e, eps, n)
matrix z;
#define zs(i, j) (z[(i)*n + (j)])
register vector d;
register vector e;
double eps;
int n;

/*	Algorithm from NUMERISCHE MATHEMATIK 11, 293-306 (1968)		*/
/*	For QL eigenvalue procedure for tridiagonal matrices.		*/
/*	This procedure provides a virtually orthonormal set of		*/
/*	eigenvectors.  When z is the output of tred2, the eigen-	*/
/*	vectors are those of the original matrix a.  When z is		*/
/*	the identity matrix the eigenvectors are those of the 		*/
/*	tridiagonal matrix whose diagonal elements are given in		*/
/*	d and subdiagonal elements in e.				*/
/*	Eigenvalues are returned in d and eigenvectors in z, by		*/
/*	columns, ordered in increasing order of eigenvalues.		*/

{
register int i, L, j;
int k, m;
register double h, p, r, b, c, f, g, s;

#ifdef DEBUG
printf("tql2( 0x%x, 0x%x, 0x%x, eps=%g n=%d)\n",
z, d, e, eps, n);
 dmatdump(z, n, n, " tql2 zin = ");

 fflush(stdout);
#endif
	for (i=1; i<n; i++)
		e[i-1] = e[i];
	e[n-1] = 0.0;
	b = 0.0;
	f = 0.0;
	for (L=0; L<n; L++) {
		j = 0;
		h = eps * (fabs(d[L]) + fabs(e[L]));
		if (b < h) b = h;

		/* look for small sub-diagonal element */
        for (m=L; fabs(e[m]) > b  && m<n ; m++)
			;
#ifdef DEBUG
       /*DEBUG*/printf("tql2 L=%d h=%g b=%g f=%g  m=%d\n", L, h, b, f, m);
   fflush(stdout);/*DEBUG*/
#endif
		if (m != L)
			do {
				if (j >= 30) {
					error("tql2:no convergence\n", j);
					return; /* no easy recovery */
				}
				j++;

				/* form shift */
				p = (d[L+1] - d[L]) / (2.0 * e[L]);
				r = sqrt(p*p + 1.0);
				if (p < 0.0)	h = p - r;
				else		h = p + r;
				h = d[L] - e[L]/h;
				for (i=L; i<n; i++)
#ifndef CMSWATC
					d[i] -= h;
				f += h;
#else
                    d[i]  = d(:i:)-h;
                f=f+ h;
#endif

				/* QL transformation */
				p = d[m];
				c = 1.0;
				s = 0.0;
				for (i=m-1; i>=L; i--) {
					g = c*e[i];
					h = c*p;
					if (fabs(p) >= fabs(e[i])) {
						c = e[i] / p;
						r = sqrt(c*c + 1.0);
						e[i+1] = s * p * r;
						s = c / r;
						c = 1.0 / r;
					}
					else {
						c = p / e[i];
						r = sqrt(c*c + 1.0);
						e[i+1] = s * e[i] * r;
						s = 1.0 / r;
#ifndef CMSWATC
						c *= s;
#else
						c =c*s;
#endif
					}
					p = c*d[i] - s*g;
					d[i+1] = h + s * ( c*g + s*d[i]);

					/* form vector */
					for (k=0; k<n; k++) {
						h = zs(k, i+1);
						zs(k, i+1) = s * zs(k, i) + c*h;
						zs(k, i) = c*zs(k, i) - s*h;
					}
				}
				e[L] = s * p;
				d[L] = c * p;
			} while (fabs(e[L]) > b);
#ifndef CMSWATC
		d[L] += f;
#else
        d[L] =d(:L:) + f;
#endif
	}
#ifdef DEBUG
dmatdump(d, 1, n, "(tql2) eigenvalues before sorting");
#endif

	/* order eigenvalues and eigenvectors */
	for (i=0; i<n; i++) {
		k = i;
		p = d[i];
		for (j=i+1; j<n; j++)
			if (d[j] < p) {
				k = j;
				p = d[j];
			}
			if (k != i) {
				d[k] = d[i];
				d[i] = p;
				for (j=0; j<n; j++) {
					p = zs(j, i);
					zs(j, i) = zs(j, k);
					zs(j, k) = p;
				}
			}
	}
#ifdef DEBUG
dmatdump(d, 1, n, "(tql2) eigenvalues after  sorting");
#endif
}

dmatdump(mat, m, n, label)
double *mat; int m, n; char *label;
{ /* dump m by n matrix of doubles */
int i, j;
printf("%s\n", label);
for(i=0;i<m;i++) {
  for(j=0;j<n;j++) printf(" %9g", mat[i*n+j]);
  printf("\n");
  }
}

eigen(a, b, c, n)
register matrix a;	/* input matrix (n x n) */
register matrix b;	/* output eigenvectors in columns */
register vector c;	/* output eigenvalues */
register int n;
{
	register vector e;
	vector calloc();
	double epsilon;
#ifdef vax
	union {
		double d;
		struct {
			long high;
			long low;
		} h;
    } eps;
	union {
		double d;
		struct {
			long high;
			long low;
		} h;
	} TAU;
	/* Set up bit patterns on VAX for optimum performance,
	 * See VAX Architecture Handbook Section 4.5 and 4.6
	 */
    eps.h.high = 0xffff24ff;
    eps.h.low  = 0xffffffff;
    epsilon = eps.d;
	TAU.h.high = 0x00000080; /* 0 x 2**-127 on VAX */
	TAU.h.low  = 0x00000000;
#else
    union { double d; } eps;
	union { double d; } TAU;
#ifdef CMSWATC
 /* the math library defines EPS and INF */
extern double EPS, INF;
 epsilon = /* EPS ==5.39E-79... too small*/ 2.6E-10;
 TAU.d = 1./INF;
#else

  epsilon =  2.6E-17; /* Set this to precision of a double-precision
				* floating point value */
	TAU.d = 1. / HUGE; /* HUGE is defined in <math.h>, TAU is
supposed to be  smallest possible positive value */
#endif /* not CMSWATC */
#endif /* not vax */

#ifdef DEBUG
printf("eigen(0x%x, 0x%x, 0x%x, %d)\n", a, b, c, n);  fflush(stdout);
printf("TAU = %g    eps = %g   INF=%g\n", TAU.d, epsilon, INF);
#endif
	e = calloc (n, sizeof e[0]);
	if (e == NULL) {
		error("alloc fail", n);
		return; /* no easy recovery */
	}
	tred2 (a, b, c, e, TAU.d, n);
	tql2 (b, c, e, epsilon, n);
	free (e);
}


