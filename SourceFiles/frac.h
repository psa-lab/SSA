#ifndef frmul

  /* frac.h: 'fraction' type definitions. Uses float */
  /* $Header: /riscfs2/mp/grip/grinch/src/RCS/frac.h,v 1.1 90/01/03 18:50:02 arvai Exp $ */
  
typedef float fraction;

#ifdef vax
typedef double fastfloat; /* the current Vax C compilers do not have
			* a fast single-precision option.
			*/
#else
typedef float fastfloat; /* but we assume other machines do, and use it.
			* this assumes the code is compiled to use it
			* (sun: -fsingle, masscomp -q)
			*/
#endif

#define frmul(a,b) (((fastfloat)(a))*((fastfloat)(b)))
#define frdiv(a,b) (((fastfloat)(a))/((fastfloat)(b)))
#define fradd(a,b) (((fastfloat)(a))+((fastfloat)(b)))
#define frsub(a,b) (((fastfloat)(a))-((fastfloat)(b)))

#define frimul(a,b) ((int)(frmul(a,b)))

#define itofr(i) ((fastfloat)(i))
#define ltofr(l,s) ((fastfloat)(l)) /* not right */
#define frtoflt(fr) (fr)
#define flttofr(flt) ((float) (flt))

double sin(), cos();
#define frsin(x) ((double)sin((double)(x)))
#define frcos(x) ((double)cos((double)(x)))
#endif frmul
