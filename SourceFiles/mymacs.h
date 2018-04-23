#ifndef MYMACS
#define MYMACS /* block multiple includes of this file */

#include <stdio.h>
#define TRUE	1
#define FALSE	0

#define MAXSHORT	((1<<15) - 1)
#define MAXUSHORT	((1<<16) - 1)
#define MAXINT		((1<<31) - 1)

#define NULLCHAR	'\0'

#define NO_STEREO	0
#define STEREO_SIDE	1
#define STEREO_TIME	2
#define STEREO_ROCK	3
#define STEREO_ORTHO	4

#define local static
#define max(a, b)	((a) > (b) ? (a) : (b))
#define min(a, b)	((a) < (b) ? (a) : (b))
#define abs(a)		((a) < 0 ? (-(a)) : (a))
#define asize(a)	((sizeof(a)) / (sizeof(a[0])))

#ifdef CMSWATC
#define DUMPHEADER	" a"
#else
#define DUMPHEADER	"in file %s, line %d: a", __FILE__, __LINE__
#endif CMSWATC

#define __DUMP(a, fmt)						\
	{							\
		register int i;					\
								\
		fprintf(stderr, DUMPHEADER);			\
		for (i = 0; i < asize(a); ++i)			\
			fprintf(stderr, "\t%fmt", a[i]);	\
		fprintf(stderr, "\n");				\
	}

#define DUMPARRAY(a)	__DUMP(a, d)
#define DUMPARRAYF(a)	__DUMP(a, g)

#define streq(s1, s2)	(!strcmp(s1, s2))

#define word short
#define boolean int

#define FOR_ONCE				\
	{					\
		static int _DONE_ONCE = FALSE;	\
						\
		if (!_DONE_ONCE)		\
		{				\
			_DONE_ONCE = TRUE;

#define FOR_END	}}

#define vcpy(a, b) (a)[0] = (b)[0], (a)[1] = (b)[1], (a)[2] = (b)[2]

#ifdef CMSWATC
#define error(mess, val)	emsg(mess, val, 0, "", 0)
#define DEBUGMSG(mess, val)	fprintf(stderr,	\
					"debug msg: %s %d\n", mess, val)
#else 
#define error(mess, val)	emsg(mess, val, 0, __FILE__, __LINE__)
#define DEBUGMSG(mess, val)	emsg(mess, val, 1, __FILE__, __LINE__)
#endif CMSWATC

#endif MYMACS
