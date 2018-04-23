#include "frac.h"
#include "q.h"
#define GRAPHHID "@(#) $Header: /riscfs2/mp/grip/grinch/src/RCS/graph.h,v 1.1 90/01/03 18:50:10 arvai Exp $" /* included in graph.c */
#ifndef MYMACS
#include "mymacs.h"
#endif

#define MAXNODES 64000

 /* Note MAXEDGES must not exceed range of MAXMARKVAL below */
 /* Note MAXEDGES must be >= MAXNODES for 'marks' to work correctly */
#ifdef MKSKEL
#   define MAXEDGES 65000 /* 65000 AT UNC */
#endif
#ifndef MKSKEL
#   define MAXEDGES 65000
#endif

#define MAXDEGREE 5
#define MAXLEVEL 31
#define COORDSCALE 100	/* coord units are 1/SCALE Angstroms */
#define DEFAULTWIDTH (10.*COORDSCALE)
typedef int NODE;
typedef int EDGE;
typedef short COORDS[3];

#define MAXMARKS 7
typedef unsigned short MARKVAL; /* Technically, union of NODE and EDGE,*/
				/*  but this saves space */
#define MAXMARKVAL (MAXUSHORT)

typedef struct MMMARK{
	MARKVAL low;
	MARKVAL high;
	short in_use; /* used as boolean only */
	MARKVAL marks[MAXEDGES];
} MMARK;
typedef MMARK *MARK;
extern MMARK secretmark[MAXMARKS];
extern MARK create_mark();

#ifndef mc500
struct flagdefs {
	unsigned mode : 1;
	unsigned exists : 1;
	unsigned is_clipped : 1;
	unsigned filler : 32-3; /* For VMS C compiler compatibility */
};
#else
	/* bypass mc500 C compiler bug : */
struct flagdefs {
	char mode;
	char exists;
	char is_clipped;
};
#endif

typedef struct {
	struct flagdefs flags[MAXEDGES];
	unsigned short ends[MAXEDGES][2]; /* Technically NODE */
	int label[MAXEDGES];
	int label2[MAXEDGES];
	char type[MAXEDGES];
	char twist[MAXEDGES];
	int num;
	int changed;
} EDGESTYPE;

typedef struct {
	unsigned short adj[MAXNODES][MAXDEGREE]; /* Technically EDGE */
	COORDS coords[MAXNODES];
	short residue[MAXNODES]; /* internal residue id, see residue.h */
	char den[MAXNODES];
#ifndef mc500
	struct {
		unsigned exists : 1;
		unsigned is_clipped : 1;
		unsigned filler : 32-2; /* for VMS C compatibility */
	} flags[MAXNODES];
#else
	/* bypass mc500 C compiler bug */
	struct {
		char exists;
		char is_clipped;
	} flags[MAXNODES];
#endif
	int num;
} NODESTYPE;

typedef struct {
	COORDS center;
	int width;
	char denlevel;
	char softlev;
	fraction matrix[3][3];
	COORDS cursor;
} VIEW;

extern EDGESTYPE *secedgep; /* See initgraph in graph.c */
extern NODESTYPE *secnodep; /* for these allocations */
#define secretedge (*secedgep)
#define secretnode (*secnodep)
extern VIEW view;

#define FOR_3(i) {for(i=0; i<3; i++) {
#define FOR_EDGES(e) {for(e=1; e<=secretedge.num; e++) { \
	if(!secretedge.flags[e].exists) continue;
#define FOR_EDGES_BEFORE(b,e) {for(e=1; e<b&&e<=secretedge.num; e++) { \
	if(!secretedge.flags[e].exists) continue;
#define FOR_EDGES_AFTER(a,e) {for(e=(a<1)?1:a; e<=secretedge.num; e++) { \
	if(!secretedge.flags[e].exists) continue;
#define FOR_NODES(n) {for(n=1; n<=secretnode.num; n++) { \
	if(!secretnode.flags[n].exists) continue;
#define FOR_ADJ_E(n,e) {register int _nb_nb_; for(_nb_nb_=0; \
	_nb_nb_<MAXDEGREE && (e=secretnode.adj[n][_nb_nb_]) != NULL; _nb_nb_++){
#define FOR_ADJ_N(e,n) {register int _nb_nbn_; for(_nb_nbn_=0; \
	_nb_nbn_<=1; _nb_nbn_++) { n=secretedge.ends[e][_nb_nbn_];
#define FOR_E_NBRS(e,e2) \
	{register int _nb_1,_nb_2, _nb_3; \
	for(_nb_1=0; _nb_1<=1; _nb_1++) \
		for(_nb_2=secretedge.ends[e][_nb_1],_nb_3=0; _nb_3<MAXDEGREE &&\
			(e2=secretnode.adj[_nb_2][_nb_3]) != NULL; _nb_3++) {\
				if (e2 == e) continue;
	
#define FOR_END }}

#define DO_BFS(n, e, n2) \
	{static QUEUE bfsq; MARK bfsmark; NODE bfsnode;\
		bfsmark = create_mark(1); \
		set_mark(bfsmark, n); \
		q_init(&bfsq, 1000); \
		q_put(bfsq, n); \
		while (! q_empty(bfsq)) { \
			bfsnode = q_get(bfsq); \
			FOR_ADJ_E(bfsnode, e) \
			    FOR_ADJ_N(e, n2) \
				if (is_markset(bfsmark, n2)) \
					continue; \
				else    set_mark(bfsmark, n2);
#define END_BFS(n2) \
				q_put(bfsq, n2); \
			    FOR_END; \
			FOR_END; \
		} \
		del_mark(bfsmark); \
	}

#define set_changed()	(secretedge.changed = 1);
#define reset_changed()	(secretedge.changed = 0);
#define is_changed()	(secretedge.changed != 0)

#define ismodel(e)	(secretedge.flags[e].mode)
#define ismap(e)	(! ismodel(e))

#define nodevalid(n) (n>=0 && n<=secretnode.num && secretnode.flags[n].exists)
#define edgevalid(e) (e>=0 && e<=secretedge.num && secretedge.flags[e].exists)

#define set_label(e, val)	(secretedge.label[e] = val)
#define get_label(e)		(secretedge.label[e])
#define set_label2(e, val)	(secretedge.label2[e] = val)
#define get_label2(e)		(secretedge.label2[e])

  /* See 'marks' comments in mark.c. */
#define set_markval(m, i, val)	if (val >= 0 && val <= m->high) \
					m->marks[i] = val+m->low; \
				else { \
		error( "mark %d not between 0 and max range\n",val); \
					exit(-1); \
				}
#define get_mark(m, i)		((m)->marks[i]-(m)->low)
#define low_mark(m)		((m)->low)
#define high_mark(m)		((m)->high-(m)->low)
#define set_mark(m, i)		((m)->marks[i] = (m)->low)
#define is_markset(m, i)	((m)->marks[i] >= (m)->low \
			&&(m)->marks[i] <= (m)->high)
#define set_end(e, end, val)	(secretedge.ends[e][end] = val)
#define get_end(e, end)		(secretedge.ends[e][end])
#define other_end(e, n)		(secretedge.ends[e][0]==n ? \
				secretedge.ends[e][1] : secretedge.ends[e][0])
#define set_residue(n, val)	(secretnode.residue[n] = val)
#define get_residue(n)		(secretnode.residue[n])
#define set_nodeden(n, val)	(secretnode.den[n] = val)
#define get_nodeden(n)		(secretnode.den[n])
#define set_twist(e, val)	(secretedge.twist[e] = val)
#define get_twist(e)		(secretedge.twist[e])

#define EDGEDEN_IS_MAX /* if this is defined, MAX den not AVERAGE den is used*/
#ifdef EDGEDEN_IS_MAX
#define get_edgeden(e) \
	max(get_nodeden(get_end(e,0)), get_nodeden(get_end(e,1)))
#else
#define get_edgeden(e)\
	( (get_nodeden(get_end(e,0)) + \
	   get_nodeden(get_end(e,1))) / 2 )
#endif
#define set_clipped(e, val)	(secretedge.flags[e].is_clipped = val)
#define get_clipped(e)		(secretedge.flags[e].is_clipped)
#define is_visible(e)		(!get_nodeclipped(get_end(e,0)) && \
				 !get_nodeclipped(get_end(e,1)) && \
				  (ismodel(e) || \
				  get_edgeden(e) >= view.softlev))
#define is_nodevis(n)		(!get_nodeclipped(n) && \
				  get_nodeden(n) >= view.softlev)
#define set_nodeclipped(n, val)	(secretnode.flags[n].is_clipped = val)
#define get_nodeclipped(n)		(secretnode.flags[n].is_clipped)
#define set_type(e, val)	(secretedge.type[e] = val)
#define get_type(e)		(secretedge.type[e])
#define set_mode(e, val)	(secretedge.flags[e].mode = val)
#define get_mode(e)		(secretedge.flags[e].mode)

#define set_coords(n, c)    {register short *cp1, *cp2; \
    cp2 = c; cp1 = secretnode.coords[n]; \
    *cp1++ = *cp2++; *cp1++ = *cp2++; *cp1 = *cp2; }

#ifndef CMSWATC
#define get_coords(n, c) 	{register short *cp1, *cp2; \
	cp1 = c; cp2 = secretnode.coords[n]; \
	*cp1++ = *cp2++; *cp1++ = *cp2++; *cp1 = *cp2; }
#else
#define get_coords(n,c) {register int d;for(d=0;d<3;d++)\
   c[d]=\
   secretnode.coords[n][d];}
#endif
#define get_xcoord(n) (secretnode.coords[n][0])
#define get_ycoord(n) (secretnode.coords[n][1])
#define get_zcoord(n) (secretnode.coords[n][2])

extern MARK freedchain;
extern MARK substratemark;
extern MARK subsurfmark;
