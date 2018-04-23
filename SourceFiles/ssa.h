/********************************************************************
   Tharuna Niranjan
   Kuhn Laboratory
   Dept. of Biochemistry
   Michigan State University

   ssa.h    Main header file for the automated SSA procedure  6.13.96

Makefile provided to compile all files belonging to this module
*********************************************************************/
#ifndef SSA_H
#define SSA_H

#define FILENAMELEN 100
#define PDBLINESIZE 81
#define NUMPOINTS 16
#define DIM 3 /* Number of dimensions */
#define ERROR 1
#define SUCCESS 0
#define BOUND 111111.1
#define THRESHOLD 0.75

/*Type definitions */
struct _sequerydata{
    char pdb_code[5];
    char chain_id;
    int from_res_no, to_res_no;
};


/*Function Prototypes */

int match_with_template(double a[NUMPOINTS][DIM], double c[NUMPOINTS][DIM], char *type);

int match_with_template_and_user(double a[NUMPOINTS][DIM], double c[NUMPOINTS][DIM], char *type, double *UserTemplateRMSD);

int process_pdbfile(struct _sequerydata *data, double a[NUMPOINTS][DIM], double
c[NUMPOINTS/2][DIM]);

int check_lr_helix(struct _sequerydata *shifted, char *type, int (*func)(double [][DIM],double [][DIM],char *));

void eigen(double *a, double *b, double *c, int n);

void emsg(char *, int, int, char *, int);

void matvec(double *, double *, double *);

int rotate(double[][DIM], double [][DIM], double [], int, double[DIM][DIM]);

static crossprod(double *, double *, double *);

int rotlsqfit (double [][DIM], double [][DIM], double [][DIM], double [], int, double [DIM][DIM], double [DIM]);

#endif
