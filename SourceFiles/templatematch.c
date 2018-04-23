/**********************************************************************
   Tharuna Niranjan
   Kuhn Laboratory
   Dept. of Biochemistry
   Michigan State University

   templatematch.c                              6.13.96

This function will take in the data from the sequery match and do a
template matching with the Beta, Helix, Turn Type 1, Turn Type 1P,
Turn Type 2, Turn Type 2P, Turn Type 8, by doing a least squares fit
and thresholding against a rmsd value. It returns the type of the data
as determined as well as an error code on encountering any abnormality
in the input data.  

Makefile provided to compile all files belonging to this module
************************************************************************/

/************************************************************************
Paul C. Sanschagrin
Kuhn Laboratory
Department of Biochemistry
Michigan State University
East Lansing MI 48824

Modifications to code described above made by PCS on 5-Dec-96

These modifications allow the function to assign the secondary structure
  type of a segment and to compare it to a user template and report the
  RMSD back.
For more info see ssa.c, README.ssa

Modifcations are as a modified version of match_with_template named 
  match_with_template_and_user

These modifications are demarcated in the code via my intitials (PCS) 
  and primarily consist of an additional received reference parameter, 
  assignemnt of RMSD to this parameter, and an alteration in the iteration
  schedule to enable the majority of the code to remain as well as to
  ensure the comparison of the user template
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ssa.h"

void get_rmsd (double a[NUMPOINTS][DIM], double b[NUMPOINTS][DIM], double c[NUMPOINTS][DIM], double w[NUMPOINTS], double *rmsd);

extern FILE *debugfile;
extern FILE *logfile;

extern char betafilename[FILENAMELEN], helixfilename[FILENAMELEN], tt1filename[FILENAMELEN], tt2filename[FILENAMELEN], tt1pfilename[FILENAMELEN], tt2pfilename[FILENAMELEN], tt8filename[FILENAMELEN] ;
/*----------Start of PCS----------*/
extern char userfilename[FILENAMELEN];
extern int INCLUDEUSERTEMPLATE;
/*-----------End of PCS-----------*/


extern int DEBUG;

int match_with_template(double a[][DIM], double x[][DIM], char *type)
{
double b[NUMPOINTS][DIM], c[NUMPOINTS][DIM], tmp[NUMPOINTS][DIM], w[NUMPOINTS];
double  rmsd, minrmsd[2];
int  i,j,k,iter =0;
char templatefilename[FILENAMELEN], pdbline[PDBLINESIZE];
char atom_type[4];
FILE *templatefile;

/*Initializing the index and value of the min rmsd value to keep tab of the turn types */
minrmsd[0] = THRESHOLD;  /*we are looking for a value that is less than 0.75*/
minrmsd[1] = 0;   /*This check is done only for turn types, not Beta or helix*/

/*Initializing weights to a value of 1*/
for(i=0; i<NUMPOINTS; i++)
    w[i] = 1;

/*Checking the validity of the input data; The anamolies could be that we could not grab the 3D coordinates for all 4 atoms of the backbone chain. If there is an anamoly, we simply return with an error message */
for(i=0;i<NUMPOINTS;i++)
    if (a[i][0] == BOUND || a[i][1] == BOUND || a[i][2] == BOUND) 
    {
	printf("ERROR: Backbone atoms from PDB file could not be matched with the template for superposition,\
probably due to missing or extra atoms. - ");
	strcpy(type,"ERROR");
	fprintf(logfile,"ERROR: Backbone atoms from PDB file could not be matched with the template for superposition,\
probably due to missing or extra atoms. - ");
	return ERROR;
    }
for(i=0;i<NUMPOINTS/2;i++)
    if(x[i][0] == BOUND || x[i][1] == BOUND || x[i][2] == BOUND)
    {
      printf("ERROR: Insufficient Input Data. Couldnt grab residues on the lower and upper boundary of the sequence. - ");
      strcpy(type,"ERROR");
      fprintf(logfile,"ERROR: Insufficient Input Data. Couldnt grab residues on the lower and upper boundary of the sequence. - ");
      return ERROR;
    }

while ( iter < 7)      /* iter is the total num of types that we are matching against. If any of the matches succeeds, we return*/
    {
    switch(iter)  {
      case 0: strcpy(templatefilename,betafilename);
              break;
      case 1: strcpy(templatefilename,helixfilename);
              break;
      case 2: strcpy(templatefilename,tt1filename);
              break;
      case 3: strcpy(templatefilename,tt2filename);
              break;
      case 4: strcpy(templatefilename,tt1pfilename);
              break;
      case 5: strcpy(templatefilename,tt2pfilename);
              break;
      case 6: strcpy(templatefilename,tt8filename);
              break;
    }

    /*Open the template file and Input the template residues' 3D coordinates*/
    if( (templatefile = fopen(templatefilename,"r")) == NULL)
    {
	printf("UNABLE to open template file %s - Returning without processing other sequery matches\n", templatefilename);
	fprintf(logfile,"UNABLE to open template file %s - Returning without processing other sequery matches\n", templatefilename);
	return ERROR;
    }

    /*Reading the 3D coordinates from the template file */
    i = 0;
    while (!(feof(templatefile)) && i<NUMPOINTS)
    {
	fgets(pdbline,PDBLINESIZE+1, templatefile);
	/*If we have managed to grab an ATOM record, we shall parse it */
	if ((pdbline[0] == 'A') && (pdbline[1] == 'T') && (pdbline[2] == 'O') && (pdbline[3] == 'M')) 
	{  
	    sscanf(pdbline,"%*12c %s", &atom_type);
	    /*Grab the backbone atoms*/
	    if( (strcmp(atom_type,"N") == 0) || (strcmp(atom_type,"CA") == 0) || (strcmp(atom_type,"C") == 0) ||  (strcmp(atom_type,"O") == 0) )
	    {
		sscanf(pdbline,"%*30c %lf %lf %lf ", &b[i][0], &b[i][1], &b[i][2]);
		i++;
	    }
	}
    }

    fclose(templatefile);

if(DEBUG)
{
    fprintf(debugfile,"Reading 3D coordinates from the template file %s\n", templatefilename);
    fprintf(debugfile,"Printing matrix b\n");
    for(i=0; i<NUMPOINTS; i++)   
	fprintf(debugfile,"%f %f %f \n", b[i][0], b[i][1], b[i][2]);
    fprintf(debugfile,"\n");
}

    get_rmsd (a,b,c,w,&rmsd);

    /*Threshold this rmsd value against a threshold of 0.75 to determine the type. 
      If preliminary determination yields a Helix type, we need to check
      on more input data - the C array included*/
    if (rmsd < THRESHOLD)
    {
	switch (iter)   {
	  case 0:  strcpy(type, "EXTE");
		   return SUCCESS;
        	   break;
	  case 1:  /*patch a with the head of x and redetermine if its a helix*/
		   for(i=0; i<4; i++)
		   {
		       tmp[i][0] = x[i][0];
		       tmp[i][1] = x[i][1];
		       tmp[i][2] = x[i][2];
		   }
        	   for(j=0;j<NUMPOINTS-4;j++,i++) 
		   {
		       tmp[i][0] = a[j][0];
		       tmp[i][1] = a[j][1];
		       tmp[i][2] = a[j][2];
		   }
		   get_rmsd (tmp,b,c,w,&rmsd);
		   if (rmsd < THRESHOLD)
		   {
		       strcpy(type, "HELI");
		       return SUCCESS;
		   }
		   else 
		   {
		   /*patch a with the tail of x and redetermine if its a helix*/
	        	  for(i=0;i<(NUMPOINTS-4);i++)
			  {
			      tmp[i][0] = a[i+4][0];
			      tmp[i][1] = a[i+4][1];
			      tmp[i][2] = a[i+4][2];
			  }
			  for(k=4;i<NUMPOINTS;i++,k++)
			  {
			      tmp[i][0] = x[k][0];
			      tmp[i][1] = x[k][1];
			      tmp[i][2] = x[k][2];
			  }
			  get_rmsd (tmp,b,c,w,&rmsd);
			  if (rmsd < THRESHOLD)
			  {
			      strcpy(type, "HELI");
			      return SUCCESS;
			  }
		    }
	           break;
	  case 2:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;
	  case 3:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;
	  case 4:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;
	  case 5:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;
	  case 6:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;

	    }
    }
    iter++;
}/*End of while loop*/

if (minrmsd[1] != 0 && minrmsd[0] != THRESHOLD) /*If the initial values have been changed, we note that we have hit one of the types */
{
    switch((int)minrmsd[1])  {
      case 2:  strcpy(type, "TUR1");
	       return SUCCESS;
      case 3:  strcpy(type, "TUR2");
	       return SUCCESS;
      case 4:  strcpy(type, "TT1P");
	       return SUCCESS;
      case 5:  strcpy(type, "TT2P");
	       return SUCCESS;
      case 6:  strcpy(type, "TUR8");
	       return SUCCESS;
    }
}
else
{   
  int firstCAIndex = 1;
  int lastCAIndex = 13;
  float thres = 7.0;
  float x1,y1,z1,x2,y2,z2,dist;

   x1=a[firstCAIndex][0];   
   y1=a[firstCAIndex][1];   
   z1=a[firstCAIndex][2];   
   x2=a[lastCAIndex][0];   
   y2=a[lastCAIndex][1];   
   z2=a[lastCAIndex][2];   
  
  dist = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
  if(dist <= thres) 
    strcpy(type,"TIRR");
  else
    strcpy(type,"IRRE");

  return SUCCESS;
}

}/*end of function match_with_template*/

/*PCS- The following is a modifed version of the above match_with_template */
/*  to allow for the addition of comparing against a user defined template */
/*  and reporting the resulting RMSD as a reference parameter to the caller*/
/*  function. Note: This does NOT replace the normal str assignment, it is */
/*  in addition to this                                                    */
int match_with_template_and_user(double a[][DIM], double x[][DIM], char *type, double *UserTemplateRMSD)
{
double b[NUMPOINTS][DIM], c[NUMPOINTS][DIM], tmp[NUMPOINTS][DIM], w[NUMPOINTS];
double  rmsd, minrmsd[2];
int  i,j,k =0;
int iter = -1;
char templatefilename[FILENAMELEN], pdbline[PDBLINESIZE];
char atom_type[4];
FILE *templatefile;

/*Initializing the index and value of the min rmsd value to keep tab of the turn types */
minrmsd[0] = THRESHOLD;  /*we are looking for a value that is less than 0.75*/
minrmsd[1] = 0;   /*This check is done only for turn types, not Beta or helix*/

/*Initializing weights to a value of 1*/
for(i=0; i<NUMPOINTS; i++)
    w[i] = 1;

/*Checking the validity of the input data; The anamolies could be that we could not grab the 3D coordinates for all 4 atoms of the backbone chain. If there is an anamoly, we simply return with an error message */
for(i=0;i<NUMPOINTS;i++)
    if (a[i][0] == BOUND || a[i][1] == BOUND || a[i][2] == BOUND) 
    {
       printf("ERROR: Backbone atoms from PDB file could not be matched with the template for\
superposition, probably due to missing or extra atoms. - ");
       strcpy(type,"ERROR");
       fprintf(logfile,"ERROR: Backbone atoms from PDB file could not be matched with the template\
for superposition, probably due to missing or extra atoms. - ");
       return ERROR;
    }
for(i=0;i<NUMPOINTS/2;i++)
    if(x[i][0] == BOUND || x[i][1] == BOUND || x[i][2] == BOUND)
    {
      printf("ERROR: Insufficient Input Data. Couldnt grab residues on the lower and upper boundary of the sequence. - ");
      strcpy(type,"ERROR");
      fprintf(logfile,"ERROR: Insufficient Input Data. Couldnt grab residues on the lower and upper boundary of the sequence. - ");
      return ERROR;
    }
while ( iter < 7)      /* iter is the total num of types that we are matching against. If any of the matches succeeds, we return*/
  {
    switch(iter)  {
      /*----------Start of PCS----------*/
      /*Following is change in iteration. The first iteration is the user template */
      /*  comparison and is always performed.                                      */
      case -1: strcpy(templatefilename,userfilename);
	       break;
      /*-----------End of PCS-----------*/
      case 0: strcpy(templatefilename,betafilename);
              break;
      case 1: strcpy(templatefilename,helixfilename);
              break;
      case 2: strcpy(templatefilename,tt1filename);
              break;
      case 3: strcpy(templatefilename,tt2filename);
              break;
      case 4: strcpy(templatefilename,tt1pfilename);
              break;
      case 5: strcpy(templatefilename,tt2pfilename);
              break;
      case 6: strcpy(templatefilename,tt8filename);
              break;
    }

    /*Open the template file and Input the template residues' 3D coordinates*/
    if( (templatefile = fopen(templatefilename,"r")) == NULL)
    {
	printf("UNABLE to open template file %s - Returning without processing other sequery matches\n", templatefilename);
	fprintf(logfile,"UNABLE to open template file %s - Returning without processing other sequery matches\n", templatefilename);
	return ERROR;
    }
    
    /*Reading the 3D coordinates from the template file */
    i = 0;
    while (!(feof(templatefile)) && i<NUMPOINTS)
    {
	fgets(pdbline,PDBLINESIZE+1, templatefile);
	/*If we have managed to grab an ATOM record, we shall parse it */
	if ((pdbline[0] == 'A') && (pdbline[1] == 'T') && (pdbline[2] == 'O') && (pdbline[3] == 'M')) 
	{  
	    sscanf(pdbline,"%*12c %s", &atom_type);
	    /*Grab the backbone atoms*/
	    if( (strcmp(atom_type,"N") == 0) || (strcmp(atom_type,"CA") == 0) || (strcmp(atom_type,"C") == 0) ||  (strcmp(atom_type,"O") == 0) )
	    {
		sscanf(pdbline,"%*30c %lf %lf %lf ", &b[i][0], &b[i][1], &b[i][2]);
		i++;
	    }
	}
    }

    fclose(templatefile);
if(DEBUG)
{
    fprintf(debugfile,"Reading 3D coordinates from the template file %s\n", templatefilename);
    fprintf(debugfile,"Printing matrix b\n");
    for(i=0; i<NUMPOINTS; i++)   
	fprintf(debugfile,"%f %f %f \n", b[i][0], b[i][1], b[i][2]);
    fprintf(debugfile,"\n");
}

    get_rmsd (a,b,c,w,&rmsd);
    /*Threshold this rmsd value against a threshold of 0.75 to determine the type. If preliminary determination yields a Helix type, we need to check on more input data - the C array included*/
    /*----------Start of PCS----------*/
    /*If we are comparing to the user template, iter = -1, we set the resultant RMSD to */
    /*  the reference parameter to be sent back to the caller function                  */
    if (iter == -1)
      {
	*UserTemplateRMSD = rmsd;
      }
    if (rmsd < THRESHOLD && iter != -1)
    {
/*-----------End of PCS-----------*/
	switch (iter)   {
	  case 0:  strcpy(type, "EXTE");
		   return SUCCESS;
        	   break;
	  case 1:  /*patch a with the head of x and redetermine if its a helix*/
		   for(i=0; i<4; i++)
		   {
		       tmp[i][0] = x[i][0];
		       tmp[i][1] = x[i][1];
		       tmp[i][2] = x[i][2];
		   }
        	   for(j=0;j<NUMPOINTS-4;j++,i++) 
		   {
		       tmp[i][0] = a[j][0];
		       tmp[i][1] = a[j][1];
		       tmp[i][2] = a[j][2];
		   }
		   get_rmsd (tmp,b,c,w,&rmsd);
		   if (rmsd < THRESHOLD)
		   {
		       strcpy(type, "HELI");
		       return SUCCESS;
		   }
		   else 
		   {
		   /*patch a with the tail of x and redetermine if its a helix*/
	        	  for(i=0;i<(NUMPOINTS-4);i++)
			  {
			      tmp[i][0] = a[i+4][0];
			      tmp[i][1] = a[i+4][1];
			      tmp[i][2] = a[i+4][2];
			  }
			  for(k=4;i<NUMPOINTS;i++,k++)
			  {
			      tmp[i][0] = x[k][0];
			      tmp[i][1] = x[k][1];
			      tmp[i][2] = x[k][2];
			  }
			  get_rmsd (tmp,b,c,w,&rmsd);
			  if (rmsd < THRESHOLD)
			  {
			      strcpy(type, "HELI");
			      return SUCCESS;
			  }
		    }
	           break;
	  case 2:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;
	  case 3:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;
	  case 4:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;
	  case 5:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;
	  case 6:  if (rmsd < minrmsd[0]) /*If this is the min rmsd that we have seen so far, we shall store the rmsd value and the index */
	           {
	              minrmsd[1] = iter; 
		      minrmsd[0] = rmsd;
		   }
	           break;

	    }
    }
    iter++;
}/*End of while loop*/

if (minrmsd[1] != 0 && minrmsd[0] != THRESHOLD) /*If the initial values have been changed, we note that we have hit one of the types */
{
    switch((int)minrmsd[1])  {
      case 2:  strcpy(type, "TUR1");
	       return SUCCESS;
      case 3:  strcpy(type, "TUR2");
	       return SUCCESS;
      case 4:  strcpy(type, "TT1P");
	       return SUCCESS;
      case 5:  strcpy(type, "TT2P");
	       return SUCCESS;
      case 6:  strcpy(type, "TUR8");
	       return SUCCESS;
    }
}
else
{
  int firstCAIndex = 1;
  int lastCAIndex = 13;
  float thres = 7.0;
  float x1,y1,z1,x2,y2,z2,dist;

   x1=a[firstCAIndex][0];   
   y1=a[firstCAIndex][1];   
   z1=a[firstCAIndex][2];   
   x2=a[lastCAIndex][0];   
   y2=a[lastCAIndex][1];   
   z2=a[lastCAIndex][2];   
  
  dist = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
  if(dist <= thres) 
      strcpy(type,"TIRR");
  else
    strcpy(type, "IRRE");

  return SUCCESS;
}

}/*end of function match_with_template_and_user*/

void get_rmsd (double a[][DIM], double b[][DIM], double c[][DIM], double w[], double *rmsd)
{
double r[3][3], t[3];
double a_t[NUMPOINTS][DIM], d_sqr[NUMPOINTS], sum_d_sqr=0;
int i;

/*call the rotlsqfit function */
rotlsqfit(a,b,c,w,NUMPOINTS,r,t);

if(DEBUG)
{
    fprintf(debugfile,"Translation Vector is: %f %f %f\n", t[0], t[1], t[2]);
    fprintf(debugfile,"Printing Rotation Matrix\n");
    for(i=0;i<3;i++)
	fprintf(debugfile," %f %f %f\n", r[i][0], r[i][1], r[i][2]);
    fprintf(debugfile,"\n");
}

/*Calculating the RMSD */
/*Transformed points of matrix: r*a + t = a_t */
for(i=0; i<NUMPOINTS; i++)
{
    a_t[i][0] = r[0][0]*a[i][0] + r[0][1]*a[i][1] + r[0][2]*a[i][2] + t[0];
    a_t[i][1] = r[1][0]*a[i][0] + r[1][1]*a[i][1] + r[1][2]*a[i][2] + t[1];
    a_t[i][2] = r[2][0]*a[i][0] + r[2][1]*a[i][1] + r[2][2]*a[i][2] + t[2];
}

if(DEBUG)
{
    /*  Output the transformed matrix a_t*/
    fprintf(debugfile,"Printing transformed matrix a_t\n");
    for(i=0; i<NUMPOINTS; i++)   
	fprintf(debugfile,"%f %f %f \n", a_t[i][0], a_t[i][1], a_t[i][2]);
    fprintf(debugfile,"\n");
}

/*RMSD between the points in a_t and b*/
for(i=0; i<NUMPOINTS; i++)
    d_sqr[i] = ( pow((a_t[i][0] - b[i][0]),2) +  pow((a_t[i][1] - b[i][1]),2) + pow((a_t[i][2] - b[i][2]),2) );

for(i=0; i<NUMPOINTS; i++)
    sum_d_sqr += d_sqr[i];

*rmsd = sqrt(sum_d_sqr/NUMPOINTS);

if(DEBUG)
{
    /*Print the final answer - the rmsd value */
    fprintf(debugfile,"The RMSD value for a and b is: %f\n\n", *rmsd);
}

}
