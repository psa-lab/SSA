/**************************************************************************
SSA, Superpositional Structural Assignment, is a template-based
algorithm for assigning secondary structure to tetrapeptides from
Protein Data Bank files that was developed, coded, and validated by
Paul Sanschagrin, Tharuna Niranjan, and Leslie Kuhn at Michigan State
University.  This software should not to be modified, distributed, or
used for profit without prior agreement from Leslie Kuhn, who can be
contacted at the address below.

Copyright 1996, Michigan State University

Leslie Kuhn
Protein Structural Analysis and Design Laboratory
Department of Biochemistry
Michigan State University
East Lansing, MI 48824-1319

(517) 353-8745
(517) 353-9334 (fax)
kuhn@agua.bch.msu.edu (e-mail)

SSA should be attributed to the above authors and institution in any
publication or presentation deriving from its use, and a copy of such
publications should be sent to the authors.  This will not only
provide useful feedback on how the program is being used, but also
help us obtain support for future software development.  Suggestions
on SSA improvements are also welcome!

------------------------------------------------------------------------------

SSA assigns structure for Protein Data Bank tetrapeptides listed in a
Sequery file provided as input to SSA.  Sequery is a search tool for
identifying matches to user-specified sequence patterns in a sequence
file created from the Protein Data Bank.  Sequery was developed by
Leslie Kuhn, Michael Pique, Elizabeth Getzoff, and John Tainer at The
Scripps Research Institute, and may be obtained by contacting:

Michael Pique
Department of Molecular Biology /MB5
The Scripps Research Institute
10666 North Torrey Pines Road
La Jolla, CA  92037
(619) 554-9775
mp@scripps.edu (e-mail)

Sequery and some of its applications are described in the following papers:

J. F. Collawn, M. Stangel, L. A. Kuhn, V. Esekogwu, S. Jing, I. S.
Trowbridge, and J. A. Tainer (1990) "Transferrin Receptor
Internalization Sequence YXRF Implicates a Tight Turn as the Structural
Recognition Motif for Endocytosis",  Cell 63, 1061-1072.

J. F. Collawn, L. A. Kuhn, L.-F. S. Liu, J. A. Tainer, and I. S.
Trowbridge  (1991) "Transplanted LDL and Mannose-6-Phosphate Receptor
Internalization Signals Promote High-Efficiency Endocytosis of the
Transferrin Receptor", EMBO Journal 10, 3247-3253.
******************************************************************************/

/***************************************************************************
V1.0
   Tharuna Niranjan
   Kuhn Laboratory
   Dept. of Biochemistry
   Michigan State University

   ssa.c    Main file for the automated SSA procedure       6.13.96

This file parses the sequery file by each line and opens the relevant
pdbfile (process_pdbfile function) and grabs the 3-D coordinate points 
of the back bone chain atoms and stores them in a 16 x 3 matrix. 
This matrix is then input to the template matching function (match_
with_template), which matches the sequery sequence against different 
template types by doing a least squares fit and thresholding the RMSD
values thus obtained. SSA determines the type of the sequery match and
also provides statistics of all the matches in that sequery file.

INPUT: The sequery file with fullpath name.

       The program looks for the PDB file in the directory path as set
       in the PDBHOME environment variable. Make sure you have this set
       before you run the SSA program. 
       ** Do This: setenv $PDBHOME <pathname where you have the pdb files>

       The program will also ask you for the template filenames and the
       filename to save the results in. 

OUTPUT: if the  DEBUG flag is entered on the command line as an 
        optional argument, the SSA will print debug data (all the 
        calculated values) in a debug file "ssa.dbg" created in the
        working directory.

        The SSA will print the output into a file specified 
        by you in <SaveFileName>. 

Makefile provided to compile all files belonging to this module
******************************************************************************/
/******************************************************************************
V2.0

Paul C. Sanschagrin
Kuhn Laboratory
Michigan State Univeristy
East Lansning MI 48824

Modifications to the SSA code described above were made by PCS on 5-Dec-96
These modifications are noted in the code by my initials (PCS)

These modifications allow the user to input a template which will be compared
  against each of the Sequery matches. The RMSD between the match and the user
  template will be reported along with an indication if the RMSD is below the
  THRESHOLD. This template will NOT be used as a structural assignment type.

The program has been modified as follows:
  1) in interactive mode, prompts the user if they wish to use a user template 
       file and if so, the filename of the template
  2) in command line mode, accepts the user template filename, if included,
       as the last command line argument
  3) calls match_with_template_and_user (instead of match_with_template) which
       acts like the original except it reports back the RMSD between the user
       template file and the Sequery match
  4) outputs the RMSD and a '*' if the RMSD is below the THRESHOLD at the end
       of the output line, following the assigned secondary structure type
  5) output the total number of Sequery matches that had RMSDs less than the
       THRESHOLD relative to the user template and the percentage of the
       Sequery matches that fit this criteria

The runssa has also been modifed to allow for the inclusion of a user template
  and now can be run as follows:

runssa [DEBUG] <SequeryFilename> [UserTemplateFilename] <SaveFilename>
where:
  DEBUG indicates to run in Debug mode and is optional
  SequeryFilename is the filename of the input Sequery file and is required
  UserTemplateFilename is the filename of the user template and is optional
  SaveFilename is the filename which to write ouput in and is required

NOTE: The UserTemplateFilename must be a PDB formatted file, but must 
  consist of only the four residues which are to be compared to the
  Sequery matches
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ssa.h"

#define SEQLINESIZE 81
#define TYPE_STR_LEN 6
/*struct _sequerydata{
    char pdb_code[5];
    char chain_id;
    int from_res_no, to_res_no;
};
*/

int process_pdbfile(struct _sequerydata *data, double a[NUMPOINTS][DIM], double c[NUMPOINTS/2][DIM]);

/*Global variables */
FILE * debugfile;
FILE * logfile;

char betafilename[FILENAMELEN], helixfilename[FILENAMELEN], tt1filename[FILENAMELEN], tt2filename[FILENAMELEN], tt1pfilename[FILENAMELEN], tt2pfilename[FILENAMELEN], tt8filename[FILENAMELEN] ;

int DEBUG = 0;
/*----------Start of PCS----------*/
int INCLUDEUSERTEMPLATE = 0;
char userfilename[FILENAMELEN];
/*-----------End of PCS-----------*/


main(int argc, char **argv) {
char sequeryfilename[FILENAMELEN];
char savefilename[FILENAMELEN];
/*----------Start of PCS----------*/
char UserTemplateResponse[1], UserTemplateRMSDString[4];
double UserTemplateRMSD; 
float percent_tot_usermatches;
int tot_usermatches = 0;
/*-----------End of PCS-----------*/
struct _sequerydata seq_data;
char seqline[SEQLINESIZE], type[TYPE_STR_LEN], outputline[SEQLINESIZE];
char x, junk[3], sequence[25];
int i, k, sign, errcode, cnt, len;
FILE *seqfile, *savefile;
double a[NUMPOINTS][DIM], b[NUMPOINTS][DIM], c[NUMPOINTS/2][DIM], rmsd;
int tot_seqlines = 0, tot_beta = 0, tot_helix = 0, tot_turn1 = 0, tot_turn2 = 0, tot_turn1p = 0, tot_turn2p = 0, tot_turn8 = 0, tot_irreg = 0, tot_errs = 0;
int tot_turns; 
int tot_irreg_turns = 0;/*Added by Vishal */
int tot_irreg_helix = 0;/*Added by Vishal */
float percent_tot_turns, percent_tot_turn1, percent_tot_turn2, percent_tot_turn1p,
  percent_tot_turn2p, percent_tot_turn8, percent_tot_helix, 
  percent_tot_beta, percent_tot_irreg, percent_tot_errs,
  percent_tot_irreg_turns,percent_tot_irreg_helix; /*last line by Vishal */

if (argc > 1 && (!strcasecmp(argv[1],"DEBUG")) )
    DEBUG = 1;

if(DEBUG) {
  if( (debugfile = fopen("ssa.dbg","w")) == NULL) {
    printf("ERROR! Unable to open the debug file ssa.dbg\n");
    exit(1);
  }
}
/*if the user chose to run SSA interactively from the command line*/
if(argc == 1) {
    printf("Please Input the Filename of the Sequery File\n");
    scanf("%s", &sequeryfilename);
    
    /*Get the filenames of all the template files*/
    printf("Please specify the name of the Beta Template file.Add path if not in curr dir\n");
    scanf("%s", betafilename);
    printf("Please specify the name of the Helix Template file. Add path if not in curr dir\n");
              scanf("%s", helixfilename);
    printf("Please specify the name of the Turn Type 1 Template file. Add path if not in curr dir\n");
    scanf("%s", tt1filename);
    printf("Please specify the name of the Turn Type 2 Template file. Add path if not in curr dir\n");
    scanf("%s", tt2filename);
    printf("Please specify the name of the Turn Type 1P Template file. Add path if not in curr dir\n");
    scanf("%s", tt1pfilename);
    printf("Please specify the name of the Turn Type 2P Template file. Add path if not in curr dir\n");
              scanf("%s", tt2pfilename);
    printf("Please specify the name of the Turn Type 8 Template file. Add path if not in curr dir\n");
    scanf("%s", tt8filename);
    /*----------Start of PCS----------*/ 
    /*Determine if the user wishes to compare the Sequery matches to a user template*/
    printf("Would you like to compare the Seqeury matches in the Sequery file to a user \n");
    printf("defined template? (Y or N) \n");
    printf("(Note: This template will not be used to classify sequery structure matches, but the \n");
    printf("RMSD between it and the Seqeury along with an indication if it is below the threshold.\n");
    scanf("%s", UserTemplateResponse);
    if ((UserTemplateResponse[0] == 'Y') || (UserTemplateResponse[0] == 'y')) {
      printf("Please specify the name of the User Defined Template file. Add path if not in curr dir\n");
      scanf("%s", userfilename); 
      INCLUDEUSERTEMPLATE = 1;
    }
    /*-----------End of PCS-----------*/
    /*Get the save filename*/
    printf("Please specify the name of the file to save the output in. Add path if not in curr dir\n");
    scanf("%s", savefilename);
}
else
/*If user enters insufficient number of command line arguments*/
if(argc > 1 && argc < 10) {
  /*PCS-Added in optional DEBUG and optional UserTemplateFilename*/
    printf("SSA USAGE: ssa [DEBUG] InputFilename BetaTemplateFilename HelixTemplateFilename Type1-turnTemplateFilename ");
    printf("Type2-turnTemplateFilename Type1'-turnTemplateFilename Type2'-turnTemplateFilename Type8-turnTemplateFilename");
    printf("OutputFilename [UserTemplateFilename]\n\n");
    exit (1);
}

else
/*if user enters the correct number of parameters and does not want to run SSA interactively, we shall grab the filenames from a scriptfile*/
/*PCS - 10 arguments indicates no usertemplate file desired and not in debug mode*/
if(argc == 10) {
strcpy(sequeryfilename,argv[1]);
strcpy(betafilename,argv[2]);
strcpy(helixfilename,argv[3]);
strcpy(tt1filename,argv[4]);
strcpy(tt2filename,argv[5]);
strcpy(tt1pfilename,argv[6]);
strcpy(tt2pfilename,argv[7]);
strcpy(tt8filename,argv[8]);
strcpy(savefilename,argv[9]);
}
else
/*----------Start of PCS----------*/
/*if the user enters 11 arguments, it could indicate either to run with standard templates */
/*  in debug mode or to run with the user defined template file, but not in debug mode     */
/*If the former, DEBUG will have benn set previously and the remainder will be assinged as */
/*  in the 12 argument case, without assignment of the userfilename                        */
if (argc == 11 && DEBUG == 1) {
  strcpy(sequeryfilename,argv[2]);
  strcpy(betafilename,argv[3]); 
  strcpy(helixfilename,argv[4]);
  strcpy(tt1filename,argv[5]);  
  strcpy(tt2filename,argv[6]);  
  strcpy(tt1pfilename,argv[7]); 
  strcpy(tt2pfilename,argv[8]); 
  strcpy(tt8filename,argv[9]);  
  strcpy(savefilename,argv[10]);
}
else
/*If the latter, i.e. user wants non-debug mode with a user template, DEBUG will not be set, */
/*  i.e. it will = 0, and the assignments will be done as in the 10 argument case with the   */
/*  addition of the user template file.                                                      */
if (argc == 11 && DEBUG == 0) {
  strcpy(sequeryfilename,argv[1]);
  strcpy(betafilename,argv[2]);
  strcpy(helixfilename,argv[3]);
  strcpy(tt1filename,argv[4]); 
  strcpy(tt2filename,argv[5]); 
  strcpy(tt1pfilename,argv[6]);
  strcpy(tt2pfilename,argv[7]);
  strcpy(tt8filename,argv[8]); 
  strcpy(savefilename,argv[9]);
  strcpy(userfilename,argv[10]);
  INCLUDEUSERTEMPLATE = 1;
}
else
/*-----------End of PCS-----------*/

/*if user enters the correct number of parameters,desires to run in debug mode, but does not want to run SSA interactively, we shall grab the filenames from a scriptfile*/
/*PCS - Changed below to 12 -- This is max # - indicates debug mode and user template*/

if(argc == 12) {
strcpy(sequeryfilename,argv[2]);
strcpy(betafilename,argv[3]);
strcpy(helixfilename,argv[4]);
strcpy(tt1filename,argv[5]);
strcpy(tt2filename,argv[6]);
strcpy(tt1pfilename,argv[7]);
strcpy(tt2pfilename,argv[8]);
strcpy(tt8filename,argv[9]);
strcpy(savefilename,argv[10]);
strcpy(userfilename,argv[11]);
/*----------Start of PCS----------*/
INCLUDEUSERTEMPLATE = 1;
/*----------Start of PCS----------*/
}
 
printf("\n\n\n"); 

/*Parse the sequery file and grab the Pdb code, chain_id and range of residue numbers */
if( (seqfile = fopen(sequeryfilename, "r")) == NULL) {
    printf("UNABLE to open sequery file %s\n", sequeryfilename);
    exit (1);
}

if( (savefile = fopen(savefilename, "w")) == NULL) {
    printf("UNABLE to open save file %s\n", savefilename);
    exit (1);
}

if( (logfile = fopen("ssa.log", "w")) == NULL) {
    printf("UNABLE to open log file \"ssa.log\"\n");
    exit (1);
}

/*Print the sequery filename in the output file so that we can relate both files*/
/*PCS- Added "Sequery Filename" to below */
fprintf(savefile,"# Sequery Filename : %s\n", sequeryfilename);

/*PCS-Added "Sequery Filename" to below*/
/*Print the sequery filename in the log file so that we can relate both files*/

fprintf(logfile,"# Sequery Filename : %s\n", sequeryfilename);

/*----------Start of PCS----------*/
if (INCLUDEUSERTEMPLATE) {
  /* Output the user template filename both the output file and the log file */
  fprintf(savefile,"# User Template : %s\n", userfilename);
  fprintf(logfile,"# User Template : %s\n", userfilename);
}  
/*-----------End of PCS-----------*/

fprintf(savefile,"\n");
fprintf(logfile,"\n");

if(DEBUG) {
  /*Print the sequery filename in the debug file so that we can relate both files*/
  /*PCS-Added "Sequery Filename" to below*/
  fprintf(debugfile,"# Sequery Filename : %s\n", sequeryfilename);
  /*----------Start of PCS----------*/
  if (INCLUDEUSERTEMPLATE)   
    {
      /* If including a user template file, output the filename to the debug file */
      fprintf(debugfile,"# User Template : %s\n", userfilename);
    }
  fprintf(debugfile,"\n");
  /*-----------End of PCS-----------*/

}

/*----------Start of PCS----------*/
/*Following prints header lines into output file */
fprintf(savefile,"#    C\n#    h\n#    a\n");
fprintf(savefile,"#    i   Residue                                         User User\n");
fprintf(savefile,"#PDB n      Range      Sequence                    Type  RMSD Match?\n");
fprintf(savefile,"#------------------------------------------------------------------\n");
/*-----------End of PCS-----------*/

/*Do For each line in the sequery file */
while (!(feof(seqfile))) {
    /*Read a line from the sequery file */
    fgets(seqline, SEQLINESIZE+1, seqfile);
    fscanf(seqfile,"\n");
    if((k=strncmp(seqline,"#",1)) != 0)  {
        /*count the number of seq lines of sequences read */
        tot_seqlines++;


        sscanf(seqline,"%s %c %d %*s %d %*s %s", &seq_data.pdb_code,
               &seq_data.chain_id, 
               &seq_data.from_res_no, 
               &seq_data.to_res_no, 
               &sequence);


        /*To determine the no of sequery residues by parsing the sequence*/
        cnt = 0;
        for(i=0; sequence[i] != '\0'; i++) {
            if (!(isupper(sequence[i])));
            else cnt++;
        }
        
        /*Error checking. The program will not process the pdb file for these sequery residues*/
        if(cnt < 4) {
            printf("Error in Sequery matching (Based on sequence \"%s\" ) for PDB code \"%s\".  \
                Cannot match less than four residues with tetrapeptide template\n",sequence,seq_data.pdb_code); 
	    printf("SeqLine: ");
	    puts(seqline);
	    tot_errs++;
            fprintf(logfile,"Error in Sequery matching (Based on sequence \"%s\" ) for PDB code \"%s\". \
               Cannot match less than four residues with tetrapeptide template\n",sequence,seq_data.pdb_code);
            continue;	    
        }
        /*Error checking, but the .dat file will be created. Upto the user whether to process the next step of
	  lsqfit and template matching.*/
        else
	  if(cnt > 4) {
	  printf("WARNING! More than four residues (Based on sequence \"%s\" ) for PDB code \"%s\". \
                 Matching the first four residues with tetrapeptide template\n",sequence,seq_data.pdb_code); 
	  printf("SeqLine: ");
	  puts(seqline);
	  tot_errs++;
	  fprintf(logfile,"WARNING! More than four residues (Based on sequence \"%s\" ) for PDB code \"%s\".\
                 Matching the first four residues with tetrapeptide template\n",sequence,seq_data.pdb_code);             
        }
	
        /*debug ; outputtting information received */
	if(DEBUG) {
	  fprintf(debugfile,"%s\t",seq_data.pdb_code);
	  fprintf(debugfile,"%c\t", seq_data.chain_id);
	  fprintf(debugfile,"%d\t%d\t", seq_data.from_res_no,seq_data.to_res_no);
	  fprintf(debugfile,"%s\n\n", sequence);
	}

        /*Initialize the arrays a and c before filling them in the next routine*/
        for(i=0; i<NUMPOINTS; i++)   
	  a[i][0] = a[i][1] = a[i][2] = BOUND; /*assuming that the actual values are not equal BOUND (defined in ssa.h) */
        for(i=0; i<NUMPOINTS/2; i++)   
	  c[i][0] = c[i][1] = c[i][2] = BOUND;
            
        /*Now that we have seq_data, open the relevant pdbfile and grab the pdb record - what is of primary 
	  importance to us is the x,y,z coordinates */
        if( (errcode = process_pdbfile(&seq_data,a,c)) != SUCCESS ) {	  
            for(i=0;i<SEQLINESIZE;i++)
                outputline[i] = '\0';
            len = strlen(seqline);
            if (seqline[len-1] == '\n')
                seqline[len-1] = '\0';
            strncpy(outputline,seqline,72);
            strcat(outputline,"  ");
            strcat(outputline, "ERROR");
            fprintf(savefile,"%s\n", outputline);
            /*Increment the error count*/
            tot_errs++;
            continue;
        }

	if(DEBUG) {
	  fprintf(debugfile,"Printing matrix a\n");
	  for(i=0; i<NUMPOINTS; i++)   
	    fprintf(debugfile,"%f %f %f \n", a[i][0], a[i][1], a[i][2]);
	  fprintf(debugfile,"\n");
	  
	  fprintf(debugfile,"Printing matrix c\n");
	  for(i=0; i<NUMPOINTS/2; i++)  
	    fprintf(debugfile,"%f %f %f \n", c[i][0], c[i][1], c[i][2]);
	  fprintf(debugfile,"\n");
	}
	
    }
    /*PCS-Following block is done of not using user template and is left intact as */
    /*  from the original code except for the initial if statement                 */
    if (!INCLUDEUSERTEMPLATE){
        /*Now that we have the data, we match the data with our templates for Beta, Helix, 
	  and Turns to determine the type of the data ; Template matching is done by least squares 
	  fit and rmsd thresholding. Given the data, this function will determine the type. */
     
      if( (errcode = match_with_template(a, c, type)) != SUCCESS) {
	
	for(i=0;i<SEQLINESIZE;i++)
	  outputline[i] = '\0';
	len = strlen(seqline);
	if (seqline[len-1] == '\n')
	  seqline[len-1] = '\0';

	  strncpy(outputline,seqline,72);
	  strcat(outputline,"  ");
	  strcat(outputline, "ERROR");

	  fprintf(savefile,"%s\n", outputline);
	  /*Increment the error count*/
	  tot_errs++;
	  printf("PDB code is \"%s\" and sequence is \"%s\" \n",seq_data.pdb_code,sequence);
	  fprintf(logfile,"PDB code is \"%s\" and sequence is \"%s\" \n",seq_data.pdb_code,sequence);
      }
      else {
	/*starting by putting 2 spaces (as a separator), I will append the type result at the end
	  of the sequery line. If the sequery line is long  (i.e. more than 80 char long), I will
	  over write the last 5 places with my result. */
            for(i=0;i<SEQLINESIZE;i++)
                outputline[i] = '\0';
            len = strlen(seqline);
            if (seqline[len-1] == '\n')
                seqline[len-1] = '\0';

            if(!strcmp(type,"EXTE"))
                tot_beta++;
            else
            if(!strcmp(type,"HELI"))
                tot_helix++;
            else
            if(!strcmp(type,"TUR1"))
                tot_turn1++;
            else
            if(!strcmp(type,"TUR2"))
                tot_turn2++;
            else
            if(!strcmp(type,"TT1P"))
                tot_turn1p++;
            else
            if(!strcmp(type,"TT2P"))
                tot_turn2p++;
            else
            if(!strcmp(type,"TUR8"))
                tot_turn8++;
            else
            if(!strcmp(type,"IRRE"))
                tot_irreg++;        
            else
            if(!strcmp(type,"TIRR")) {
	      switch (check_lr_helix(&seq_data,type,match_with_template)) {
	      case 1: tot_irreg_helix++;break;
	      case 0: tot_irreg_turns++;break;
	      case -1: tot_errs++;break;
	      }        
	    }
	    if(DEBUG) {
	      fprintf(debugfile,"The type as determined for this sequery match is: %s\n\n ", type);
	    }
	    
	    /* write type in output file */
	    strncpy(outputline,seqline,73);
	    strcat(outputline,"  ");
	    strcat(outputline,type);
	    fprintf(savefile,"%s\n", outputline);
      }
    }
	
    /*PCS-End of original block                                                  */
    /*The following is added code and is executed of the user wishes to include  */
    /*  a user template rmsd data in the output. This is modified from the above */
    /*  block of code                                                            */
    
    else if (INCLUDEUSERTEMPLATE) {
      /*Now that we have the data, we match the data with our templates for Beta, Helix, and Turns
	to determine the type of the data ; Template matching is done by least squares fit and rmsd
	thresholding. Given the data, this function will determine the type. */
      if( (errcode = match_with_template_and_user(a, c, type, &UserTemplateRMSD)) != SUCCESS) {               
	for(i=0;i<SEQLINESIZE;i++)
	  outputline[i] = '\0';
	len = strlen(seqline);
	if (seqline[len-1] == '\n')
	  seqline[len-1] = '\0';
	strncpy(outputline,seqline,72);
	strcat(outputline,"  ");
	strcat(outputline, "ERROR");
	fprintf(savefile,"%s\n", outputline);
	/*Increment the error count*/
	tot_errs++;
	printf("PDB code is \"%s\" and sequence is \"%s\" \n",seq_data.pdb_code,sequence);
	fprintf(logfile,"PDB code is \"%s\" and sequence is \"%s\" \n",seq_data.pdb_code,sequence);
      }
      else {
	/*starting by putting 2 spaces (as a separator), I will append the type result at the end
	  of the sequery line. If the sequery line is long  (i.e. more than 80 char long), I will 
	  over write the last 5 places with my result. */
	for(i=0;i<SEQLINESIZE;i++)
	  outputline[i] = '\0';
	len = strlen(seqline);
	if (seqline[len-1] == '\n')
	  seqline[len-1] = '\0';
	
            if(!strcmp(type,"EXTE"))
                tot_beta++;
            else
            if(!strcmp(type,"HELI"))
                tot_helix++;
            else
            if(!strcmp(type,"TUR1"))
                tot_turn1++;
            else
            if(!strcmp(type,"TUR2"))
                tot_turn2++;
            else
            if(!strcmp(type,"TT1P"))
                tot_turn1p++;
            else
            if(!strcmp(type,"TT2P"))
                tot_turn2p++;
            else
            if(!strcmp(type,"TUR8"))
                tot_turn8++;
            else
            if(!strcmp(type,"IRRE"))
                tot_irreg++;
            else
	    if(!strcmp(type,"TIRR")) {
	      switch (check_lr_helix(&seq_data,type,match_with_template)) {
	      case 1: tot_irreg_helix++;break;
	      case 0: tot_irreg_turns++;break;
	      case -1: tot_errs++;break;
	      }        
	      if(DEBUG) {
                fprintf(debugfile,"The type as determined for this sequery match is: %s\n\n ", type);
	      }
	    }

            /* write to file */
            strncpy(outputline,seqline,67);
            strcat(outputline,"  ");
            strcat(outputline,type);
            /*----------Start of PCS----------*/
            sprintf(UserTemplateRMSDString, " %1.2f", UserTemplateRMSD);
            strcat(outputline," ");
            strcat(outputline,UserTemplateRMSDString);
            strcat(outputline," ");
            /*Addtion of User Template match indicator*/
            if (UserTemplateRMSD < THRESHOLD) {
                tot_usermatches ++;
                strcat(outputline, "*");
	    }
            else
              strcat(outputline, " ");
            /*-----------End of PCS-----------*/
            fprintf(savefile,"%s\n", outputline);

      }
    }
    /*PCS-End of additional modified code                                                */
}

/*Compute statistics for this sequery file*/
tot_seqlines -= tot_errs;
if(tot_seqlines>0){
  tot_turns = tot_turn1 + tot_turn2 + tot_turn1p + tot_turn2p + tot_turn8 +tot_irreg_turns;
  tot_helix += tot_irreg_helix;
  percent_tot_turns = tot_turns * 100.0/tot_seqlines;
  percent_tot_turn1 = tot_turn1 * 100.0/tot_seqlines;
  percent_tot_turn2 = tot_turn2 * 100.0/tot_seqlines;
  percent_tot_turn1p = tot_turn1p * 100.0/tot_seqlines;
  percent_tot_turn2p = tot_turn2p * 100.0/tot_seqlines;
  percent_tot_turn8 = tot_turn8 * 100.0/tot_seqlines;
  percent_tot_beta = tot_beta * 100.0/tot_seqlines;
  percent_tot_helix = tot_helix * 100.0/tot_seqlines;
  percent_tot_irreg = tot_irreg * 100.0/tot_seqlines;
  percent_tot_irreg_turns = tot_irreg_turns*100.0/tot_seqlines;
  percent_tot_irreg_helix = tot_irreg_helix*100.0/tot_seqlines;
  percent_tot_errs = tot_errs * 100.0/(tot_seqlines+tot_errs);
  /*----------Start of PCS----------*/
  if (INCLUDEUSERTEMPLATE)
    {
      /*Calculate pecentage of the sequery mathces which matched the user template*/
      percent_tot_usermatches = tot_usermatches * 100.0/tot_seqlines;
    }
/*-----------End of PCS-----------*/
}
else{
  printf("All sequences failed.\n");
  exit(-1);
}
    

/*Now write the statistics to the output file */
/*PCS-I Slightly altered some of the below to attempt to improve output appearence*/
fprintf(savefile,"\n\n####################################STATISTICS##############################\n");
fprintf(savefile,"# Total number of sequery lines processed is:           %d\n", tot_seqlines+tot_errs);
fprintf(savefile,"# Total number of Errors (excluded from statistics):    %d  (%.2f %%)\n",tot_errs,percent_tot_errs);
fprintf(savefile,"# # of Regular Helix Types: %d\t (%.2f %%)\n",tot_helix-tot_irreg_helix, 
	                                                                             (tot_helix-tot_irreg_helix)*100.0/tot_seqlines);
fprintf(savefile,"# # of Irreg Helix Types:   %d\t (%.2f %%)\n",tot_irreg_helix, percent_tot_irreg_helix);
fprintf(savefile,"# Total # of Helix Types:   \t\t%d\t (%.2f %%)\n#\n", tot_helix, percent_tot_helix);

fprintf(savefile,"# # of Extended Types:      \t\t%d\t (%.2f %%)\n#\n", tot_beta, percent_tot_beta);
fprintf(savefile,"# # of Type1 Turns:         %d\t (%.2f %%)\n", tot_turn1, percent_tot_turn1);
fprintf(savefile,"# # of Type2 Turns:         %d\t (%.2f %%)\n", tot_turn2, percent_tot_turn2);
fprintf(savefile,"# # of Type1' Turns:        %d\t (%.2f %%)\n", tot_turn1p, percent_tot_turn1p);
fprintf(savefile,"# # of Type2' Turns:        %d\t (%.2f %%)\n", tot_turn2p, percent_tot_turn2p);
fprintf(savefile,"# # of Type8 Turns:         %d\t (%.2f %%)\n", tot_turn8, percent_tot_turn8);
fprintf(savefile,"# # of Irreg Turns:         %d\t (%.2f %%)\n",tot_irreg_turns, percent_tot_irreg_turns);
fprintf(savefile,"# Total # Of Turns:         \t\t%d\t (%.2f %%)\n#\n", tot_turns, percent_tot_turns);
fprintf(savefile,"# # of Irreg Types *:        \t\t%d\t (%.2f %%)\n", tot_irreg, percent_tot_irreg);
fprintf(savefile,"#\n#* # of Irreg Types excludes Irreg Helices and Irreg Turns\n");

/*----------Start of PCS----------*/
if (INCLUDEUSERTEMPLATE)
{
  fprintf(savefile,"# Number of User Matches:  %d\t and %% that matched is:  %.2f\n", tot_usermatches, percent_tot_usermatches);
}
/*-----------End of PCS-----------*/



fclose(seqfile);
fclose(savefile);
fclose(logfile);
if(DEBUG)
   fclose(debugfile);

}/*End main loop */

int  process_pdbfile(struct _sequerydata *data, double a[][DIM], double c[][DIM]){
char pdbfilename[FILENAMELEN];
char pdbline[PDBLINESIZE];
FILE *pdbfile;
int residue_number, prev_res_no, prev_residue_number; 
int i,j, k, n_i = 0, ca_i = 1, c_i = 2, o_i = 3;
int resflag, in_range;
char atom_type[4], *pdbhome;
 
/*Initializing variables to keep track of no. of residues grabbed; need to grab exactly 4*/
resflag = 0;
prev_res_no = 0;
prev_residue_number = 0;
/*The path name for the PDBfile should be read as an environment variable*/
if( (pdbhome = getenv("PDBHOME")) == NULL)
{
    printf("UNABLE to find the path for the PDBfile. Please check to see if you have set the environment variable\
           PDBHOME to the path of the directory that contains the PDBfiles. Returning from the function call\n");
    fprintf(logfile,"UNABLE to find the path for the PDBfile. Please check to see if you have set the\
           environment variable PDBHOME to the path of the directory that contains the PDBfiles. Returning from the function call\n");
    return ERROR;
}

strcpy(pdbfilename,pdbhome);
strcat(pdbfilename,"/pdb");
strcat(pdbfilename,data->pdb_code);
strcat(pdbfilename,".");
strcat(pdbfilename,"ent");

/*Open the pdbfile and read the relevant records */
if( (pdbfile = fopen(pdbfilename,"r")) == NULL)
{
    printf("UNABLE to open pdb file %s - Continuing to process other sequery matches\n", pdbfilename);
    fprintf(logfile,"UNABLE to open pdb file %s - Continuing to process other sequery matches\n", pdbfilename);
    return ERROR;
}
/*flag to indicate whether we are in the range of 4 residues specified as a sequence*/
in_range = 0;  /*0 means not in range, which is true initially, anyway*/
/*Read the pdb record and store it in the datafile*/
while(!(feof(pdbfile)))
{
    fgets(pdbline, PDBLINESIZE+1, pdbfile);
    /*If we have managed to grab an ATOM record, we shall parse it */
    if ((pdbline[0] == 'A') && (pdbline[1] == 'T') && (pdbline[2] == 'O') && (pdbline[3] == 'M')) 
    {  
        sscanf(pdbline,"%*22c %4d", &residue_number);
        if ( pdbline[21] == data->chain_id || pdbline[21] == ' ' ) /*then grab the backbone atoms*/
        {   
            /*Not in range; Actually Below range*/
            if(!(in_range))
            {
                if(residue_number == data->from_res_no)
                    in_range = 1;
                else
                {
                    /*store the previous residue number encountered */
                    if(prev_residue_number != residue_number)
                    {
                        prev_residue_number = residue_number;
                    }
                    sscanf(pdbline,"%*12c %s", &atom_type);
                    /*store previous info in the first half of matrix c*/
                    /*Storing the N atoms in order */
                    if (strcmp(atom_type,"N") == 0) 
                        sscanf(pdbline, "%*30c %lf %lf %lf ", &c[n_i][0], &c[n_i][1], &c[n_i][2]);
                   /*storing the CA atoms in order */ 
                    else                
                        if (strcmp(atom_type,"CA") == 0)
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[ca_i][0], &c[ca_i][1], &c[ca_i][2]);
                    /*Storing the C atoms in order */
                    else
                        if (strcmp(atom_type,"C") == 0)   
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[c_i][0], &c[c_i][1], &c[c_i][2]);
                    /*storing the O atoms in order */
                    else
                        if (strcmp(atom_type,"O") == 0)    
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[o_i][0], &c[o_i][1], &c[o_i][2]);
                 }
            }
            /*Within range of the tetrapeptide sequence*/
            if(in_range)
              {  /*Narrowing down to the four residues that we are interested in*/              
                if (residue_number >=(data->from_res_no) && residue_number <=(data->to_res_no) && resflag <= 5) 
                  {
                    sscanf(pdbline,"%*12c %s", &atom_type);

                    if (prev_res_no != residue_number)
                    {
                        resflag++;
                        prev_res_no = residue_number;
                    }

                    /*Storing the N atoms in order */
                    /*4 is the stride because we are considering 4 atoms in the backbone chain */
                    if (strcmp(atom_type,"N") == 0) 
                    {   
                        if (resflag == 5)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[n_i+4][0], &c[n_i+4][1], &c[n_i+4][2]); 
                        else
                        if (resflag == 1)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[n_i][0], &a[n_i][1], &a[n_i][2]);  
                        else
                        if (resflag == 2)     
                                 sscanf(pdbline, "%*30c %lf %lf %lf ", &a[n_i+4][0], &a[n_i+4][1], &a[n_i+4][2]);
                        else
                        if (resflag == 3)     
                                 sscanf(pdbline, "%*30c %lf %lf %lf ", &a[n_i+8][0], &a[n_i+8][1], &a[n_i+8][2]);
                        else
                        if (resflag == 4)     
                                 sscanf(pdbline, "%*30c %lf %lf %lf ", &a[n_i+12][0], &a[n_i+12][1], &a[n_i+12][2]);    
                    }
                    /*storing the CA atoms in order*/
                    else
                    if (strcmp(atom_type,"CA") == 0) 
                    {
                        if (resflag == 5)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[ca_i+4][0], &c[ca_i+4][1], &c[ca_i+4][2]); 
                        else
                        if (resflag == 1)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[ca_i][0], &a[ca_i][1], &a[ca_i][2]);       
                        else
                        if (resflag == 2)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[ca_i+4][0], &a[ca_i+4][1], &a[ca_i+4][2]);
                        else
                        if (resflag == 3)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[ca_i+8][0], &a[ca_i+8][1], &a[ca_i+8][2]);
                        else
                        if (resflag == 4)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[ca_i+12][0], &a[ca_i+12][1], &a[ca_i+12][2]);      
                        }
                    /*Storing the C atoms in order */
                    else
                    if (strcmp(atom_type,"C") == 0) 
                    {
                        if (resflag == 5)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[c_i+4][0], &c[c_i+4][1], &c[c_i+4][2]);  
                        else
                        if (resflag == 1)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[c_i][0], &a[c_i][1], &a[c_i][2]);  
                        else
                        if (resflag == 2)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[c_i+4][0], &a[c_i+4][1], &a[c_i+4][2]);
                        else
                        if (resflag == 3)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[c_i+8][0], &a[c_i+8][1], &a[c_i+8][2]);
                        else
                        if (resflag == 4)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[c_i+12][0], &a[c_i+12][1], &a[c_i+12][2]); 
                        }
                    /*storing the O atoms in order */
                    else
                    if (strcmp(atom_type,"O") == 0) 
                    {
                        if (resflag == 5)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[o_i+4][0], &c[o_i+4][1], &c[o_i+4][2]);    
                        else
                        if (resflag == 1)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[o_i][0], &a[o_i][1], &a[o_i][2]);  
                        else
                        if (resflag == 2)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[o_i+4][0], &a[o_i+4][1], &a[o_i+4][2]);
                        else
                        if (resflag == 3)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[o_i+8][0], &a[o_i+8][1], &a[o_i+8][2]);
                        else
                        if (resflag == 4)     
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &a[o_i+12][0], &a[o_i+12][1], &a[o_i+12][2]); 
                        }

                }
                else
                {
                    if((residue_number > data->to_res_no) && (resflag >= 4))
                    {
                        sscanf(pdbline,"%*12c %s", &atom_type);

                        if (prev_res_no != residue_number)
                        {
                             resflag++;

                            if(resflag >= 6)
                                break; /*break out of the while loop;this pdbfile has been processed*/                     
                            prev_res_no = residue_number;                      
                        }

                        /*store upperbound info in the secondhalf of matrix c*/
                        /*Storing the N atoms in order */
                        if (strcmp(atom_type,"N") == 0) 
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[n_i+4][0], &c[n_i+4][1], &c[n_i+4][2]);
                        /*storing the CA atoms in order */ 
                        else            
                        if (strcmp(atom_type,"CA") == 0)
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[ca_i+4][0], &c[ca_i+4][1], &c[ca_i+4][2]);
                        /*Storing the C atoms in order */
                        else
                        if (strcmp(atom_type,"C") == 0)   
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[c_i+4][0], &c[c_i+4][1], &c[c_i+4][2]);
                        /*storing the O atoms in order */
                        else
                        if (strcmp(atom_type,"O") == 0)    
                            sscanf(pdbline, "%*30c %lf %lf %lf ", &c[o_i+4][0], &c[o_i+4][1], &c[o_i+4][2]);
                    }
                }
            }/*Within range of tetrapeptide sequence*/
        }/*end of chain_id if loop*/
   }
 }
fclose(pdbfile);

return SUCCESS;
}

/* Test to see whether it is HIRR */
/* Algo: we will modify seq_data, read arrays, a and c,
 * and restore seq_data to original state 
 * shift left */
/* returns 1 for a HIRR assignment
 * returns -1 for an error (ERROR)
 * returns 0 for a TIRR assignment
 */
int check_lr_helix(struct _sequerydata *shifted, char *type, int (*func)(double [][DIM],double [][DIM],char *)) {
  
  double atmp[NUMPOINTS][DIM], 
    btmp[NUMPOINTS][DIM], 
    ctmp[NUMPOINTS/2][DIM], rmsd;               
  char new_type[10];
  int lefthelix =0, righthelix=0;
  int i,errcode;
    
  /*1. left shift*/
  shifted->from_res_no--;
  shifted->to_res_no--;

  /*2. Initialize the arrays a and c before filling them in the next routine*/
  for(i=0; i<NUMPOINTS; i++)   
    atmp[i][0] = atmp[i][1] = atmp[i][2] = BOUND; 
  for(i=0; i<NUMPOINTS/2; i++)   
    ctmp[i][0] = ctmp[i][1] = ctmp[i][2] = BOUND;

  /*3. Read file */
  if( (errcode = process_pdbfile(shifted,atmp,ctmp))==SUCCESS) {
    /*errcode = match_with_template(atmp,ctmp,new_type);*/
    errcode = (*func)(atmp,ctmp,new_type);
    if(errcode == SUCCESS) {         
      if(!strcmp(new_type,"HELI"))
	lefthelix=1;
    }else {
      printf("Could not match left-shifted sequence\n");
      strcpy(type,"ERROR");
      return -1;
    }
  } else {
    printf("Could not read shifted left-shiftedsequence from pdb file \n");
    strcpy(type,"ERROR");
    return -1;
  }
  
  /*4. restore seq_data structure */
  shifted->from_res_no++;
  shifted->to_res_no++;
  
  /*Right Shift*/
  if(!lefthelix) {
    /* rightshift */
    shifted->from_res_no++;
    shifted->to_res_no++;
    
    /*Initialize the arrays a and c before filling them in the next routine*/
    for(i=0; i<NUMPOINTS; i++)   
      atmp[i][0] = atmp[i][1] = atmp[i][2] = BOUND; 
    for(i=0; i<NUMPOINTS/2; i++)   
      ctmp[i][0] = ctmp[i][1] = ctmp[i][2] = BOUND;
    
    if( (errcode = process_pdbfile(shifted,atmp,ctmp))==SUCCESS) {
      /*errcode = match_with_template(atmp,ctmp,new_type);*/
      errcode = (*func)(atmp,ctmp,new_type);
      if(errcode == SUCCESS) {       
	if(!strcmp(new_type,"HELI"))
	  righthelix=1;
      }else {
	printf("Could not match right-shifted sequence\n");
	strcpy(type,"ERROR");
	return -1;
      }
    } else {
      printf("Could not read shifted right-shifted sequence from pdb file \n");
      strcpy(type,"ERROR");
      return -1;
    }
  }
  if(lefthelix || righthelix) {
    strcpy(type,"HIRR");
    return 1;
  }
  else
    strcpy(type,"TIRR");
  return 0;
}
  
