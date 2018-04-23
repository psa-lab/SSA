# SSA

## Introduction

SSA is a tool to assign the secondary structure of a peptide from its atomic coordinates in PDB format based on their superposition with sequences of ideal secondary structure. SSA was developed by members of the Protein Structural Analysis and Design Laboratory at Michigan State University, and was supported by funding from the MSU Vice Provost for Libraries, Computing, & Technology and from the [Department of Biochemistry](https://bmb.natsci.msu.edu/).

SSA is designed to work on 4-residue peptide fragments (approximately reflecting the periodicity of helices and reverse turns and two repeats within beta-strand) and is most often coupled for use with [Sequery](https://github.com/psa-lab/Sequery), which performs sequence pattern searches of the [Protein Data Bank](http://rcsb.org) (PDB). For more information and software availability on Sequery, see https://github.com/psa-lab/Sequery.

## References

References related to the applications of SSA include the following:

- L. Craig, P. C. Sanschagrin, A. Rozek, S. Lackie, L. A. Kuhn, and J. K. Scott (1998) The Role of Structure in Antibody Cross-Reactivity Between Peptides and Folded Proteins (pdf), J. Mol. Biol., 281, 183-201. [Link to PDF](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/CraigJMB1998.pdf)

- M. A. Kron, S. Kutlesa, A. Hendrick, A. Liu, J. Leykam, S. Cichenowicz, and L. A. Kuhn (2008) Using Structural Analysis to Generate Parasite-Selective Monoclonal Antibodies (pdf), Protein Science 17, 983-989. [Link to PDF](http://www.kuhnlab.bmb.msu.edu/publication_papers/pdf/Kron_et_al.ProtSci08.pdf)

For specific information related to the SSA algorithm, please see the section on *Algorithmic Details* below.

Usage information for SSA can be found in the section on *Running SSA* below.

## Installation

SSA was designed and written for use on a UNIX workstation with a C-compiler and tested on Sun and Silicon Graphics workstations, though it should work on any system with a C-compiler installed and access to the [Protein Data Bank](http://rcsb.org). Implicit declaration errors may occur during compilation with the GNU C-compiler (gcc), but these errors do not affect program compilation or execution. Additional errors may be reported during compilation with the SGI C-compiler, but these will also have no effect on compilation or execution.

To install SSA 3.0, perform the following steps:

- Download this repository to your local machine (Note: Not tested on 64-bit machines)
- Place this file in the directory into which you wish to install the software and unzip it
- To compile the software, type the following command from the `SourceFiles` directory: `make install`
- To finish installation, edit the `SourceFiles/runssa` script in the following manner:
    + Set `PDBHOME` to the directory containing the PDB database on your system
    + Set `ssahome` to the directory where you installed SSA using the full and explicit path.
- To run SSA, use the `runssa` script in the `installation/bin` directory. For example, if SSA is installed in `/usr/local/biochem/SSA`, then you would use the command: `/usr/local/biochem/SSA/bin/runssa` to run the software. Further information on the use of SSA can be found in the section on *Running SSA* below.

**Troubleshooting Shell Scripts**

`ssa/SourceFiles/runssa` and `ssa/share/genbcl` are C-shell scripts. They assume the C-shell is installed in the standard location as `/bin/csh`. If you have difficulty running this script, try editing the first line of each script to point to the correct shell location for your system. Also verify these files are set with the proper executable permissions.

## Running SSA

SSA can be run in two different modes. The standard mode takes an output file produced by Sequery, compares the structure of each sequence match to a set of standard templates to assign a secondary structure type, then outputs the resulting assignments and summary information for the complete Sequery file.

**Syntax**: 

    runssa SequeryOutputFile SSAOutputFile


| **Argument**        | **Description**                                                                                         |
|---------------------|---------------------------------------------------------------------------------------------------------|
| `SequeryOutputFile` | file resulting from running Sequery; each line corresponds to a sequence of residues within a PDB file. |
| `SSAOutputFile`     | file where you wish to place the SSA results                                                            |


**Example of SSA Input (Sequery Output) File:**

    1inp _  166 to  169 -> ikgsADITpnqg matching .DIT
    2bpa 1  230 to  233 -> tsydADNRpllv matching A.NR
    1cpc A   12 to   15 -> avaaADSQgrfl matching .D5Q

**Example of SSA Standard Mode Output File** (described in detail below in Output Files):

    1inp _  166 to  169 -> ikgsADITpnqg matching .DIT  IRRE
    2bpa 1  230 to  233 -> tsydADNRpllv matching A.NR  IRRE
    1cpc A   12 to   15 -> avaaADSQgrfl matching .D5Q  HELI


SSA can also be run with an additional, user-supplied template file in the secondary structure template set, for instance, when you want to compare Sequery matches to a nonstandard structure. Running in user-template mode will result in the RMSD between each of the Sequery matches and the user template being reported, in addition to the standard template matching assignments and summary.

**Syntax**:

    runssa SequeryOutputFile UserTemplateFile SSAOutputFile

| **Argument**        | **Description**                                                                                                |
|---------------------|----------------------------------------------------------------------------------------------------------------|
| `SequeryOutputFile` | file resulting from running Sequery; each line corresponds to a sequence of residues within a PDB file.        |
| `UserTemplateFile`  | file which contains a user-supplied template in PDB format (see ssa/Examples/usertemplate.pdb for an example.) |
| `SSAOutputFile`     | file where you wish to place the SSA results                                                                   |

**Example of SSA User-template Mode Output File** (Note the addition of and RMSD in the second to last column, which is the main-chain (N, C-alpha, C, O) RMSD between the match and the user-template tetrapeptide) and the * in final column, indicating it was a structural match to the supplied template.)

    1inp _  166 to  169 -> ikgsADITpnqg matching .DIT  IRRE  2.13
    2bpa 1  230 to  233 -> tsydADNRpllv matching A.NR  IRRE  1.96
    1cpc A   12 to   15 -> avaaADSQgrfl matching .D5Q  HELI  0.64 *

    SSA can also be run in **Debug-mode** by giving the argument `DEBUG` before other arguments. This mode will output all distance and matching matrices and is probably of limited use for the user.

## Example Files

Example input, output, and log files have been included in `ssa/Examples`directory.

| **File**                     | **Description**                                                                                                                                          |
|------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| Example.sequery              | A Sequery output file                                                                                                                                    |
| Example.standardmode.out     | SSA output file generated by running in standard mode (runssa Example.sequery Example.standardmode.out)                                                  |
| ssa.standardmode.log         | Logfile associated with the above run, showing error messages that occur (explained further in Diagnostics)                                              |
| Example.usertemplatemode.out | SSA output file generated by running in user-template mode with the usertemplate.pdb. (runssa Example.sequery usertemplate.pdb Example.standardmode.out) |
| ssa.usertemplatemode.log     | Logfile for the above run                                                                                                                                |
| usertemplate.pdb             | An example user-template file                  

## Output Files


**Standard Mode** -- SSA will generate an output file containing a line for each Sequery match, listing the PDB code, chain identifier, numerical residue range, 4-residue peptide matching sequence (in upper case) along with 4 residues of sequence on either side (in lower case), the sequence pattern matched, and the assigned secondary structure. The secondary structure codes are as follows:


| **Secondary Structure Code** | **Description**                                                                                                                                              |
|------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| HELI                         | Alpha-helix                                                                                                                                                  |
| HIRR                         | Irregular Helix -- the first or last turn of a helix, or a location of a helical kink                                                                        |
| EXTE                         | Extended -- matching a beta-strand template (though not necessarily associated with a beta-sheet tertiary structure, since hydrogen bonding is not checked). |
| IRRE                         | Irregular -- tetrapeptide that is not a good structural match with any of the template structures, including turns                                           |
| **Turns**                    |                                                                                                                                                              |
| TUR1                         | Type 1 Turn                                                                                                                                                  |
| TT1P                         | Type 1' Turn                                                                                                                                                 |
| TUR2                         | Type 2 Turn                                                                                                                                                  |
| TT2P                         | Type 2' Turn                                                                                                                                                 |
| TUR8                         | Type 8 Turn                                                                                                                                                  |
| TIRR                         | Irregular Turn -- Does not match very closely with a turn template, but i and i+3 C-alpha's are within 7.0 Å                                                 |
| **Misc**                     |                                                                                                                                                              |
| ERROR                        | Error -- an error occurred during processing of this sequence match. Check ssa.log for details.                                                              |


In addition to individual lines for each sequence match, a summary of results is posted at the bottom. This summary includes basic information about the processing and statistics for the different types of secondary structures assigned, including sums and percentages for each general class, excluding processing errors. If errors occured, check the ssa.log file for more details.


**Example Full Standard Mode Output File:** (taken from `ssa/Examples/Example.standardmode.out`, with most of the output lines removed for brevity)

    # Sequery Filename : Example.sequery
     
    #    C
    #    h
    #    a
    #    i   Residue                                         User User
    #PDB n      Range      Sequence                    Type  RMSD Match?
    #------------------------------------------------------------------
    1pxt B   40 to   43 -> vvivAANRsaig matching A.NR  EXTE
    1tiv _   21 to   24 -> qpktACNRchck matching A6NR  EXTE
    1pox A  324 to  327 -> iavlADAQktla matching .D5Q  IRRE
    2bpa 1  230 to  233 -> tsydADNRpllv matching A.NR  IRRE
    1cpc A   12 to   15 -> avaaADSQgrfl matching .D5Q  HELI
    1hpm _  131 to  134 -> mkeiAEAYlgkt matching .E.Y  HELI
    1aor A  602 to  605 -> elgiAEFY     matching .E.Y  ERROR
    1aor A  407 to  410 -> syrlAESYghpe matching .E.Y  HELI
    1mxa _  213 to  216 -> pilpAEWLtsat matching .EWL  TUR1
    1fnc _  223 to  226 -> mkekAPDNfrld matching AP..  TUR1
     
    ####################################STATISTICS##############################
    # Total number of sequery lines processed is:           10
    # Total number of Errors (excluded from statistics):    1  (10.00 %)
    # # of Regular Helix Types: 3    (33.33 %)
    # # of Irreg Helix Types:   0    (0.00 %)
    # Total # of Helix Types:           3    (33.33 %)
    #
    # # of Extended Types:              2    (22.22 %)
    #
    # # of Type1 Turns:         2    (22.22 %)
    # # of Type2 Turns:         0    (0.00 %)
    # # of Type1' Turns:        0    (0.00 %)
    # # of Type2' Turns:        0    (0.00 %)
    # # of Type8 Turns:         0    (0.00 %)
    # # of Irreg Turns:         0    (0.00 %)
    # Total # Of Turns:                 2    (22.22 %)
    #
    # # of Irreg Types *:               2    (22.22 %)
    #
    #* # of Irreg Types excludes Irreg Helices and Irreg Turns


**User-Template Mode** -- In addition to the above results, SSA will also output the main-chain RMSD between each of the sequence matches and the user-template and whether it matches within the threshold RMSD of 0.75 Å, which is denoted by an asterisk (*) under the ``User Match?'' column heading. A summary of the number and percentage of matching sequences is also included.

**Example Full User-Template Mode Output File:** (taken from ssa/Examples/Example.usertemplate.mode,out, with most of the output lines removed for brevity)

    # Sequery Filename : Example.sequery
    # User Template : usertemplate.pdb
     
    #    C
    #    h
    #    a
    #    i   Residue                                         User User
    #PDB n      Range      Sequence                    Type  RMSD Match?
    #------------------------------------------------------------------
    1pxt B   40 to   43 -> vvivAANRsaig matching A.NR  EXTE  2.71
    1tiv _   21 to   24 -> qpktACNRchck matching A6NR  EXTE  2.55
    1pox A  324 to  327 -> iavlADAQktla matching .D5Q  IRRE  2.25
    2bpa 1  230 to  233 -> tsydADNRpllv matching A.NR  IRRE  1.96
    1cpc A   12 to   15 -> avaaADSQgrfl matching .D5Q  HELI  0.64 *
    1hpm _  131 to  134 -> mkeiAEAYlgkt matching .E.Y  HELI  0.57 *
    1aor A  602 to  605 -> elgiAEFY     matching .E.Y  ERROR
    1aor A  407 to  410 -> syrlAESYghpe matching .E.Y  HELI  0.66 *
    1mxa _  213 to  216 -> pilpAEWLtsat matching .EWL  TUR1  0.75
    1fnc _  223 to  226 -> mkekAPDNfrld matching AP..  TUR1  0.57 *
     
    ####################################STATISTICS##############################
    # Total number of sequery lines processed is:           10
    # Total number of Errors (excluded from statistics):    1  (10.00 %)
    # # of Regular Helix Types: 3    (33.33 %)
    # # of Irreg Helix Types:   0    (0.00 %)
    # Total # of Helix Types:               3        (33.33 %)
    #
    # # of Extended Types:                  2        (22.22 %)
    #
    # # of Type1 Turns:         2    (22.22 %)
    # # of Type2 Turns:         0    (0.00 %)
    # # of Type1' Turns:        0    (0.00 %)
    # # of Type2' Turns:        0    (0.00 %)
    # # of Type8 Turns:         0    (0.00 %)
    # # of Irreg Turns:         0    (0.00 %)
    # Total # Of Turns:                     2        (22.22 %)
    #
    # # of Irreg Types *:                   2        (22.22 %)
    #
    #* # of Irreg Types excludes Irreg Helices and Irreg Turns
    # Number of User Matches:  4     and % that matched is:  44.44


## Diagnostics


**Compilation Errors**

As stated above in the Installation instructions, SSA may cause compilation warnings when compiled with GNU's gcc compiler. These errors have no effect on execution.

**Run-Time Errors**

- `ERROR: Insufficient Input Data`. Couldnt grab residues on the lower and upper boundary of the sequence. - PDB code is "1aor" and sequence is "elgiAEFY"
Explanation -- This error was generated with a sequence match to the sequence AEFY at the C-terminus of the structure with PDB code 1aor. The error given indicates that there are not sufficient residues to test the structure of flanking regions for helicity.

- `ERROR: Backbone atoms from PDB file could not be matched with the template for superposition, probably due to missing or extra atoms. - PDB code is "1efg" and sequence is "vwrqAEKYkvpr"`
Explanation -- This error was generated with a sequence match to the sequence AEKY in the structure with PDB code 1efg. This error indicates that there are insufficient atoms to perform the required superposition to calculate the RMSD. SSA uses a superposition of all 4 backbone atoms (nitrogen, alpha-carbon, carboyl-carbon & carbonyl-oxygen). Since superposition requires a one-to-one atom correspondance, it will not work correctly if there are fewer atoms or additional atoms in one of the sets. Most cases of this error occur for alpha-carbon-only PDB structures.

## Algorithmic Details

SSA uses superpositional analysis to assign secondary structure to short peptidyl segments (4 residues in its current implementation). The basic algorithmic procedure is as follows:

#### 1)

Superimpose the tetrapeptide in question onto each of the 7 template teptrapeptides (alpha-helix, beta-strand, type I turn, type II turn, type I' turn, type II' turn, type VIII turn). The alpha-helix and beta-strand templates were generated using the phi-psi values in InsightII (MSI, San Diego, CA) and the reverse turn templates were generated using the phi-psi angles for average turns as reported in: 

Hutchinson EG, Thornton JM
"A revised set of potentials for beta-turn formation in proteins"
Prot. Sci., 3(12) (Dec): 2207-2216 (1994)

#### 2)

Assign the structure type based on a threshold of < 0.75 Å main-chain (N,C-alpha, C, O) RMSD using the following algorithm: 

- If the tetrapeptide matches the beta-strand template, assign as beta-strand (EXTE).
- If the tetrapeptide matches the helix template, check the structure of the tetrapeptide found by shifting one residue upstream and the structure of the tetrapeptide found by shifting one residue downstream against the helix template. 
    - If one OR the other of the shifted peptides also match the helix template, assign as helix (HELI).
    - If NEITHER of the shifted peptides match the helix template, check futher

- If the tetrapeptide matches one of the turn templates, assign as turn type with the lowest RMSD (type 1 = TUR1; type 1' = TT1P; type 2 = TUR2; type 2' = TT2P; or type 8 = TUR8)
- If the tetrapeptide matches none of the templates, calculate the C-alpha(i) to C-alpha(i+3) distance. 
    + If this distance is LESS than 7.0 Å, check the structure of the tetrapeptide found by shifting one residue upstream and the structure of the tetrapeptide found by shifting one residue downstream against the helix template. 
        - If one OR the other of the shifted peptides also match the helix template, assign as irregular helix (HIRR).
        - If NEITHER match the helix template, assign as irregular turn (TIRR).
    - If this distance is GREATER than 7.0 Å, assign as irregular (IRRE).


## More Information

Scientific inquries concerning SSA should be directed to Leslie Kuhn at: kuhnlab@msu.edu