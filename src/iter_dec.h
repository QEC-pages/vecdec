#ifndef ITER_DEC_H
#define ITER_DEC_H
/**
 * @file vecdec.h
 *
 * @brief vecdec - a simple vectorized decoder
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2022 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 */
#ifdef __cplusplus
extern "C"{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

  /** structure to hold global variables */

  typedef struct PARAMS_T {
    int nrows; /** rows in `H` or `r` (set in the `input file`) */
    int ncws;  /** how many codewords `k` (set in the `input file`) */
    int n;     /** columns in H `n` (set in the `input file`) */
    int steps; /** number of BP steps, default: `50` */
    int ntot;  /** total number of syndromes to generate (default: `1`) */
    int nfail; /** when non-zero, num of fails to terminate the run (default: `0`, do not terminate) */
    int lerr;  /** OSD level (default: `0`) */
    int mode;  /** operation mode bitmap, see help */
    int debug; /** `debug` information */ 
    char *finH; /** `input file` name for H (if input separately or a classical code) */
    char *finP; /** `input file` name for P (if input separately or a classical code) */
    char *fdem; /** `input file` name for detector error model (`DEM`) */
    char *fdet; /** `input file` name for detector events */
    char *fobs; /** `input file` name for observables */
    int seed;  /** rng `seed`, set=0 for automatic */
    double *vP; /** probability vector (total of `n`) */
    double *vLLR; /** vector of LLRs (total of `n`) */
    int nzH, nzL; /** count of non-zero entries in `H` and `L` */
    csr_t *mH; /** sparse version of H (by rows) */
    csr_t *mHt; /** sparse version of H (by columns) */
    csr_t *mL; /** sparse version of L (by rows) */
    csr_t *mLt; /** sparse version of L (by columns) */
    csr_t *mG; /** sparse version of generator matrix `G` (by rows) */
    /** rows of `G` orthogonal to rows of both `H` and `L` */
    int maxJ;  /** memory to initially allocate for local storage */
    double LLRmin;
    double LLRmax;
  } params_t;

  extern params_t prm;
  
  /** 
   * @brief The help message. 
   * 
   * @todo: This has to be checked and updated, especially the `debug` options.
   */
#define USAGE                                                           \
  "%s:  iter_dec - BP decoder for quantum and classical codes\n"	\
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Command line arguments are processed in the order given.\n"	\
  "\t Supported parameters:\n"						\
  "\t --help\t: give this help (also '-h' or just 'help')\n"            \
  "\t fdem=[string]\t: name of the input file with detector error model\n" \
  "\t finP=[string]\t: input file for probabilities (mtx or a column of doubles)\n" \
  "\t finH=[string]\t: input file with parity check matrix (mtx format)\n" \
  "\t fdet=[string]\t: input file with detector events (01 format)\n"   \
  "\t fobs=[string]\t: file with observables (01 matching lines in fdet)\n" \
  "\t\t\t (space is OK in front of file name to enable shell completion)\n" \
  "\t steps=[integer]\t: num of BP decoding steps (default: 50)\n"	\
  "\t lerr =[integer]\t: OSD level (0, no OSD)\n"			\
  "\t ntot =[integer]\t: total syndromes to generate (default: 1)\n"	\
  "\t nfail=[integer]\t: total fails to terminate (0, do not terminate)\n" \
  "\t seed= [integer]\t: RNG seed or use time(NULL) if 0 (default)\n"	\
  "\t mode= [integer]\t: BP operation mode bitmap (default: 0)\n"	\
  "\t\t* 0: use basic parallel BP\n"					\
  "\t\t\t read detector events from file 'fdet' if given, otherwise\n"  \
  "\t\t\t generate 'ntot' detector events and matching observable flips;\n" \
  "\t\t\t read observable flips from file 'fobs' if given\n"            \
  "\t debug=[integer]\t: bitmap for aux information to output (default: 1)\n" \
  "\t\t*   0: clear the entire debug bitmap to 0.\n"                    \
  "\t\t*   1: output misc general info (on by default)\n"		\
  "\t\t*   2: output matrices for verification\n"                       \
  "\t See program documentation for input file syntax.\n"               \
  "\t Multiple `debug` parameters are XOR combined except for 0.\n"	\
  "\t Use debug=0 as the 1st argument to suppress all debug messages.\n"

#ifdef __cplusplus
}
#endif


#endif /* ITER_DEC_H */
