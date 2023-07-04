#ifndef VECDEC_H
#define VECDEC_H
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

/**< minimum probability for LLR calculations */
#define MINPROB (1e-7)  

/** @brief structure to read in one line of data */
typedef struct ONE_PROB_T {
  double p; /**< probability */
  int n1;   /**< number of entries in H column */
  int n2;   /**< number of entries in L column */
  int idx[]; /**< flexible array to store `n1+n2` entries */
} one_prob_t;

/** @bried helper structure to sort by probabilities */
typedef struct IPPAIR_T {int index; double prob; } ippair_t; 

/** structure to hold global variables */

typedef struct PARAMS_T {
  int nrows; /** rows in `H` or `r` (set in the input file) */
  int ncws;  /** how many codewords `k` (set in the input file) */
  int n;     /** columns in H `n` (set in the input file) */
  int colw;  /** max column weight (default: `10`, increase at the command line if needed) */ 
  int steps; /** how many random decoding steps, default: `1` */
  int nvec;  /** max number of syndromes to process in one bunch */
  int mode;  /** mode information */
  int debug; /** `debug` information */ 
  char *fin; /**< `input file` name for error model */
  int seed;  /**< rng `seed`, set=0 for automatic */
  double *vP; /**< probability vector (total of `n`) */
  double *vLLR; /**< vector of LLRs */
  int nzH, nzL; /**< count of non-zero entries in `H` and `L` */
  csr_t *mH; /**< sparse version of H (by rows) */
  csr_t *mHt; /**< sparse version of H (by rows) */
  csr_t *mLt; /**< sparse version of L (by columns) */
  int maxJ;  /** memory to initially allocate for local storage */
} params_t;

extern params_t prm;

/** 
 * @brief The help message. 
 * 
 * @todo: This has to be checked and updated, especially the `debug` options.
 */
#define USAGE                                                           \
  "%s: a simple vectorized random information set decoder.\n"           \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Command line arguments are processed in the order given.\n"	\
  "\t Supported parameters:\n"						\
  "\t --help\t: give this help (also '-h' or just 'help')\n"            \
  "\t f=[string]\t: name of the input file with the error model\n"      \
  "\t steps=[integer]\t: how many random window decoding steps (default: 1)\n" \
  "\t nvec =[integer]\t: max vector size for decoding (default: 16)\n"  \
  "\t colw =[integer]\t: max column weight in stacked H+L (default: 10)\n" \
  "\t seed= [integer]\t: RNG seed or use time(NULL) if 0 (default)\n"	\
  "\t mode= [integer]\t: bitmap for operation mode (default: 0)\n"      \
  "\t\t*   0: clear the entire mode bitmap to 0.\n"                     \
  "\t\t*   1: run as standalone program (no stdin)\n"                   \
  "\t\t*   2: generate random errors and syndrome vectors to stdout\n"  \
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

#endif /* VECDEC_H */
