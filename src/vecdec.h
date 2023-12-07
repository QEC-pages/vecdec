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
#include "uthash.h" /** hashing storage macros */
#include  "qllr.h" 

  /**< @brief structure to hold sparse vectors in a hash */
  typedef struct ONE_VEC_T {
    UT_hash_handle hh;
    double energ; /**< sum of LLRs */
    int weight; /**< number of integers in the list */
    int cnt; /** how many times this vector was encountered */
    //  size_t len; /** `weight*sizeof(int)` (is this really needed?) */
    int arr[0]; /** array of `weight` integers, the actual key  */
  } one_vec_t;

  /** @brief structure to read in one line of data */
  typedef struct ONE_PROB_T {
    double p; /**< probability */
    int n1;   /**< number of entries in H column */
    int n2;   /**< number of entries in L column */
    int idx[]; /**< flexible array to store `n1+n2` entries */
  } one_prob_t;

  /** @bried helper structure to sort by probabilities */
  typedef struct IPPAIR_T {int index; double prob; } ippair_t; 


typedef enum EXTR_T { TOTAL, CONV_TRIVIAL, CONV_BP, CONV_BP_AVG,
  SUCC_TRIVIAL, SUCC_BP, SUCC_OSD0, SUCC_OSD1, SUCC_OSD2,
  SUCC_TOT, EXTR_MAX } extr_t;

/** various success counters */
extern long long int cnt[EXTR_MAX];
extern long long int iter1[EXTR_MAX]; /** sums of BP iteration numbers */
extern long long int iter2[EXTR_MAX]; /** sums of BP iteration numbers squared */

  
  /** structure to hold global variables */

  typedef struct PARAMS_T {
    int nchk; /** rows in `H` or `r` (set in the `input file`) */
    int ncws;  /** how many codewords `k` (set in the `input file`) */
    int nvar;     /** columns in H `n` (set in the `input file`) */
    int steps; /** number of random window or BP decoding steps, default: `50` */
    int nvec;  /** max number of syndromes to process in one bunch (default: `16`) */
    int ntot;  /** total number of syndromes to generate (default: `1`) */
    int nfail; /** when non-zero, num of fails to terminate the run (default: `0`, do not terminate) */
    int swait; /** gauss decoding steps with no vectors changed to stop (default: `0`, do not stop) */
    int lerr;  /** local search after gauss up to this weight (default: `0`) */
    int mode;  /** operation mode, see help */
    int submode; /** additional options, see help */
    int d1, d2, d3; /** QLLR parameters for BP */
    int use_stdout; /** with mode=3 */
    int debug; /** `debug` information */ 
    char *finH; /** `input file` name for Hx=H (if input separately or a classical code) */
    char *finL; /** `input file` name for Lx=L (if input separately or a classical code) */
    char *finG; /** `input file` name for Hz=G (must use separate input) */
    char *finP; /** `input file` name for P (if input separately or a classical code) */
    char *fout; /** `output file name`  header for files creaded with `mode=3` */
    char *fdem; /** `input file` name for detector error model (`DEM`) */
    char *fdet; /** `input file` name for detector events */
    char *fobs; /** `input file` name for observables */
    char *ferr; /** `input file` name for error vectors */
    int classical; /** `1` if this is a classical code? */
    int internal; /** `1` to generate obs/det internally, `2` to generate from `err` file */
    int seed;  /** rng `seed`, set=0 for automatic */
    double useP; /** global error probability `overriding` values in the `DEM` file (default: 0, no override) */
    double *vP; /** probability vector (total of `n`) */
    double *vLLR; /** vector of LLRs (total of `n`) */
    int nzH, nzL; /** count of non-zero entries in `H` and `L` */
    csr_t *mH; /** sparse version of H (by rows) */
    csr_t *mHt; /** sparse version of H (by columns) */
    csr_t *mL; /** sparse version of L (by rows) */
    csr_t *mLt; /** sparse version of L (by columns) */
    csr_t *mG; /** sparse version of generator matrix `G` (by rows) */
    /** rows of `G` orthogonal to rows of both `H` and `L` */
    //    int maxJ;  /** memory to initially allocate for local storage */
    qllr_t LLRmin;
    qllr_t LLRmax;
    one_vec_t *codewords; /** `hash table` with found codewords */
    long int num_cws; /** `number` of codewords in the `hash` */
    FILE *file_err;
    FILE *file_det;
    FILE *file_obs;
    int line_err; /** current line of the err file */
    int line_det; /** current line of the det file */
    int line_obs; /** current line of the obs file */
    mzd_t *mE;
    mzd_t *mHe;
    mzd_t *mLe;
    mzd_t *mHeT;
    mzd_t *mLeT;
  } params_t;

  extern params_t prm;

  /** @brief helper function to sort `ippair_t`
   *  use `qsort(array, len, sizeof(ippair_t), cmp_ippairs);`
   */
  static inline int cmp_ippairs(const void *a, const void *b){
    const double pa=((ippair_t *) a) -> prob;
    const double pb=((ippair_t *) b) -> prob;
    if (pa<pb)
      return +1;
    else if (pa>pb)
      return -1;
    return 0;
  }

  /** @brief helper function to sort `int`
   *  use `qsort(array, len, sizeof(rci_t), cmp_rci_t);`
   */
  static inline int cmp_rci_t(const void *a, const void *b){
    const rci_t va= *((rci_t *) a);
    const rci_t vb= *((rci_t *) b);
    return va-vb;
#if 0
    if (va<vb)
      return +1;
    else if (va>vb)
      return -1;
    return 0;
#endif
  }
  
  void cnt_out(int print_banner);
  void cnt_update(extr_t which, int iteration);
  void out_llr(const char str[], const int num, const qllr_t llr[]);
  
  int syndrome_check(const qllr_t LLR[], const mzd_t * const syndrome,
		     const csr_t * const H,
		     [[maybe_unused]] const params_t * const p);
  void cnt_update(extr_t which, int iteration);
  
  int do_parallel_BP(qllr_t * outLLR, const mzd_t * const srow,
		   const csr_t * const H, const csr_t * const Ht,
		     const qllr_t LLR[], const params_t * const p);

  /** 
   * @brief The help message. 
   * 
   * @todo: This has to be checked and updated, especially the `debug` options.
   */
#define USAGE                                                           \
  "%s:  vecdec - vectorized decoder and LER estimator\n"                \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Command line arguments are processed in the order given.\n"	\
  "\t Supported parameters:\n"						\
  "\t --help\t\t: give this help (also '-h' or just 'help')\n"		\
  "\t --morehelp\t: give more help\n"					\
  "\t fout=[string]\t: header for output file names ('tmp', see 'mode=3')\n" \
  "\t finP=[string]\t: input file for probabilities (mm or a column of doubles)\n" \
  "\t finH=[string]\t: file with parity check matrix (mm or alist)\n"	\
  "\t finG=[string]\t: file with dual check matrix (mm or alist)\n"	\
  "\t finL=[string]\t: file with logical dual check matrix (mm or alist)\n" \
  "\t fdem=[string]\t: name of the input file with detector error model\n" \
  "\t ferr=[string]\t: file with error vectors (01 format), or\n"	\
  "\t fdet=[string]\t: input file with detector events (01 format)\n"   \
  "\t fobs=[string]\t: file with observables (01 matching lines in fdet)\n" \
  "\t\t\t (space is OK in front of file name to enable shell completion)\n" \
  "\t steps=[integer]\t: num of random window decoding steps (default: 1)\n" \
  "\t lerr =[integer]\t: local search level after gauss (0, no search)\n" \
  "\t swait=[integer]\t: steps w/o new errors to stop (0, do not stop)\n" \
  "\t nvec =[integer]\t: max vector size for decoding (default: 16)\n"  \
  "\t\t\t (list size in distance or energy calculations)\n"             \
  "\t ntot =[integer]\t: total syndromes to generate (default: 1)\n"	\
  "\t nfail=[integer]\t: total fails to terminate (0, do not terminate)\n" \
  "\t seed= [integer]\t: RNG seed or use time(NULL) if 0 (default)\n"	\
  "\t mode= [int.int]\t: operation mode.submode (default: 0.0)\n"	\
  "\t\t* 0: use basic vectorized (random information set) decoder\n"	\
  "\t\t\t read detector events from file 'fdet' if given, otherwise\n"  \
  "\t\t\t generate 'ntot' detector events and matching observable flips;\n" \
  "\t\t\t read observable flips from file 'fobs' if given\n"            \
  "\t\t* 1: differrent BP versions\n"					\
  "\t\t\t .0 basic parallel BP\n"					\
  "\t\t* 2: generate most likely fault vectors, estimate Prob(Fail)\n"  \
  "\t\t\t generate up to 'ntot' unique min-energy fault vectors\n"	\
  "\t\t\t use up to 'steps' random window decoding steps unless no new\n" \
  "\t\t\t fault vectors have been found for 'swait' steps.\n"           \
  "\t\t\t Keep vectors of weight up to 'nfail' above min weight found\n" \
  "\t\t* 3: Read in the DEM file and output the corresponding \n"	\
  "\t\t\t H, G, and L matrices and the probability vector P.\n"		\
  "\t\t\t Use 'fout=' command line argument to generate file names\n"	\
  "\t\t\t ${fout}H.mmx, ${fout}G.mmx, ${fout}L.mmx, and ${fout}P.mmx\n"	\
  "\t\t\t with 'fout=stdout' all output is sent to 'stdout'\n"		\
  "\t debug=[integer]\t: bitmap for aux information to output (default: 1)\n" \
  "\t\t*   0: clear the entire debug bitmap to 0.\n"                    \
  "\t\t*   1: output misc general info (on by default)\n"		\
  "\t\t*   2: output matrices for verification\n"                       \
  "\t\t\t see the source code for more options\n"			\
  "\t See program documentation for input file syntax.\n"               \
  "\t Multiple `debug` parameters are XOR combined except for 0.\n"	\
  "\t Use debug=0 as the 1st argument to suppress all debug messages.\n"

#define MORE_HELP							\
  "   Matrices used by %s:\n"						\
  "\t Think about a CSS code with binary generator matrices Hx=H, Hz=G,\n" \
  "\t and logical-operator generating matrices Lx=L and Lz=K.  These\n"	\
  "\t matrices have n columns each and satisfy orthogonality properties\n" \
  "\t\t Hx*Hz^T=0, Lx*Hz^T=0, Lx*Lz^T=0, and Lx*Lz^T=I (identity matrix).\n" \
  "\t We are trying to correct binary Z errors, whose probabilities\n"	\
  "\t are given by the error vector P (or can be overridden with\n"	\
  "\t the command-line parameter 'Perr').\n"				\
  "   A detector error model (DEM) file, in the format produced by Stim,\n" \
  "\t contains matrices H and Lx, and the error probability vector P.\n" \
  "\t The same code can be obtained by specifying the matrices \n"	\
  "\t independently, via separate files.\n"				\
  "   The dual CSS matrix Hz can be specified instead of Lx.\n"	\
  "\t In such a case, the internal error generator must be used\n"	\
  "\t (an attempt to specify 'fdet' and 'fobs' files will result in an error).\n" \
  "   For a classical code, specify only the parity check matrix Hx=H.\n" \
  "\t In this case G matrix is trivial (has zero rank), and\n"		\
  "\t Lx is set to identity matrix.  \n"				\
  "\t Only the internal error generator can be used for the classical code\n"
  
  
#ifdef __cplusplus
}
#endif

#endif /* VECDEC_H */
