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
#include  "qllr.h" 


  typedef enum EXTR_T { TOTAL, CONV_TRIVIAL, CONV_BP, CONV_BP_AVG, CONV_BP_TOT,
    SUCC_TRIVIAL, SUCC_BP, SUCC_OSD, SUCC_TOT, EXTR_MAX } extr_t;
  
  /** various success counters */
  extern long long int cnt[EXTR_MAX];
  extern long long int iter1[EXTR_MAX]; /** sums of BP iteration numbers */
  extern long long int iter2[EXTR_MAX]; /** sums of BP iteration numbers squared */
  
  
  /** structure to hold global variables */
  typedef struct PARAMS_T {
    int nchk; /** rows in `H` or `r` (set in the `input file`) */
    int ncws;  /** how many codewords `k` (set in the `input file`) */
    int nvar;     /** columns in H `n` (set in the `input file`) */
    int pads; /** if `1`, pad detector events (syndrome vectors) read from
		     `fdet` file to required number of columns; if `0`, do
		     no padding */
    int rankH;
    int rankG;
    int rankL;
    int steps; /** number of random window or BP decoding steps, default: `50` */
    int nvec;  /** max number of syndromes to process in one bunch (default: `16`) */
    long long int ntot;  /** total number of syndromes to generate (default: `1`) */
    long long int nfail; /** when non-zero, num of fails to terminate the run (default: `0`, do not terminate) */
    long long int maxC; /** when non-zero, max number of codewords (default: `0`, no maximum) */
    int dmin; /** if non-zero, terminate distance calculation immediately when a vector 
		  of weight `d`<=`dmin` is found, return `-d` (default: 0) */
    int swait; /** gauss decoding steps with no vectors changed to stop (default: `0`, do not stop) */
    int lerr;  /** local search after gauss up to this weight (default: `-1`, no OSD) */
    int maxosd;  /** max column for OSD2 and above (default: `100`) */
    double bpalpha; /** average LLR parameter; multiply old LLR by `bpalpha` new by `1-bpalpha` (default: `0.5`) */
    int bpretry; /** for each syndrome try BP up to this many times  (default: `1`) */
    int mode;  /** operation mode, see help */
    int submode; /** additional options, see help */
    int d1, d2, d3; /** QLLR parameters for BP */
    int use_stdout; /** with mode=3 */
    int debug; /** `debug` information */ 
    char *finH; /** `input file` name for Hx=H (if input separately or a classical code) */
    char *finL; /** `input file` name for Lx=L (if input separately or a classical code) */
    char *finK; /** `input file` name for Lz=K (not used much) */
    char *finG; /** `input file` name for Hz=G (must use separate input) */
    char *finP; /** `input file` name for P (if input separately or a classical code) */
    char *finC; /** `input file` name for `C` (list of non-trivial CWs for decoding) */
    char *outC; /** `output file` name for `C` (list of non-trivial CWs for decoding) */
    char *fout; /** `output file name`  header for files creaded with `mode=3` */
    char *fdem; /** `input file` name for detector error model (`DEM`) */
    char *fdet; /** `input file` name for detector events */
    char *fobs; /** `input file` name for observables */
    char *ferr; /** `input file` name for error vectors */
    //    char *gdet; /** `output file` name for `generated` detector events */
    //    char *gobs; /** `output file` name for `generated` observables */
    //    char *gerr; /** `output file` name for `generated` error vectors */
    char *pdet; /** `output file` name for `predicted` detector events */
    char *pobs; /** `output file` name for `predicted` observables */
    char *perr; /** `output file` name for `predicted` error vectors */
    int classical; /** `1` if this is a classical code? */
    int internal; /** `1` to generate obs/det internally, `2` to generate from `err` file */
    long long int seed;  /** rng `seed`, set<=0 for automatic */
    double useP; /** global error probability `overriding` values in the `DEM` file (default: 0, no override) */
    double *vP; /** probability vector (total of `n`) */
    qllr_t *vLLR; /** vector of LLRs (total of `n`) */
    int minW; /** minimum weight of a codeword or error vector found */
    int dW; /** weight over `minW` to keep the CW or error vector in a hash (-1: no limit; default `0`) */
    int maxW; /** if non-zero, skip any vectors above this weight (default `0`, no upper limit) */
    qllr_t minE; /** minimum energy of a codeword or error vector found */
    qllr_t dE; /** energy over `minE` to keep the CW or error vector in a hash (default: -1, no limit on `E`) */
    double dEdbl; /** temp value */
    int nzH, nzL; /** count of non-zero entries in `H` and `L` */
    csr_t *mH; /** sparse version of `H`=`Hx` (by rows) */
    csr_t *mHt; /** sparse version of H (by columns) */
    csr_t *mL; /** sparse version of `L`=`Lx` (by rows) */
    csr_t *mK; /** sparse version of `K`=`Lz` (by rows) */
    csr_t *mLt; /** sparse version of `L` (by columns) */
    csr_t *mG; /** sparse version of generator matrix `G=Hz` (by rows) */
    /** rows of `G` orthogonal to rows of both `H` and `L` */
    /** rows of `H` orthogonal to rows of both `G` and `K` */
    /** `rank L` = `rank K` = `k`, the number of encoded qubits.  Any non-zero linear combination of rows of `L` 
	gives a non-zero product with (some) rows of `K` and similarly, any such combination of rows of `K` 
	gives a non-zero product with some rows of `L` (we do not require `L Kt=Identity`) */
    //    int maxJ;	/** memory to initially allocate for local storage */
    qllr_t LLRmin;
    qllr_t LLRmax;
    double epsilon; /** probability to ignore, default `1e-8` */
    one_vec_t *codewords; /** `hash table` with found codewords */
    long long int num_cws; /** `number` of codewords in the `hash` */
    FILE *file_err;
    FILE *file_det;
    FILE *file_obs;
    //    FILE *file_err_g; /** out file, `generated` errors */
    //    FILE *file_det_g;
    //    FILE *file_obs_g;
    FILE *file_err_p; /** out file, `predicted` errors */
    FILE *file_det_p;
    FILE *file_obs_p;
    long long int line_err; /** current line of `file_err` */
    long long int line_det; /** current line of `file_det` */
    long long int line_obs; /** current line of `file_obs` */
    mzd_t *mE;
    mzd_t *mHe;
    mzd_t *mLe;
    mzd_t *mHeT;
    mzd_t *mLeT;
  } params_t;

  extern params_t prm;

  /** @bried helper structure to sort by probabilities (inv by LLR) */
  typedef struct IPPAIR_T {int index; qllr_t llr; } ippair_t; 
  
  /** @brief helper function to sort `ippair_t`
   *  use `qsort(array, len, sizeof(ippair_t), cmp_ippairs);`
   */
  static inline int cmp_ippairs(const void *a, const void *b){
    const qllr_t pa=((ippair_t *) a) -> llr;
    const qllr_t pb=((ippair_t *) b) -> llr;
    if (pa<pb)
      return -1; /** was `+1` with probabilities */
    else if (pa>pb)
      return +1; /** was `-1` */
    return 0;
  }

  /** @brief return permutation = decreasing probabilities (increasing LLR) */
  mzp_t * sort_by_llr(mzp_t *perm, const qllr_t vLLR[], params_t const * const p);  

  /** @brief prepare an ordered pivot-skip list of length `n-rank` */
  mzp_t * do_skip_pivs(const size_t rank, const mzp_t * const pivs);
  
  /** functions defined in `iter_dec.c` ******************************************** */
  void cnt_out(int print_banner);
  void cnt_update(extr_t which, int iteration);
  void out_llr(const char str[], const int num, const qllr_t llr[]);
  
  int syndrome_check(const qllr_t LLR[], const mzd_t * const syndrome,
		     const csr_t * const H,
		     _maybe_unused const params_t * const p);
  
  int do_osd_start(qllr_t * LLR, const mzd_t * const srow,
		   const csr_t * const H, const params_t * const p);
  
  int do_parallel_BP(qllr_t * outLLR, const mzd_t * const srow,
		   const csr_t * const H, const csr_t * const Ht,
		     const qllr_t LLR[], const params_t * const p);

  int do_serialC_BP(qllr_t * outLLR, const mzd_t * const srow,
		   const csr_t * const H, const csr_t * const Ht,
		     const qllr_t LLR[], const params_t * const p);

  int do_serialV_BP(qllr_t * outLLR, const mzd_t * const srow,
		    const csr_t * const H, const csr_t * const Ht,
		    const qllr_t LLR[], const params_t * const p);
  
  /** function defined in `star_poly.c` ********************************************* */
  
  /** @brief replace the DEM (matrices Ht, Lt, and LLR vector) with star-triangle transformed */
  int star_triangle(csr_t * Ht, csr_t * Lt, qllr_t *LLR, const one_vec_t * const codewords,
		  _maybe_unused const long int debug);
  
  /** @brief create K=Lz matrix with minimum weight rows from a list of codewords in hash */
  csr_t * do_K_from_C(const csr_t * const mLt, const one_vec_t * const codewords,
		      const int k, const int n, const int minW, const int maxW,
		      _maybe_unused const int debug);

  /** @brief create G=Hz matrix with minimum weight rows from a list of codewords in hash */
  csr_t * do_G_from_C(const csr_t * const mLt, const one_vec_t * const codewords,
		      const int num_need, const int minW, int maxW,
		      _maybe_unused const int debug);  
  
  csr_t * do_G_matrix(const csr_t * const mHt, const csr_t * const mLt, const qllr_t LLR[], 
		      const int rankG, const int debug);

  /** 
   * @brief The help message. 
   * 
   * @todo: This has to be checked and updated, especially the `debug` options.
   */
#define USAGE                                                           \
  "%s:  vecdec - vectorized decoder and LER estimator\n"                \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Command line arguments are processed in the order given except\n"	\
  "\t for 'mode' and 'debug' (these are processed first).\n"		\
  "\t Supported parameters:\n"						\
  "\t --help\t\t: give this help (also '-h' or just 'help')\n"		\
  "\t mode=[integer] help\t\t: help for specific mode\n"		\
  "\t --morehelp\t: give more help on program conventions\n"		\
  "\t fdem=[string]\t: name of the input file with detector error model\n" \
  "\t finH=[string]\t: file with parity check matrix Hx (mm or alist)\n" \
  "\t finG=[string]\t: file with dual check matrix Hz (mm or alist)\n"	\
  "\t finL=[string]\t: file with logical dual check matrix Lx (mm or alist)\n" \
  "\t finK=[string]\t: file with logical check matrix Lz (mm or alist)\n" \
  "\t finP=[string]\t: input file for probabilities (mm or a column of doubles)\n" \
  "\t finC=[string]\t: input file name for codewords in `nzlist` format\n" \
  "\t\t (space is OK in front of file names to enable shell completion)\n" \
  "\t outC=[string]\t: output file name for codewords in `nzlist` format\n" \
  "\t\t\t (if same as finC, the file will be updated)\n"		\
  "\t maxC=[long long int]\t: max number of codewords to read/write/store\n" \
  "\t epsilon=[double]\t: small probability cutoff (default: 1e-8)\n"	\
  "\t useP=[double]\t: fixed probability value (override values in DEM file)\n"	\
  "\t\t for a quantum code specify 'fdem' OR 'finH' and ( 'finL' OR 'finG' );\n" \
  "\t\t for classical just 'finH' (and optionally the dual matrix 'finL')\n" \
  "\t ferr=[string]\t: input file with error vectors (01 format)\n"	\
  "\t fobs=[string]\t: input file with observables (01 matching lines in fdet)\n" \
  "\t fdet=[string]\t: input file with detector events (01 format)\n"   \
  "\t\t specify either 'ferr' OR a pair of 'ferr' and 'fdet' (or none for internal)\n" \
  "\t perr, pobs, pdet=[string]\t: out file for predicted vectors (01 format)\n" \
  "\t pads=[integer]\t: if non-zero, pad vectors from `fdet` file with zeros (0)\n" \
  "\t fout=[string]\t: header for output file names ('tmp', see 'mode=3')\n" \
  "\t steps=[integer]\t: num of RIS or BP decoding steps (default: 50)\n" \
  "\t lerr =[integer]\t: OSD search level (-1, only implemented with `mode=0`, `1`)\n" \
  "\t maxosd=[integer]\t: max column for OSD2 and above (100)\n"	\
  "\t bpalpha=[float]\t: average LLR scaling coefficient for BP (default 0.5)\n" \
  "\t bpretry=[integer]\t: retry BP up to this many times per syndrome (1)\n" \
  "\t swait=[integer]\t: Gauss steps w/o new errors to stop (0, do not stop)\n" \
  "\t nvec =[integer]\t: max vector size for decoding (default: 1024)\n" \
  "\t\t\t (list size for distance or energy calculations)\n"		\
  "\t ntot =[long long int]\t: total syndromes to generate (default: 1)\n"	\
  "\t nfail=[long long int]\t: total fails to terminate (0, do not terminate)\n" \
  "\t dW=[integer]\t: if 'dW>=0', may keep vectors of weight up to 'minW+dW' (0)\n" \
  "\t maxW=[integer]\t: if non-zero, skip any vectors above this weight (0)\n" \
  "\t dE=[double]\t: if 'dE>=0', may keep vectors of energy up to 'minE+dE'\n" \
  "\t\t\t (default value: -1, no upper limit on energy)\n"		\
  "\t dmin=[integer]\t: terminate distance calculation immediately when\n" \
  "\t\t\t a vector of weight 'W<=dmin' is found, return '-w' (default: 0)\n" \
  "\t seed= [long long int]\t: RNG seed or automatic if <=0 (default: 0)\n"	\
  "\t qllr1=[integer]\t: if 'USE_QLLR' is set, parameter 'd1' (12)\n"	\
  "\t qllr2=[integer]\t: if 'USE_QLLR' is set, parameter 'd2' (300)\n"	\
  "\t qllr3=[integer]\t: if 'USE_QLLR' is set, parameter 'd3' (7)\n"	\
  "\t\t These are used to speed-up LLR calculations, see 'qllr.h'\n"	\
  "\t\t Use 'qllr2=0' for min-sum.\n"					\
  "\t mode=int[.int]\t: operation mode[.submode] (default: 0.0)\n"	\
  "\t\t* 0: use basic vectorized (random information set) decoder\n"	\
  "\t\t* 1: Belief Propagation decoder.\n"		\
  "\t\t* 2: Generate most likely fault vectors, estimate Prob(Fail).\n"  \
  "\t\t* 3: Read in the DEM file and optionally write the corresponding \n" \
  "\t\t\t G, K, H, and L matrices and the probability vector P.\n"	\
  "\t debug=[integer]\t: bitmap for aux information to output (default: 1)\n" \
  "\t\t*   0: clear the entire debug bitmap to 0.\n"                    \
  "\t\t*   1: output misc general info (on by default)\n"		\
  "\t\t*   2: output matrices for verification\n"                       \
  "\t\t\t see the source code for more options\n"			\
  "\t See program documentation for input file syntax.\n"               \
  "\t Multiple 'debug' parameters are XOR combined except for 0.\n"	\
  "\t Use debug=0 as the 1st argument to suppress all debug messages.\n"

#define MORE_HELP							\
  "   Matrices used by %s:\n"						\
  "\t We have a CSS code with binary generator matrices Hx=H, Hz=G,\n" \
  "\t and logical-operator generating matrices Lx=L and Lz=K.  These\n"	\
  "\t matrices have 'n' columns each and satisfy orthogonality properties\n" \
  "\t\t Hx*Hz^T=0, Lx*Hz^T=0, Lx*Lz^T=0, and Lx*Lz^T=I (identity matrix).\n" \
  "\t We are trying to correct binary Z errors, whose probabilities\n"	\
  "\t are given by the n-component double vector P\n"	\
  "\t  (or can be overridden with the command-line parameter 'useP').\n" \
  "\t\t   A detector error model (DEM) file, in the format produced by Stim,\n" \
  "\t contains matrices Hx and Lx, and the error probability vector P.\n" \
  "\t The same code can be obtained by specifying the matrices \n"	\
  "\t independently, via separate files.\n"				\
  "   The dual CSS matrix Hz can be specified instead of Lx.\n"	\
  "\t In such a case, the internal error generator must be used\n"	\
  "\t (an attempt to specify 'fdet' and 'fobs' files will result in an error).\n" \
  "   For a classical code, just give the parity check matrix Hx=H.\n"	\
  "\t In this case G matrix is trivial (has zero rank), and\n"		\
  "\t Lx has all rows of weight '1'.  \n"				\
  "\t Only the internal error generator can be used for classical codes\n" \
  "\t                                                       \n"		\
  "\t Parameter(s) used by all modes:                     \n"		\
  "\t seed=[integer] : when negative or zero, combine provided value\n" \
  "\t\t with 'time(null)' and 'pid()' for more randomness.\n"		\
  "\t                                                       \n"	
  

#define HELP0 /** help for `mode=0` */  \
  " mode=0 : use basic vectorized (random information set) decoder\n"	\
  "\t No 'submode' can be used with this mode. \n"	\
  "\t This decoder is exponentially slow but it is not specific to \n" \
  "\t quantum LDPC codes.  Accuracy and performance are \n" \
  "\t determined by parameters 'steps' (number of RIS rounds) and 'lerr'.\n" \
  "\t Long codes may require exponentially large number of 'steps'.\n" \
  "\t Values 'lerr>1' can be slow for long codes.\n" \
  "\t Specify a single DEM file 'fdem', or 'finH', 'finL', and 'finP'\n"\
  "\t separately (either 'finL' or 'finG' is needed for a quantum code).\n" \
  "\t Use 'useP' to override error probability values in DEM file.   \n" \
  "\t Errors can be generated internally or read from 01 file 'ferr'.\n" \
  "\t Alternatively, files with detector events and observables \n"	\
  "\t can be specified via 'fdet' and 'fobs'. \n"			\
  "\t Long lines in these files may be silently truncated. \n"		\
  "\t Use 'pads=1' to pad lines in 'fdet' file with zeros.\n"	\
  "\t Set 'nfail' and/or 'swait' for early termination.\n"		\
  "\t Total of 'ntot' errors will be read or generated in chunks of 'nvec'.\n" \
  "\t                                                       \n"	

#define HELP1 /** help for `mode=1` */  \
  " mode=1.[submode] : use one of several iterative decoder versions\n"	\
  "\tSubmode bitmap values:\n"						\
  "\t\t\t .1 (bit 0) use regular LLR\n"					\
  "\t\t\t .2 (bit 1) use average LLR - these take precendence\n"	\
  "\t\t\t .4 (bit 2) use serial BP schedule (not parallel)\n"		\
  "\t\t\t .8 (bit 3) with serial, use V-based order (not C-based)\n"	\
  "\t\t\t .16 (bit 4) with serial, randomize node order in each round \n" \
  "\t\t\t     (by default randomize only once per run)\n"		\
  "\t For convenience, 'submode=0' is equivalent to 'submode=3'. \n"	\
  "\t This decoder may experience convergence issues.\n" \
  "\t Accuracy and performance are determined by parameters \n" \
  "\t 'steps' (number of BP rounds), 'lerr' (OSD level, defaul=-1, on OSD).\n" \
  "\t and 'maxosd', the number of columns for OSD in levels 2 and above.\n" \
  "\t Using 'steps' not higher than 50 is recommended.\n"		\
  "\t Use 'bpalpha' to specify how averaging is done (default: 0.5).\n"	\
  "\t Use 'bpretry' to specify how many times to retry BP (default: 1, do not retry) \n" \
  "\t Use 'qllr' parameters to set LLR quantization\n"\
  "\t   or compile with 'VER=\"\"' option for double LLR values. \n" \
  "\t With 'qllr2=0' the Sum-Product algorithm reduces to a 'Min-Sum'.\n" \
  "\t Specify a single DEM file 'fdem', or 'finH', 'finL', and 'finP'\n"\
  "\t separately (either 'finL' or 'finG' is needed for a quantum code).\n" \
  "\t Use 'useP' to override error probability values in DEM file.   \n" \
  "\t Errors can be generated internally or read from 01 file 'ferr'.\n" \
  "\t Alternatively, files with detector events and observables \n"	\
  "\t can be specified via 'fdet' and 'fobs'. \n"			\
  "\t Long lines in these files may be silently truncated. \n"		\
  "\t Use 'pads=1' to pad lines in 'fdet' file with zeros.\n"	\
  "\t Set 'nfail' and/or 'swait' for early termination.\n"		\
  "\t Total of 'ntot' errors will be read or generated in chunks of 'nvec'.\n" \
  "\t                                                       \n"	

#define HELP2 /** help for `mode=2` */  \
  " mode=2 : Generate most likely fault vectors, estimate Prob(Fail).\n" \
  "\t (no submode is currently used)\n"					\
  "\t Use up to 'steps' random information set (RIS) steps\n"		\
  "\t unless no new codewords (fault vectors) have been found for 'swait' steps.\n" \
  "\t Use 'steps=0' to just use the codewords from the file \n"		\
  "\t Keep vectors of weight up to 'dW' above min weight found.\n"	\
  "\t       and energy up to 'dE' above minimum E found (sum of LLRs).\n" \
  "\t When 'maxC' is non-zero, generate up to 'maxC' unique codewords.\n" \
  "\t If 'outC' is set, write full list of CWs to this file.\n"		\
  "\t If 'finC' is set, read initial set of CWs from this file.\n"	\
  "\t Accuracy and performance are determined by parameters \n"		\
  "\t 'steps' (number of BP rounds), 'lerr' (OSD level, defaul=-1, on OSD).\n" \
  "\t and 'maxosd', the number of columns for OSD in levels 2 and above.\n" \
  "\t Specify a single DEM file 'fdem', or 'finH', 'finL', and 'finP'\n" \
  "\t separately (either 'finL' or 'finG' is needed for a quantum code).\n" \
  "\t Use 'useP' to override error probability values in DEM file.   \n" \
  "\t                                                       \n"	

#define HELP3 /** help for `mode=3` */  \
  " mode=3 : Export matrices associated with the code.\n" \
  "\t Read in the DEM file and optionally write the corresponding \n"	\
  "\t G, K, H, and L matrices and the probability vector P.\n"	\
  "\t By default (submode&31=0) output everything, otherwise\n"	\
  "\t .1 (bit 0) write G=Hz matrix\n"				\
  "\t .2 (bit 1) write K=Lz matrix\n"				\
  "\t .4 (bit 2) write H=Hx matrix\n"				\
  "\t .8 (bit 3) write L=Lx matrix\n"					\
  "\t .16 (bit 4) write P vector\n"					\
  "\t Codewords file 'finC', if given, will be used to create 'G' and 'K'\n" \
  "\t   matrices with rows of smallest possible weight.\n"		\
  "\t In addition, mode=3.32 (just bit 5 set) in combination with\n"	\
  "\t  codewords file 'finC' forces code transformation mode.\n"	\
  "\t Similarly, use mode=3.64 (just bit 6 set) to create DEM file '${fout}D.dem'\n" \
  "\t Use 'fout=' command line argument to generate file names\n"	\
  "\t ${fout}H.mmx, ${fout}G.mmx, ${fout}L.mmx, ${fout}K.mmx, and ${fout}P.mmx\n" \
  "\t with 'fout=stdout' all output is sent to 'stdout'\n"		\

  
#ifdef __cplusplus
}
#endif

#endif /* VECDEC_H */
