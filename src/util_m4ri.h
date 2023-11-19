#ifndef UTIL_M4RI_H
#define UTIL_M4RI_H

/************************************************************************ 
 * helper functions for use with m4ri library, including binary sparse
 * matrices and conversion utilities 
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu> 
 * some code borrowed from various sources 
 ************************************************************************/

#define SWAPINT(a,b) do{ int t=a; a=b; b=t; } while(0)


/**
 * macros from nauty.h
 * SETWD(pos) gives the setword in which pos is located
 * SETBT(pos) gives the location of bit pos in a setword
 */
#define SETWD(pos) ((pos)>>6)
#define SETBT(pos) ((pos)&0x3F)
#define TIMESWORDSIZE(w) ((w)<<6)    /* w*WORDSIZE */

#define FIRSTBIT(x) __builtin_ctzll(x) // number of trailing zeros 

// #ifdef __POPCNT__
#if 0

/* 
 * optimized code copied verbatim from https://danluu.com/assembly-intrinsics/ 
 */
uint32_t builtin_popcnt_unrolled_errata_manual(const uint64_t* buf, int len) {
  assert(len % 4 == 0);
  uint64_t cnt[4];
  for (int i = 0; i < 4; ++i) {
    cnt[i] = 0;
  }

  for (int i = 0; i < len; i+=4) {
    __asm__(
        "popcnt %4, %4  \n\t"
        "add %4, %0     \n\t"
        "popcnt %5, %5  \n\t"
        "add %5, %1     \n\t"
        "popcnt %6, %6  \n\t"
        "add %6, %2     \n\t"
        "popcnt %7, %7  \n\t"
        "add %7, %3     \n\t" // +r means input/output, r means intput
        : "+r" (cnt[0]), "+r" (cnt[1]), "+r" (cnt[2]), "+r" (cnt[3])
        : "r"  (buf[i]), "r"  (buf[i+1]), "r"  (buf[i+2]), "r"  (buf[i+3]));
  }
  return cnt[0] + cnt[1] + cnt[2] + cnt[3];
}

static inline int m4ri_bitcount(word w){
  return __builtin_popcountll(w);  
}

#else /* no __POPCNT__ */

#define MASK(c)    (((uint64_t)(-1)) / (__M4RI_TWOPOW(__M4RI_TWOPOW(c)) + 1))
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (__M4RI_TWOPOW(c))) & MASK(c))

static inline int m4ri_bitcount(word w)  {
   uint64_t n = __M4RI_CONVERT_TO_UINT64_T(w);
   n = COUNT(n, 0);
   n = COUNT(n, 1);
   n = COUNT(n, 2);
   n = COUNT(n, 3);
   n = COUNT(n, 4);
   n = COUNT(n, 5);
   return (int)n;
}

#endif /* __POPCNT__ */


/**
 * sparse binary matrix in compressed-row form (CSR, nz=-1) or 
 * List-Of-Pairs (nz pairs).
 * use mzp_compress() to convert from LOP to CSR. 
 */
typedef struct{    /*  */
  int rows ;	    /* number of rows */
  int cols ;	    /* number of columns */
  int nz ;	    /* # of entries in triplet matrix */
  int nzmax ;	    /* # allocated size */
  int *p ;	    /* row pointers (size rows+1) OR row indices */
  int *i ;	    /* col indices, size nzmax */
} csr_t ;


typedef struct { int a; int b; } int_pair;


#if defined(__cplusplus) && !defined (_MSC_VER)
extern "C" {
#endif


/**
 * @brief initialize a CSR matrix from list of pairs (row,col)
 * @param mat use existing handle if sufficient size (allocate if NULL)
 * @param nz number of non-zero pairs in the list 
 * @param prs array of pairs (will be sorted)
 * @param nrows matrix dimension
 * @param ncols matrix dimension
 * @return the created matrix 
 */
csr_t * csr_from_pairs(csr_t *mat, int nz, int_pair * const prs, int nrows, int ncols);

/** 
 * @brief convert `m4ri` dense matrix to `csr`
 * 
 * Optimized for sparse matrices.  
 * 
 * @param mat use existing handle if sufficient size (allocate if NULL)
 * @param orig the existing matrix to be converted
 * @return the constructed matrix 
 *
 */
csr_t * csr_from_mzd(csr_t *mat, const mzd_t * const orig);

/**
 * @brief Compute logical generator matrix Lx for a CSS code
 *
 * Given a pair of binary CSS generator matrices in CSR format,
 * Hx*Hz^T=0, compute a sparse matrix Lx s.t. Lx*Hz^T=0
 * and rows of Lx be linearly independent from those of Hx.
 * TODO: see if sparsity of Lx can be improved.
 */

  csr_t * Lx_for_CSS_code(const csr_t * const Hx, const csr_t *const Hz);

/** 
 * number of set bits in a matrix.  
 * TODO: what is the problem with built-in version?  
 */
size_t mzd_weight(const mzd_t *A);

/**
* nextelement(set1,m,pos) = the position of the first element in set set1   
* which occupies a position greater or equal than pos.  If no such element exists,   
* the value is -1.  pos can have any value less than n, including negative  
* values.                                                                   
*  
* near verbatim copy from naututil.c (Nauty library by Brendan McKay)
*/
//  int nextelement(word *set1, int m, int pos);

/** 
 * return first non-zero bit in raw set-word vector set1 
 * of length m, starting with position pos.
 * with all zero bits, return -1 or number outside the range
 */
static inline int nextelement(const word * const set1, const int m, const int pos){
  word setwd;
  int w;
  w = SETWD(pos);
  setwd = set1[w] & (m4ri_ffff<< SETBT(pos));
  for (;;){
    if (setwd != 0) return  TIMESWORDSIZE(w) + FIRSTBIT(setwd);
    if (++w == m) return -1;
    setwd = set1[w];
  }
}


/**
 * Copy of mzd_gauss_delayed from mzd.c (m4ri package) except additionally 
 * returns the list of pivot columns in second argument 
 */
rci_t mzd_gauss_naive(mzd_t *M, mzp_t *q, int full);

/** 
 * return max row weight of CSR matrix p
 * TODO: add code for List of Pairs 
 */
  int csr_max_row_wght(const csr_t *p);
  
/** 
 * transpose compressed CSR matrix, 
 * (re) allocate if needed 
 * return resulting matrix
 * TODO: add code for List of Pairs 
 */
  csr_t * csr_transpose(csr_t *dst, const csr_t *p);

  
/**
 * Convert CSR sparse binary matrix to MZD
 * allocate dst if needed (must be correct size or NULL)
 */
mzd_t *mzd_from_csr(mzd_t *dst, const csr_t *p);

/**
 * Convert a sparse binary matrix CSR into a standard form [ I C ],
 * with some col permutations if needed, create the dense generator
 * [ CT I ], and permute the cols back.  
 * (re)allocate G if needed.
 */
mzd_t *mzd_generator_from_csr(mzd_t *G, const csr_t * const H);

/**
 * sparse-S by dense B multiplication
 * C=C+S*B; allocate C if needed.
 */
mzd_t * csr_mzd_mul(mzd_t *C, const csr_t *S, const mzd_t *B, int clear);

/**
 * helper function to compute the weight of the product 
 * A*B (transpose == 0) or A*B^T (transpose == 1)
 * with A sparse, B dense binary matrices
 */
size_t product_weight_csr_mzd(const csr_t *A, const mzd_t *B, int transpose);

/**
 * return uniformly distributed random number in the range [0,...,max-1] 
 * uses RAND internally 
 * \todo Replace by a better generator 
 */
int rand_uniform(const int max);

/**
 * @brief replace pivot q with a random pivot of same length, 
 * *** WARNING: LAPACK style pivot permutations! ***
 * `pivots=[p0,p1,p2,...]` with `p0<p1<p2...`
 * requires pair permutations `(0,p0),(1,p1),(2,p2), ...`
 * @return pointer to q.
 * input: perm -- existing permutation
 */ 
mzp_t * mzp_rand(mzp_t *q);

/**
 * print out the permutation (only needed under windows)
 */
void mzp_out(mzp_t const *p);

/**
 * apply pivot p to permutation q in place from start; 
 * initialize q to identity permutation if NULL
 * return q 
 */
mzp_t *perm_p(mzp_t *q, const mzp_t *p,rci_t start);

/**
 * apply pivot p (transposed) to permutation q in place from start; 
 * initialize q to identity permutation if NULL
 * return q 
 */
mzp_t *perm_p_trans(mzp_t *q, const mzp_t *p,const rci_t start);


/**
 * kill a CSR matrix 
 */
csr_t *csr_free(csr_t *p);

/**
 * initialize a CSR matrix 
 * check existing size and (re)allocate if  needded 
 */
csr_t *csr_init(csr_t *mat, int rows, int cols, int nzmax);

/**
 *  compress a CSR matrix  
 */ 
void csr_compress(csr_t *mat);

/**
 *  output a CSR matrix  
 */ 
  void csr_out(const csr_t *mat);

/**
 * read sparse matrix into a (binary) CSR (all entries default to 1)
 * (re)allocate mat if needed
 * use transpose=1 to transpose.
 */
  csr_t *csr_mm_read(char *fin, csr_t *mat, int transpose);

  /** write out CSR matrix.  Use fout="stdout" for `stdout` output */
void csr_mm_write( char * const fout, const char fext[], const csr_t * const mat,
		  const char comment[]);
  
/** 
 * Permute columns of a CSR matrix with permutation perm.
 */
csr_t *csr_apply_perm(csr_t *dst, const csr_t * const src, const mzp_t * const perm);

/**
 * \brief Flip the bit at position M[row,col].
 *
 * \param M Matrix
 * \param row Row index
 * \param col Column index
 *
 * \note No bounds checks whatsoever are performed.
 *
 */

static inline void mzd_flip_bit(mzd_t const *M, rci_t const row, rci_t const col ) {
  __M4RI_FLIP_BIT(M->rows[row][col/m4ri_radix], col%m4ri_radix);
}

  
/** 
 * Check whether syndrome is zero or not 
 */ 
int syndrome_bit_count(mzd_t *row, csr_t *spaQ);

/** 
 * Check if row is linearly dependent with the rows of matP0
 * which is assumed to be in standard form.
 * rankP0 is the number of non-zero rows in rankP0.
 */

int do_reduce(mzd_t *row, const mzd_t *matP0, const rci_t rankP0);

/**
 * generate binary error vector with error probability p 
 */
void make_err(mzd_t *row, double p);

#if defined(__cplusplus) && !defined (_MSC_VER)
}
#endif
  
#endif /* UTIL_M4RI_H */
