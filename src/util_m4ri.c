/************************************************************************ 
 * helper functions for use with m4ri library, 
 * including binary sparse matrices and conversion utilities 
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu> 
 * some code borrowed from various sources 
 ************************************************************************/
#include <m4ri/m4ri.h>
#include <m4ri/mzd.h>
#include <stdio.h>
// #include <"copy_m4ri.h"
#include "mmio.h"
#include <ctype.h>

#include "utils.h"
#include "util_m4ri.h"


size_t mzd_weight_naive(const mzd_t *A){
  size_t count = 0;
  for(rci_t i = 0; i < A->nrows; ++i)
    for(rci_t j = 0; j < A->ncols; ++j)
      if(mzd_read_bit(A, i, j))
	++count;
  return (count);
}

void mzd_row_print_sparse(const mzd_t * const A, const int row){
  //  for(rci_t i = 0; i < A->nrows; ++i){
  int i = row;
    printf("# row=%d [",i);
    for(rci_t j = 0; j < A->ncols; ++j)
      if(mzd_read_bit(A, i, j))
	printf(" %d",j);
    printf(" ]\n");
    //  }
}


size_t mzd_weight(const mzd_t * const A){
  size_t count = 0;
  if(A->width == 1) {
    for(rci_t i = 0; i < A->nrows; ++i)
      for(rci_t j = 0; j < A->ncols; ++j)
	if(mzd_read_bit(A, i, j))
	  ++count;
    return (count);
  }
  for(rci_t i = 0; i < A->nrows; ++i) {
    const word * const truerow = mzd_row_const(A,i);
    for(wi_t j = 0; j < A->width - 1; j ++)
      count += m4ri_bitcount(truerow[j]);

    for(int j = 0; j < A->ncols % m4ri_radix; ++j)
      if(mzd_read_bit(A, i, m4ri_radix * (A->ncols / m4ri_radix) + j))
	++count;
  }

  return count ;
}

int mzd_row_is_zero(const mzd_t  * const A, const int i) {
  const word mask_end = A->high_bitmask;
  const word * const truerow = mzd_row_const(A,i);
  for (wi_t j = 0; j < A->width - 1; ++j)
    if(truerow[j])
      return 0;
  if(truerow[A->width - 1] & mask_end)
    return 0;
  return 1;
}



/**
 * Copy of mzd_gauss_delayed from mzd.c (m4ri package) except additionally 
 * returns the list of pivot columns in second argument 
 * permutation q must be initialized!
 */
rci_t mzd_gauss_naive(mzd_t *M, mzp_t *q, int full) {
  rci_t startcol=0, startrow = 0;
  rci_t pivots = 0;
#if 0  
  if (q==NULL) 
    q=mzp_init(M->ncols);
  else {
    if ((q->length)!=M->ncols){ 
      mzp_free(q);
      q=mzp_init(M->ncols);
    }
    else 
      mzp_set_ui(q,1);
  }
#endif 

  for (rci_t i = startcol; i < M->ncols ; ++i) {
    for(rci_t j = startrow ; j < M->nrows; ++j) {
      if (mzd_read_bit(M, j, i)) {
	mzd_row_swap(M, startrow, j);
	q->values[pivots]=i; /* request cols swap: pivots <--> i */
	++pivots;
	for(rci_t ii = full ? 0 : startrow + 1;  ii < M->nrows; ++ii) {
	  if (ii != startrow) {
	    if (mzd_read_bit(M, ii, i)) {
	      mzd_row_add_offset(M, ii, startrow, i);
	    }
	  }
	}
	startrow = startrow + 1;
	break;
      }
    }
  }
  return pivots;
}

/** 
 * return max row weight of CSR matrix p
 * TODO: add code for List of Pairs 
 */
int csr_max_row_wght(const csr_t *p){
  int m=p->rows;
  int wmax=0;
  for(int i=0;i<m;i++){
    int w=(p->p[i+1])-(p->p[i]);
    if (w>wmax)
      wmax=w;
    }
  return wmax;
}

long long int csr_min_max_blk(csr_t *mat, int r1, int r2){
  if((r1 > r2) || (r2 >= mat->rows))
    ERROR("r1=%d r2=%d; must be a non-empty range of valid matrixs rows < %d\n",r1,r2,mat->rows);
  long long int max=-1, min=mat->cols;
  for(int j=r1; j<= r2; j++){
    if (mat->p[j+1] > mat->p[j] ){
      int p1 = mat->p[j];      /* index of 1st column in this row */
      int p2 = mat->p[j+1] - 1;/* index of last col in this row */
      if(max < mat->i[p2])
	max = mat->i[p2];
      if(min > mat->i[p1])
	min = mat->i[p1];
    }
  }
  if(max<0)
    ERROR("specified row interval r1=%d r2=%d is empty \n",r1,r2);
  return min + (mat->cols) * max;
}

/** 
 * transpose compressed CSR matrix, 
 * (re) allocate if needed 
 * return resulting matrix
 * TODO: add code for List of Pairs 
 */
csr_t * csr_transpose(csr_t *dst, const csr_t * const org){
  int rows=org->rows, cols=org->cols, nz=org->p[rows];
  if (dst == NULL) 
    dst = csr_init(NULL,cols,rows,nz);
  else if ((dst->cols != rows) || (dst->rows != cols) || (dst->nzmax < MAX(nz, cols+1))) 
    ERROR("Wrong size for return matrix.\n");
  else
    dst->nz=0; /* clear matrix */
  for(int i=0;i<rows;i++)    
    for(int j=org->p[i]; j < org->p[i+1] ; j++){
      dst->p[j]=org->i[j]; // pair format, to be compressed later
      dst->i[j]=i ;
    }
  dst->nz=nz;
  csr_compress(dst);
  return dst;
}

/** 
 * @brief create a submatrix of a CSR matrix
 * return resulting matrix
 * TODO: add more detailed error messages
 * TODO: add code for List of Pairs 
 */

csr_t * csr_submatrix(const csr_t * const p, int minR, int maxR, int minC, int maxC){
  const int rows=p->rows, cols=p->cols;
  if(p->nz != -1)
    ERROR("nz=%d, List of Pairs format not supported\n",p->nz);
  if((minR<0) || (minR >= maxR) || (maxR > rows)||
     (minC<0) || (minC >= maxC) || (maxC > cols))
    ERROR("invalid parameters specified");	  
  int new_nz=0;
  for(int j=minR; j< maxR; j++){
    for(int i = p->p[j]; i < p->p[j+1]; i++){
      int col = p->i[i];
      if ((col >= minC) && (col < maxC))
	new_nz++;
    }
  }
  csr_t *submat = csr_init(NULL,maxR-minR,maxC-minC, new_nz);
  
  int idx=0, j, j0;
  for(j0=minR, j=0; j0< maxR; j0++, j++){
    submat->p[j]=idx; 
    for(int i = p->p[j0]; i < p->p[j0+1]; i++){
      int col0 = p->i[i];
      if ((col0 >= minC) && (col0 < maxC))
	submat->i[idx++] = col0 - minC;
    }
  }
  submat->p[j]=idx;
  submat->nz=-1;
  
  return submat;
}


/**
 * Convert CSR sparse binary matrix to MZD
 * allocate dst if needed (must be correct size or NULL)
 */
mzd_t *mzd_from_csr(mzd_t *dst, const csr_t *p) {
  int m=p->rows, n=p->cols, nz=p->nz;
  if (dst == NULL) 
    dst = mzd_init(m,n);
  else if ((dst->nrows != m) || (dst->ncols != n)) 
    ERROR("Wrong size for return matrix.\n");
  else
    mzd_set_ui(dst,0); /* clear bits */

  if(nz==-1){ // binary CSR matrix in Compressed Row Form
    for(int i=0;i<m;i++)
      for(int j=p->p[i]; j < p->p[i+1] ; j++)
	mzd_write_bit(dst, i, p->i[j], 1);
  }
  else{  // binary CSR matrix in Pair Form
    for(int i=0;i<p->nz;i++)
      mzd_write_bit(dst, p->p[i],p->i[i],1);
  }
  return dst;
}

/**
 * Convert a sparse binary matrix CSR into a standard form [ I C ],
 * with some col permutations if needed, create the dense generator
 * [ CT I ], and permute the cols back.  
 * (re)allocate G if needed.
 */
mzd_t *mzd_generator_from_csr(mzd_t *G, const csr_t * const H){  
  rci_t n=H->cols;
  mzd_t *mat=mzd_from_csr(NULL,H); /* convert to dense matrix */
  mzp_t * pivots = mzp_init(n);    /* initialize the permutation */
  rci_t ra = mzd_gauss_naive(mat, pivots, 1);
#if 0  
  if ((prm.debug & 2048) && (n<=80)){ 
    mzd_print(mat);
    mzp_out(pivots);    
    printf("\n");
  }
#endif
  /* now we can check if G is the right size and reallocate if needed */
  if ((G!=NULL) && ((G->ncols!=n) || (G->nrows!=n-ra))){
    mzd_free(G); 
    G=NULL;
  }
  if (G==NULL)
    G=mzd_init(n-ra,n);    
  else
    mzd_set_ui(G,0); /* clear the matrix */
  mzd_apply_p_right_trans(mat,pivots); /* permute columns to make std form [ I C ] */
  mzd_t *matC = mzd_submatrix(NULL, mat, 0, ra, ra, n); /* C size ra by n-ra */
  mzd_t * winCT=mzd_init_window(G,0,0,n-ra,ra);
  mzd_transpose(winCT,matC);
  mzd_free(winCT); /* window is no longer needed */
  mzd_free(matC);  /* \todo see why window does not seem to work here */
  for(int i=0; i<n-ra ; i++)
    mzd_write_bit(G,i,i+ra,1); /* identity matrix on the right */

  mzd_apply_p_right(G,pivots); /* permute columns back */
  mzp_free(pivots);
  mzd_free(mat);
  return G;
}


/**
 * sparse-S by dense B multiplication
 * C=C+S*B; allocate C if needed.
 */
mzd_t * csr_mzd_mul(mzd_t *C, const csr_t *S, const mzd_t *B, int clear){
  if (S->cols != B->nrows)
    ERROR("column dim of S and row dim of B should match, %d=%d",S->cols,B->nrows);
  if (C == NULL) 
    C = mzd_init(S->rows, B->ncols);
  else {
    if (C->nrows != S->rows || C->ncols != B->ncols) 
      ERROR("Provided return matrix has wrong dimensions.\n");
    if(clear)
      mzd_set_ui(C,0); 
  }
  rci_t const m = S->rows;
  for(rci_t i = 0; i < m; ++i)
    for(int j=S->p[i]; j < S->p[i+1] ; j++)
      mzd_combine(C,i,0, C,i,0, B,S->i[j],0); /* combine the two rows */
  return C;
}

#if 0
/** @brief calculate `C = C + A*B`, optimized to the case where `A` is a single row */
mzd_t * mzd_csr_mul(mzd_t *C, const mzd_t *A, const csr_t *B, int clear){
  /** use mzd_row_csr_mul() */
  return NULL;  
}
#endif 


/** 
 * Check whether syndrome is zero or not 
 */ 
int syndrome_bit_count(mzd_t *row, csr_t *spaQ){
  int wei=0;
  int m=spaQ->rows;
  //  mzd_print(row);
  for(int i=0;i<m;i++){
    rci_t bit=0;
    for (int j=spaQ->p[i]; j<spaQ->p[i+1]; j++){
      //      printf("##### here ######## m=%d j=%d %d\n",m,j,spaQ->i[j]);      
      /* printf("checking %d %d: %d\n",i,spaQ->i[j],   */
      /* 	     mzd_read_bit(row, 0,spaQ->i[j]) ? 1 : 0);   */
      bit ^= mzd_read_bit(row, 0,spaQ->i[j]); 
    }
    wei+=bit;
    //    printf("wei=%d i=%d  bit=%d \n",wei,i,bit);
  }
  return wei;
}

/**
 * helper function to compute the weight of the product 
 * A*B (transpose == 0) or A*B^T (transpose == 1)
 * with A sparse, B dense binary matrices
 */
size_t product_weight_csr_mzd(const csr_t *A, const mzd_t *B, int transpose){
  mzd_t *Prod;
  mzd_t *mA=mzd_from_csr(NULL,A);
  if (transpose){
    mzd_t *BT=mzd_transpose(NULL,B);
    Prod=mzd_mul_naive(NULL,mA,BT);
    mzd_free(BT);
  }
  else 
    Prod=mzd_mul_naive(NULL,mA,B);
  mzd_free(mA);
  //  mzd_print(Prod); printf("\n");
  size_t wei=mzd_weight(Prod);
  mzd_free(Prod);
  return wei;
}

/**
 * return uniformly distributed random number in the range [0,...,max-1] 
 */
int rand_uniform(const int max){
#if 0  
  int divisor = RAND_MAX/(max);
  int retval;
  do  
    retval = rand() / divisor;
  while (retval >= max);
  return retval;
#else
  return (int) floor(max * tinymt64_generate_double(&tinymt));
#endif /* 0 */
}

/**
 * replace pivot q with a random pivot of same length, 
 * *** note: LAPACK style pivot permutations! ***
 * return pointer to q.
 * input: perm -- existing permutation
 */ 
mzp_t * mzp_rand(mzp_t *q){
  if (q==NULL)
    ERROR("permutation must be initialized!");
  rci_t length=q->length;
  for(int i=0;i<=length-2;i++){
    q->values[i]= i+rand_uniform(length-i);
    //    printf("%d %d %d\t",i,q->values[i],length-i);
  }
  q->values[length-1]=length-1; /* no swap for the last column */
  //  printf("\n");
  return q;
}

/**
 * print out the permutation (only needed under windows)
 */
void mzp_out(mzp_t const *p){
  printf("[");
  for(rci_t i=0;i<p->length;i++)
    printf(" %d",p->values[i]);
  printf(" ]\n");
}

/**
 * apply pivot p to permutation q in place from start; 
 * initialize q to identity permutation if NULL
 * return q 
 */
mzp_t *perm_p(mzp_t *q, const mzp_t *p,rci_t start){
  assert(p!=NULL);
  if (q==NULL)
    q=mzp_init(p->length);
  else if  (p->length != q->length)
    ERROR("mzp length mismatch %d and %d !\n", p->length , q->length);
  rci_t length = p->length;
  rci_t *perm = q->values; 
  for(rci_t i = start; i < length; ++i) 
    SWAPINT(perm[i], perm[p->values[i]]);
  return q;
}

/**
 * apply pivot p (transposed) to permutation q in place from start; 
 * initialize q to identity permutation if NULL
 * return q 
 */
mzp_t *perm_p_trans(mzp_t *q, const mzp_t *p,const rci_t start){
  assert(p!=NULL);
  if (q==NULL)
    q=mzp_init(p->length);
  else if  (p->length != q->length)
    ERROR("mzp length mismatch %d and %d !\n", p->length , q->length);
  rci_t length = p->length;  
  rci_t *perm = q->values; 
  for(rci_t i = start; i < length; ++i)
    SWAPINT(perm[length-i-1], perm[p->values[length-i-1]]);
  return q;
}


/**
 * kill a CSR matrix 
 */
csr_t *csr_free(csr_t *p){
  if(p!=NULL){
    free(p->i);
    free(p->p);
    p->nzmax=0;
    p->nz=0;
    p->rows=0;
    p->cols=0;
    free(p);
  }
  return NULL;
}

/**
 * initialize a CSR matrix 
 * check existing size and (re)allocate if  needded 
 */
csr_t *csr_init(csr_t *mat, int rows, int cols, int nzmax){
  const int max = MAX(nzmax,(rows+1));
  if((nzmax<0)||(rows<0)||(cols<0))
    ERROR("invalid parameters rows=%d cols=%d nzmax=%d)\n",rows,cols,nzmax);
  if ((mat!=NULL)&&((mat->nzmax < max)||(mat->nzmax < rows+1))){
    /* allocated size was too small */
    /** keep allocated `mat` */
    mat->p = realloc(mat->p, max*sizeof(int));
    mat->i = realloc(mat->i, nzmax*sizeof(int));
    if ((mat->p==NULL) || (mat->i==NULL))
      ERROR("csr_init: failed to reallocate CSR rows=%d cols=%d nzmax=%d",
            rows,cols,nzmax);
    mat->nzmax=nzmax;
    mat->p[0]=0;
  }
  else if(mat==NULL){
    mat=malloc(sizeof(csr_t));  
    mat->p=calloc(max,sizeof(int));
    mat->i=calloc(nzmax,sizeof(int)); 
    if ((mat == NULL) || (mat->p==NULL) || (mat->i==NULL))
      ERROR("csr_init: failed to allocate CSR rows=%d cols=%d nzmax=%d",
            rows,cols,nzmax);
    mat->nzmax=nzmax;
  }
  mat->rows=rows;
  mat->cols=cols;
  mat->nz=0; /* empty */
  return mat;
}

/** 
 * @brief return sparse identity matrix 
 */
csr_t *csr_identity(int rows, int cols){
  const int nzmax= rows < cols ? rows : cols ;
  csr_t * mat = csr_init(NULL, rows, cols, nzmax);
  int i;
  for(i=0; i < nzmax; i++){
    mat->p[i]=i;
    mat->i[i]=i;    
  }
  for(int j=i ; j<=rows; j++) /** including the final value */
    mat->p[j]=i; 
  mat->nz=-1; /* indicate compressed form */
  return mat;
}

/* helper function */
static int cmp_int_pairs(const void *p1, const void *p2){
  if ((((int_pair *) p1)->a)!=(((int_pair *) p2)->a))
    return  ((((int_pair *) p1)->a)-(((int_pair *) p2)->a));    
  return ((((int_pair *) p1)->b)-(((int_pair *) p2)->b));
}

/**
 * @brief initialize a CSR matrix from list of pairs (row,col)
 * @param mat use existing handle if sufficient size (allocate if NULL)
 * @param nz number of non-zero pairs in the list 
 * @param prs array of pairs (will be sorted)
 * @param nrows matrix dimension
 * @param ncols matrix dimension
 */
csr_t * csr_from_pairs(csr_t *mat, const int nz, int_pair * const prs, const int nrows, const int ncols){
  mat = csr_init(mat, nrows, ncols, nz);
  qsort(prs, nz, sizeof(int_pair), cmp_int_pairs);
  int i, j=0;
  for(i=0; i < nrows; i++){
    mat->p[i]=j;
    while ((j<nz) && (prs[j].a == i)){
      mat->i[j]=prs[j].b;
      j++;
    }
  }
  mat->p[i]=j; /* final value */
  mat->nz=-1; /* indicate compressed form */
  return mat;
}


/** 
 * @brief convert `m4ri` dense matrix to `csr`
 * 
 * Optimized for sparse matrices.  
 * 
 * @param mat pointer to existing matrix or NULL 
 *
 */
csr_t * csr_from_mzd(csr_t *mat, const mzd_t * const orig){
  int nz=mzd_weight_naive(orig);
  mat = csr_init(mat, orig->nrows, orig->ncols, nz);/** will reallocate if needed */
  int i, j=0;
  for(i=0;i < mat->rows; i++){
    mat->p[i]=j;
#if 1 /** optimized version */
    int idx=0;
    const word * const rawrow = mzd_row_const(orig,i);
    while(((idx=nextelement(rawrow,orig->width,idx))!=-1)&&
	  (idx>=0)&&
	  (idx<orig->ncols)){
      if(j>=nz)
	ERROR("j=%d nz=%d unexpected",j,nz);
      if(idx>=0)
	mat->i[j++]=idx++;
      if(idx >= orig->ncols)
	break;      
    }
#else /** naive version */
    for(int idx=0; idx< orig->ncols; idx++)
      if(mzd_read_bit(orig,i,idx)){
	mat->i[j++]=idx;
	//	printf("i=%d j=%d idx=%d\n",i,j,idx);
      }
#endif /* 0 */    
  }      
  mat->p[i]=j; /* final value */
  assert(j==nz);
  mat->nz=-1; /* indicate compressed form */
  return mat;
}

/**
 * @brief Compute logical generator matrix Lx for a CSS code
 *
 * Given a pair of binary CSS generator matrices in CSR format,
 * `Hx*Hz^T=0`, compute a sparse matrix `Lx` s.t. `Lx*Hz^T=0`
 * and rows of `Lx` be linearly independent from those of `Hx`.
 * TODO: see if sparsity of Lx can be improved.
 */

csr_t * Lx_for_CSS_code(const csr_t * const Hx, const csr_t *const Hz){
  rci_t n=Hx->cols;
  if(n!=Hz->cols)
    ERROR("unequal number of columns in matrices Hx[%d,%d] and Hz[%d,%d]",
	  Hx->rows, Hx->cols, Hz->rows, Hz->cols);  
  mzd_t *Mx = mzd_from_csr(NULL,Hx); /* convert to dense matrix */
  mzp_t *pivots = mzp_init(n);     /* initialize the permutation */
  rci_t rank = mzd_gauss_naive(Mx, pivots, 1);
  mzd_apply_p_right_trans(Mx,pivots); /** permute columns to make std form [ I C ] */
  mzp_t *perm  = perm_p(      NULL, pivots,0);
  mzp_t *permT = perm_p_trans(NULL, pivots,0);
  mzp_free(pivots);
  csr_t* HzPerm=csr_apply_perm(NULL,Hz,permT);
  mzp_free(permT);
#ifndef NDEBUG /** `verify` orthogonality */
  mzd_t *MxT = mzd_transpose(NULL,Mx); 
  if(product_weight_csr_mzd(HzPerm,MxT,0))
    ERROR("rows of Hx and Hz should be orthogonal \n");
  mzd_free(MxT);
#endif
  mzd_t *MzStar = mzd_generator_from_csr(NULL, HzPerm);
  csr_free(HzPerm);

  mzd_t *mat = mzd_stack(NULL,Mx,MzStar);
  mzd_free(MzStar);
  //  printf("stacked:\n");  mzd_print(mat);
  int rank1=mzd_echelonize(mat,0);
  //  printf("after gauss %d -> %d:\n", rank, rank1); mzd_print(mat);

  mzd_t *window = mzd_init_window(mat,rank,0,rank1,n);
  csr_t *Lx=csr_from_mzd(NULL, window);
  csr_t *ans = csr_apply_perm(NULL,Lx,perm);
  csr_free(Lx);  
  mzd_free(window);
  mzp_free(perm);

  return ans;
}
  


/**
 *  compress a CSR matrix  
 */ 
void csr_compress(csr_t *mat){
  int nz=mat->nz;
  if(nz==-1)
    ERROR("matrix already compressed");
  int_pair *pairs=calloc(nz,sizeof(int_pair));
  for(int i=0;i<nz;i++){
    pairs[i].a=mat->p[i];
    pairs[i].b=mat->i[i];
  }
  qsort(pairs,nz,sizeof(int_pair),cmp_int_pairs);
  int i, j=0;
  for(i=0;i<mat->rows;i++){
    mat->p[i]=j;
    while ((j<nz)&&(pairs[j].a == i)){
      mat->i[j]=pairs[j].b;
      j++;
    }
  }
  mat->p[i]=j; /* final value */
  mat->nz=-1; /* indicate compressed form */
  free(pairs);
}

/**
 *  output a CSR matrix  
 */ 
void csr_out(const csr_t *mat){
  int nz=mat->nz;
  if(nz==-1){
    printf("# binary CSR matrix (%d x %d) in Compressed Row Form: "
	   "nzmax= %d, nz=%d\n",
	   mat->rows,mat->cols,mat->nzmax,mat->p[mat->rows]);
    for(int i=0;i<mat->rows;i++){
      printf("%d:",i+1);
      for(int j=mat->p[i]; j < mat->p[i+1] ; j++)
	printf(" %d",mat->i[j]+1);
      printf("\n");      
    }
  }
  else{ 
    printf("# binary CSR matrix (%d x %d) in Pair Form: "
	   "nzmax= %d, nz=%d\n",
	   mat->rows,mat->cols,mat->nzmax,mat->nz);
    for(int i=0;i<mat->nz;i++)
      printf("%d %d\n",mat->p[i]+1,mat->i[i]+1);
  }
}

/** 
 * @brief write a CSR matrix in MMX integer form 
 * @param fout output file name (may be `stdout` in which case ignore `fext`)
 * @param fext file extension [concatenate `name=${fout}${fext}`]
 * @param mat the matrix `must be in Compressed Row form.`

 * @param comment if not `NULL`, insert the comment line(s) with an 
 *   extra `#` at the start of each row 
 */ 
void csr_mm_write( char * const fout, const char fext[], const csr_t * const mat,
		  const char comment[]){
  if(mat->nz>0) /** ensure this is CSR matrix */
    ERROR("list-of-pairs matrix, compress it first");
  int result=0; /**< non-zero if write error */
  size_t len=strlen(fout)+strlen(fext)+1;
  char *str;
  
  FILE *f;
  if(strcmp(fout,"stdout")!=0){
    str=calloc(len, sizeof(char));
    sprintf(str,"%s%s",fout,fext);
    f=fopen(str,"w");
  }
  else{
    str=fout;
    f=stdout;
  }
  
  if(!f)
    ERROR("can't open file '%s' for writing",str);
  if(fprintf(f,"%%%%MatrixMarket matrix coordinate integer general\n")<0)
    result++;
  if(comment!=NULL){
    if(fprintf(f,"%% %s\n",comment)<0)
      result++;
  }
  if(fprintf(f,"%d %d %d\n",mat->rows,mat->cols, mat->p[mat->rows])<3)
    result++;

  for(int i=0;i<mat->rows;i++)
    for(int j=mat->p[i]; j < mat->p[i+1] ; j++)
      if(fprintf(f,"%d %d 1\n",i+1, mat->i[j]+1)<0)
	result++;
  if(result)
    ERROR("error writing to file '%s'",str);
  
  if(strcmp(fout,"stdout")!=0){
    fclose(f);
    free(str);
  }
}

/** @brief read sparse binary matrix in `alist` format.
 * This is the format created by David MacKay, see 
 * (https://www.inference.org.uk/mackay/CodesFiles.html)
 * 
 * @Returns the corresponding CSR matrix `mH` (first in the file)
 */
csr_t *csr_alist_read(const char fnam[], csr_t * mat, int transpose, int debug){
  int N, M; /** rows, columns */
  int *num_mlist; /** num of non-zero entries per row */
  int biggest_n, biggest_m; /** max nz entries per col, row */ 
  int num_pairs=0;
  long int linenum=0;
  FILE *f = fopen(fnam, "r");
  if(f==NULL)
    ERROR("can't open the (alist) file %s for reading\n",fnam);
  
  linenum++;
  if(2!=fscanf(f," %d %d \n",&N,&M))
    ERROR("\n%s:%ld: invalid alist file",fnam, linenum);
  
  linenum++;
  if(2!=fscanf(f," %d %d \n",&biggest_m,&biggest_n))
    ERROR("\n%s:%ld: invalid alist file",fnam, linenum);

  linenum++;
  num_mlist = malloc(N*sizeof(int));
  if(!num_mlist) ERROR("memory allocation fail");
  for(int i=0; i<N; i++){
    int val=0;
    if((1!=fscanf(f," %d ",&val))||(val>biggest_m))
      ERROR("\n%s:%ld: invalid alist file",fnam, linenum);
    else{
      num_pairs += val;
      num_mlist[i] = val;
    }
  }

  linenum++;
  int nz=0;
  for(int i=0; i<M; i++){
    int val=0;
    if((1!=fscanf(f," %d ",&val))||(val>biggest_n))
      {
	printf("val=%d i=%d N=%d M=%d biggest_n=%d\n",val,i,N,M, biggest_n);	
	ERROR("\n%s:%ld: invalid alist file",fnam, linenum);
      }
    
    else
      nz+=val;
    /** otherwise `ignore` column weight entries */
  }
  assert(nz==num_pairs);
  
  int_pair * inH = malloc(nz*sizeof(int_pair));
  if (!inH)
    ERROR("memory allocation fail");
  
  int idx=0;
  for(int ir=0; ir<N; ir++){
    linenum++;
    for(int j=0; j<biggest_m; j++){
      int val=0;
      if((1!=fscanf(f," %d ",&val))||(val>M))
	ERROR("\n%s:%ld: invalid alist file",fnam, linenum);
      if(j<num_mlist[ir])
	{
	  if(!transpose)
	    inH[idx++] = (int_pair) {val, ir};
	  else 
	    inH[idx++] = (int_pair) {ir, val};
	}
      else{
	if(val!=0)
	  ERROR("\n%s:%ld: invalid alist file",fnam, linenum);
      }
    }      
  }
  if(!transpose)
    mat = csr_from_pairs(mat, idx, inH, M, N);
  else
    mat = csr_from_pairs(mat, idx, inH, N, M);
    
  free(inH);
  free(num_mlist);

  if(debug &1)
    printf("# read alist file %s %s: rows=%d cols=%d nz=%d\n",
	   fnam,!transpose?"(transposed)":"",mat->rows,mat->cols,nz);

  return mat;    
}


/**
 * read sparse matrix into a (binary) CSR (all entries default to 1)
 * (re)allocate mat if needed
 * will also try to read matrix in `alist` format 
 * use transpose=1 to transpose.
 */
csr_t *csr_mm_read(char *fnam, csr_t *mat, int transpose, int debug){
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;   
 
  if ((f = fopen(fnam, "r")) == NULL) 
    ERROR("can't open file %s",fnam);

  if (mm_read_banner(f, &matcode) != 0){
    /** try to read in `alist` format */
    if(debug&1)
      printf("# Could not process Matrix Market banner; try 'alist' format\n");
    fclose(f);
    mat = csr_alist_read(fnam,mat,transpose,debug);    
    if(!mat)
      ERROR("Could not process Matrix Market banner or load 'alist' file.");
    else
      return mat;    
  }
  
  if (!(mm_is_matrix(matcode) && mm_is_sparse(matcode) && 
	mm_is_integer(matcode) && mm_is_general(matcode) )){
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    ERROR("input file %s",fnam);
    exit(1);
  }

  /* find out size of sparse matrix .... */
  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
    ERROR("Cannot read size in input file %s",fnam);

  if(transpose) {int tmp=M;M=N;N=tmp;} /* swap M and N */
  mat=csr_init(mat,M,N,nz); /* at this point mat will fit the data */  

  if(transpose)
    for (int i=0;i<nz;i++){
      int ir,ic,iv; 
      ret_code = fscanf(f, "%d %d %d\n", &ir , &ic, &iv); 
      if((ret_code != 3)||(iv!=1))
	ERROR("at i=%d: %d %d -> %d",i,ir,ic,iv);
      mat->p[i]=ic-1; 
      mat->i[i]=ir-1; 
    }
  else 
    for (int i=0;i<nz;i++){
      int ir,ic,iv; 
      ret_code = fscanf(f, "%d %d %d\n", &ir , &ic, &iv); 
      if((ret_code != 3)||(iv!=1))
	ERROR("at i=%d: %d %d -> %d",i,ir,ic,iv);
      mat->p[i]=ir-1; 
      mat->i[i]=ic-1; 
    }    
  mat->nz=nz;  // csr_out(mat);
  csr_compress(mat); /* sort entries by row */
  // csr_out(mat);
  fclose(f);
  return mat;
}

/** verify that line has only space */
int all_space(const char * str) {
  while (*str) 
    if (!isspace(*str++)) 
      return 0;    
  return 1;
}

/** @brief read up to `lmax` lines from a file in `01` format

 * read up to `lmax` binary vectors of length `m` from a `01` file
 * `fin` open for reading.  Place the vectors as rows (or columns) of
 * matrix `M` of size `lmax` rows by `m` columns (`m` rows by `lmax`
 * columns if `by_col` is non-zero).  Lines starting with `#` are
 * silently ignored; a non-`01` line, or a `01` line of an incorrect
 * length will give an error.
 *
 * @param M initialized output matrix with `lmax` rows and `m` columns
 * @param fin file with 01 data open for reading
 * @param[input,output] lineno current line number in the file.
 * @param fnam file name (for debugging purposes)
 * @param by_col read the file into columns of the matrix 
 * @param debug if `(debug&8 !=0)` print some additonal information.
 * @return the number of rows actually read.
 *
 */
int read_01(mzd_t *M, FILE *fin, long long int *lineno, const char* fnam,
	    const int by_col, const int debug){
  if(!M)
    ERROR("expected initialized matrix 'M'!\n");
  else
    mzd_set_ui(M,0);
  int m   = by_col ? M->nrows : M->ncols;
  int lmax= by_col ? M->ncols : M->nrows;
  int il=0;
  if(!fin)
    ERROR("file 'fin' named '%s' must be open for reading\n",fnam);
  if(debug&8) /** file io */
    printf("# about to read 01 data from line %lld in file '%s'\n",
           *lineno,fnam);

  char *buf=NULL;
  size_t bufsiz=0;

  ssize_t linelen;
  while((il<lmax) && (!feof(fin)) &&
        ((linelen = getline(&buf, &bufsiz, fin))>=0)){ /* including terminating '\n' */
    (*lineno)++;
    switch(buf[0]){
    case '0': case '1':
      if(linelen != m+1)
	ERROR("incorrect line length=%zu, expected %d 01 characters\n"
	      "%s:%lld:1: '%s'\n", linelen, m,fnam,*lineno,buf);
      else{
	int len = linelen-1 < m ? linelen-1 : m;
	if(by_col){
	  for(int i=0; i<len; i++){
	    if (buf[i]=='1')	    
	      mzd_write_bit(M,i,il,1); /** row `i`, col `il` */
	  
	    else if (buf[i]!='0')
	      ERROR("invalid 01 line\n"
		    "%s:%lld:%d: '%s'\n", fnam,*lineno,i+1,buf);
	  }
	}
	else{ /** by rows */
	  for(int i=0; i<len; i++){
	    if (buf[i]=='1')	    
	      mzd_write_bit(M,il,i,1); /** row `il`, col `i` */
	    
	    else if (buf[i]!='0')
	      ERROR("invalid 01 line\n"
		    "%s:%lld:%d: '%s'\n", fnam,*lineno,i+1,buf);
	  }
	}
	(il)++; /** success */
      }
      break;
    case '#':       /** do nothing - skip this line */
      break;
    default:
      if (!all_space(buf))
	ERROR("invalid 01 line\n"
	      "%s:%lld:1: '%s'\n", fnam,*lineno,buf);
      break;
    }
  }
  
  if(debug&8) /** file io */
    printf("# read %d 01 rows from file '%s'\n",il,fnam);
  if(buf)
    free(buf);
  return il;
}


/** @brief write rows or cols of `M` to file in `01` format
 *
 * file `fout` named `fnam` should already be open for writing.
 *
 * @param M initialized input matrix with at least `lmax` = `lmin` + `count` cols 
 * @param by_cols if non-zero, write one column per line, otherwise one row per line 
 * @param fout file  open for writing
 * @param fnam file name (for debugging purposes)
 *
 */
void mzd_write_01(FILE *fout, const mzd_t * const M, const int by_cols, const char* fnam, const int debug){
  if(!M)
    ERROR("expected initialized matrix 'M'!\n");
  if(!fout)
    ERROR("file 'fout' named %s must be open for writing\n",fnam);
  const int n1 = by_cols ? M->ncols : M->nrows ;
  const int n2 = by_cols ? M->nrows : M->ncols ;  
  for(int i = 0; i < n1; i++){
    for (int j=0; j < n2; j++){
      int bit = by_cols ? mzd_read_bit(M,j,i) : mzd_read_bit(M,i,j) ;
      unsigned char ch = bit ? '1' : '0';
	if(ch != fputc(ch,fout)){
	  int err=ferror(fout);
	  printf("file write error %d:\n%s\n",err,strerror(err));
	  ERROR("error writing to file %s\n",fnam); 
	};
      }
      if('\n' != fputc('\n',fout)){
	int err=ferror(fout);
	printf("file write error %d:\n%s\n",err,strerror(err));
	ERROR("error writing to file %s\n",fnam); 
      };
    }
  if (debug&8)
    printf("# wrote %d rows to 01 file %s\n", n1, fnam);
}

/** @brief write a line of `count` zeros to open file `fout` named `fnam` */
void write_01_zeros(FILE *fout, const int count, const char * fnam){
  for (int j=0; j < count; j++){
    if('0' != fputc('0',fout)){
      int err=ferror(fout);
      printf("file write error %d:\n%s\n",err,strerror(err));
      ERROR("error writing to file %s\n",fnam); 
    };
  }
  if('\n' != fputc('\n',fout)){
    int err=ferror(fout);
    printf("file write error %d:\n%s\n",err,strerror(err));
    ERROR("error writing to file %s\n",fnam); 
  };
}

/** 
 * Permute columns of a CSR matrix with permutation perm.
 */
csr_t *csr_apply_perm(csr_t *dst, const csr_t * const src, const mzp_t * const perm){
  int m=src->rows, n=src->cols;
  if (src->nz!=-1)
    ERROR("pair format unsupported, nz=%d; expected \"-1\"",src->nz);
  int nz=src->p[m];
  if(src->cols != perm-> length)
    ERROR("perm length=%d should match cols=%d",perm->length,src->cols);
  if (dst == NULL) 
    dst = csr_init(dst,m,n,nz);
  else if ((dst->rows != m) || (dst->cols != n) || (nz > dst->nzmax)) 
    ERROR("Wrong size for return matrix.\n");
  for(int i=0;i<m;i++)    
    for(int j=src->p[i]; j < src->p[i+1] ; j++){
      dst->p[j]=i; // pair format, to be compressed later
      dst->i[j]= perm->values[src->i[j]];
      //      printf("row %d: mapping %d to %d\n",i,src->i[j],perm->values[src->i[j]]);
    }
  dst->nz=nz;
  //  csr_out(dst);
  csr_compress(dst);
  //csr_out(dst);
  return dst;
}



/** 
 * Check if row is linearly dependent with the rows of `matP0`
 * which is assumed to be in standard form.
 * `rankP0` is the number of non-zero rows in `matP0`.
 *
 * @return -1 if empty line after simplification otherwise position of
 * 	the first non-zero bit in the simplified vector
 * TODO: verify that it works as claimed 
 */

int do_reduce(mzd_t *row, const mzd_t *matP0, const rci_t rankP0){
  word * rawrow = mzd_row(row,0);  
  rci_t j=0;
  rci_t n=row->ncols;
  do{
    j=nextelement(rawrow,row->width,j);
    if(j==-1) // empty line after simplification
      return j; 
    else if (j<rankP0)
      mzd_combine_even_in_place(row,0,0,matP0,j,0);
    j++;
  } while (j < rankP0);
  if (j<n)
    return j;
  return -1;
}

/**
 * generate binary error vector with error probability p 
 * TODO: use a different generator; also double in `[0,1)`.
 */
void make_err(mzd_t *row, double p){
  const int max=1000000;
  int cutoff=p*max;
  if (row==NULL)
    ERROR("make_err: row vector must be allocated to correct size!");
  rci_t n=row->ncols, m=row->nrows;
  for(rci_t i=0;i<m;i++){
    mzd_row_clear_offset(row,i,0);
    for(rci_t j=0;j<n;j++){
      if (rand_uniform(max)<cutoff)
	mzd_write_bit(row,i,j,1);
    }
  }
}

/** @create column error vector(s) */
void make_err_vec(mzd_t *errors, double *vP){
  assert(errors!=NULL);
  const rci_t n=errors->nrows;
  const rci_t m=errors->ncols;
  for(rci_t i=0; i<n; i++){
    mzd_row_clear_offset(errors,i,0);
    const double prob=vP[i];
    for(rci_t j=0; j<m; j++){
      if (drandom()<prob)
	mzd_write_bit(errors,i,j,1);
    }
  }  
}

/** @brief create a sparse error vector (indices of set bits) 
 * @param siz pointer to the initial allocated size of `vec`
 * @param num on return, pointer to the number of entries in vec
 * @param vec vector of indices to create (content destroyed)
 * @param nvec the length of the error vector 
 * @param prob the probability of set bit 
 * @return the pointer to the created vector 
 */
int *do_sparse_rnd_vec(int *siz, int *num, int *vec, const int nvec, const double prob){
    int ivec=0;
    const double onebyL = -1.0/log(1.0 - prob);
    assert(onebyL>0);
    int j =  (int )floor(onebyL * rnd_exponential());
    while(j < nvec){
      if(ivec >= *siz){
	*siz *= 2;
	vec=realloc(vec,(*siz) * sizeof(int));
      }
      vec[ivec++]=j;
      j += (int )ceil(onebyL * rnd_exponential());
    }
    *num=ivec;
    return vec;
}

/** @brief create a sample of errors to play with.
 *  @param mHe matrix with `nvec` columns to return the syndrome `H*e`
 *  @param mLe matrix with `nvec` columns for logical error `L*e`
 *  @param Ht, Lt the initial matrices (transposed)
 * @return 0
 */
int do_errors(mzd_t *mHe, mzd_t *mLe, const csr_t * const Ht, const csr_t * const Lt,
	      const double vP[]){
  
  assert((mHe!=NULL) && (mLe!=NULL)); /** sanity check */
  assert(mHe->ncols == mLe->ncols);   /** how many errors to produce */
  assert((Lt!=NULL) && (Lt->cols == mLe->nrows));     /** rows `L` */
  assert((Ht!=NULL) && (Ht->cols == mHe->nrows));     /** rows `H` */
  assert(vP!=NULL);
  
  int max = 100;  /** initial size of `vec` */
  int * vec = malloc(max * sizeof(int));

  mzd_set_ui(mHe,0); /** zero matrix */
  mzd_set_ui(mLe,0); /** zero matrix */
  int nvec = mHe->ncols;
  /** for each error type (column of `H` = row of `Ht`) */
  for(int i=0; i < Ht->rows; i++){
    int ivec=0;
    vec = do_sparse_rnd_vec(&max, &ivec, vec, nvec, vP[i]);

    /** flip the bits in `mHe` row by row to speed it up */
    for(int ir = Ht->p[i]; ir < Ht->p[i+1]; ir++){
      const int irow = Ht->i[ir]; /** column index in `H` = row in `Ht` */
      for(int j=0; j < ivec; j++)
	mzd_flip_bit(mHe,irow,vec[j]);
    }

    /** flip the bits in `mLe` row by row */
    for(int ir = Lt->p[i]; ir < Lt->p[i+1]; ir++){
      int irow = Lt->i[ir]; /** column index in `mL` */
      for(int j=0; j<ivec; j++)
	mzd_flip_bit(mLe,irow,vec[j]);
    }
  }
  free(vec);
  return 0;
}

/** @brief calculate the rank of the csr matrix `M` */
int rank_csr(const csr_t * const M){
  assert(M);
  mzd_t *mzd_M = mzd_from_csr(NULL, M);
  int rank=mzd_gauss_delayed(mzd_M,0,0);
  mzd_free(mzd_M);
  return rank;
}
