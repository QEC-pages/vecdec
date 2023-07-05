/* 	$Id: util_m4ri.c,v 1.1 2017/10/31 21:56:50 leonid Exp leonid $	 */
/************************************************************************ 
 * helper functions for use with m4ri library, 
 * including binary sparse matrices and conversion utilities 
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu> 
 * some code borrowed from various sources 
 ************************************************************************/
#include <m4ri/m4ri.h>
#include <stdio.h>
// #include <"copy_m4ri.h"
#include "mmio.h"

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

size_t mzd_weight(const mzd_t *A){
  size_t count = 0;
  if(A->width == 1) {
    for(rci_t i = 0; i < A->nrows; ++i)
      for(rci_t j = 0; j < A->ncols; ++j)
	if(mzd_read_bit(A, i, j))
	  ++count;
    return (count);
  }
  for(rci_t i = 0; i < A->nrows; ++i) {
    word *truerow = A->rows[i];
    for(wi_t j = 0; j < A->width - 1; j ++)
      count += m4ri_bitcount(truerow[j]);

    for(int j = 0; j < A->ncols % m4ri_radix; ++j)
      if(mzd_read_bit(A, i, m4ri_radix * (A->ncols / m4ri_radix) + j))
	++count;
  }

  return count ;
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


/** 
 * transpose compressed CSR matrix, 
 * (re) allocate if needed 
 * return resulting matrix
 * TODO: add code for List of Pairs 
 */
csr_t * csr_transpose(csr_t *dst, const csr_t *p){
  int m=p->rows, n=p->cols, nz=p->p[m];
  if (dst == NULL) 
    dst = csr_init(NULL,n,m,nz);
  else if ((dst->cols != m) || (dst->rows != n) || (dst->nzmax < MAX(nz, n+1))) 
    ERROR("Wrong size for return matrix.\n");
  else
    dst->nz=0; /* clear matrix */
  for(int i=0;i<m;i++)    
    for(int j=p->p[i]; j < p->p[i+1] ; j++){
      dst->p[j]=p->i[j]; // pair format, to be compressed later
      dst->i[j]=i ;
    }
  dst->nz=nz;
  csr_compress(dst);
  return dst;
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
mzd_t *mzd_generator_from_csr(mzd_t *G, csr_t *H){  
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
  }
  if(clear)
    mzd_set_ui(C,0); 
  rci_t const m = S->rows;
  
  for(rci_t i = 0; i < m; ++i)
    for(int j=S->p[i]; j < S->p[i+1] ; j++)
      mzd_combine(C,i,0, C,i,0, B,S->i[j],0); /* combine the two rows */
  return C;
}


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
  if ((mat!=NULL)&&((mat->nzmax < nzmax)||(mat->nzmax < rows+1)))
    mat=csr_free(mat);  /* allocated size was too small */  
  if(mat==NULL){
    mat=malloc(sizeof(csr_t));  
    mat->p=calloc(MAX(nzmax,(rows+1)),sizeof(int));
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
csr_t * csr_from_pairs(csr_t *mat, int nz, int_pair *prs, int nrows, int ncols){
  mat = csr_init(mat, nrows, ncols, nz);
  qsort(prs, nz, sizeof(int_pair), cmp_int_pairs);
  int i, j=0;
  for(i=0;i < mat->rows; i++){
    mat->p[i]=j;
    while ((prs[j].a == i) && (j<nz)){
      mat->i[j]=prs[j].b;
      j++;
    }
  }
  mat->p[i]=j; /* final value */
  mat->nz=-1; /* indicate compressed form */
  return mat;
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
    while ((pairs[j].a == i)&&(j<nz)){
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
 * read sparse matrix into a (binary) CSR (all entries default to 1)
 * (re)allocate mat if needed
 * use transpose=1 to transpose.
 */
csr_t *csr_mm_read(char *fin, csr_t *mat, int transpose){
  int ret_code;
  MM_typecode matcode;
  FILE *f;
  int M, N, nz;   
 
  if ((f = fopen(fin, "r")) == NULL) 
    ERROR("can't open file %s",fin);

  if (mm_read_banner(f, &matcode) != 0)
    ERROR("Could not process Matrix Market banner.");

  if (!(mm_is_matrix(matcode) && mm_is_sparse(matcode) && 
	mm_is_integer(matcode) && mm_is_general(matcode) )){
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    ERROR("input file %s",fin);
    exit(1);
  }

  /* find out size of sparse matrix .... */
  if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
    ERROR("Cannot read size in input file %s",fin);

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
  return mat;
}

/** 
 * Permute columns of a CSR matrix with permutation perm.
 */
csr_t *csr_apply_perm(csr_t *dst, csr_t *src, mzp_t *perm){
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
 * Check if row is linearly dependent with the rows of matP0
 * which is assumed to be in standard form.
 * rankP0 is the number of non-zero rows in rankP0.
 * return: -1 if empty line after simplification
 * 	otherwise position of the first non-zero bit in the simplified vector 
 */

int do_reduce(mzd_t *row, const mzd_t *matP0, const rci_t rankP0){
  word * rawrow = row->rows[0];  
  rci_t j=-1;
  rci_t n=row->ncols;
  do{
    j=nextelement(rawrow,row->width,j);
    if(j==-1) // empty line after simplification
      return j; 
    else if (j<rankP0)
      mzd_combine_even_in_place(row,0,0,matP0,j,0);
  } while (j < rankP0);
  if (j<n)
    return j;
  return -1;
}

/*
 * generate binary error vector with error probability p 
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

