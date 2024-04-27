/**
 *  @file vecdec.c
 *
 * @brief vecdec - a simple vectorized decoder
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2023 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 *
 */

#include <inttypes.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <m4ri/m4ri.h>
#include <m4ri/mzd.h>
#include "utils.h"
#include "util_m4ri.h"
#include "vecdec.h"
#include "qllr.h"

params_t prm={ .nchk=0, .nvar=0, .ncws=0, .steps=50,
  .lerr=-1, .maxosd=100, .bpalpha=0.5, .bpretry=1, .swait=0, .maxC=0,
  .dW=0, .minW=INT_MAX, .dE=-1, .dEdbl=-1, .minE=INT_MAX,
  .nvec=1024, .ntot=1, .nfail=0, .seed=0, .useP=0, .dmin=0,
  .debug=1, .fdem=NULL, .fdet=NULL, .fobs=NULL, .fout="tmp", .ferr=NULL,
  .mode=0, .submode=0, .use_stdout=0, 
  .LLRmin=0, .LLRmax=0, .codewords=NULL, .num_cws=0,
  .finH=NULL, .finL=NULL, .finG=NULL, .finK=NULL, .finP=NULL,
  .finC=NULL, .outC=NULL, 
  .vP=NULL, .vLLR=NULL, .mH=NULL, .mHt=NULL,
  .mL=NULL, .mLt=NULL, .internal=0, 
  .file_det=NULL, .file_obs=NULL, .line_det=0, .line_obs=0,
  .mE=NULL, .mHe=NULL, .mLe=NULL, .mHeT=NULL, .mLeT=NULL,
  .nzH=0, .nzL=0
};

/** various success counters */
long long int cnt[EXTR_MAX];
long long int iter1[EXTR_MAX]; /** sums of BP iteration numbers */
long long int iter2[EXTR_MAX]; /** sums of BP iteration numbers squared */



/** @brief calculate the energy of the row `i` in `A` */
static inline qllr_t mzd_row_energ_naive(qllr_t *coeff, const mzd_t *A, const int i){
  qllr_t ans=0;
  //  mzd_print(A);
  for(rci_t j = 0; j < A->ncols; ++j)
    if(mzd_read_bit(A, i, j)){
      //      printf("non-zero i=%d j=%d coeff=%d\n",i,j,coeff[j]);
      ans += coeff[j];
    }
  return (ans);
}

#if 0
static inline qllr_t mzd_row_energ(qllr_t *coeff, const mzd_t *A, const int i){
  return mzd_row_energ_naive(coeff, A, i);
}
#else /** substantial speed-up of `mode=0` calculation for large matrices */
/** @brief calculate the energy of the row `i` in `A` */
qllr_t mzd_row_energ(qllr_t *coeff, const mzd_t *A, const int i){
#ifndef NDEBUG  
  if (mzd_is_windowed(A))
    ERROR("this does not work on a windowed matrix, use `mzd_row_energ_naive()`");  
#endif   
  qllr_t ans=0;
  word const * truerow = mzd_row(A, i);
  //  mzd_print(A);
  //  for(rci_t j = 0; j < A->ncols; ++j)
  int j=-1;
  while((j=nextelement(truerow,A->width,j+1)) != -1){
    if(j >= A->ncols)
      break;
    ans += coeff[j];
  }
#ifndef NDEBUG
#  ifdef USE_QLLR
  assert(ans==mzd_row_energ_naive(coeff,A,i));
#  else
  assert(fabs(ans-mzd_row_energ_naive(coeff,A,i)) < 0.001 * prm.LLRmin);
#  endif 
#endif   
  return (ans);
}
#endif /* if 1 */


/** @brief calculate the probability of codeword in `one_vec_t` */
static inline double do_prob_one_vec(const one_vec_t * const pvec, const params_t * const p){
  double ans=1;
  const int * idx = pvec->arr;
  for(int i=0; i< pvec->weight; i++, idx++){
    const double vP = p->vP[*idx];
    ans *= 2*sqrt(vP*(1 - vP));
  }
  return ans;
}

/**
 * @brief one step of gauss on column `idx` of matrix `M`
 * @param M the matrix
 * @param idx index of the column of `M` to deal with
 * @param begrow row to start with
 * @return number of pivot points found, `0` or `1` only
 */
static inline int gauss_one(mzd_t *M, const int idx, const int begrow){
  /** note: force-inlining actually slows it down (`???`) */
  rci_t startrow = begrow;
  rci_t pivots = 0;
  const rci_t i = idx;
  //  for (rci_t i = startcol; i < endcol ; ++i) {
  for(rci_t j = startrow ; j < M->nrows; ++j) {
    if (mzd_read_bit(M, j, i)) {
      mzd_row_swap(M, startrow, j);
      ++pivots;
      for(rci_t ii = 0 ;  ii < M->nrows; ++ii) {
        if (ii != startrow) {
          if (mzd_read_bit(M, ii, i)) {
            mzd_row_add_offset(M, ii, startrow,0);
          }
        }
      }
      startrow = startrow + 1;
      break;
    }
  }
  //  }
  return pivots; /** 0 or 1 only */
  // if one, need to update the current pivot list
}

/**
 * @brief one step of gauss on column `idx` of two-block matrix `[M|S]`
 * @param M the first block (check matrix)
 * @param S the second block (syndromes)
 * @param idx index of the column of `M` to deal with
 * @param begrow row to start with
 * @return number of pivot points found, `0` or `1` only
 */
static inline int twomat_gauss_one(mzd_t *M, mzd_t *S, const int idx, const int begrow){
  /** note: force-inlining actually slows it down (`???`) */
  rci_t startrow = begrow;
  rci_t pivots = 0;
  const rci_t i = idx;
  //  for (rci_t i = startcol; i < endcol ; ++i) {
  for(rci_t j = startrow ; j < M->nrows; ++j) {
    if (mzd_read_bit(M, j, i)) {
      mzd_row_swap(M, startrow, j);
      mzd_row_swap(S, startrow, j);
      ++pivots;
      for(rci_t ii = 0 ;  ii < M->nrows; ++ii) {
        if (ii != startrow) {
          if (mzd_read_bit(M, ii, i)) {
            mzd_row_add_offset(M, ii, startrow,0);
            mzd_row_add_offset(S, ii, startrow,0);
          }
        }
      }
      startrow = startrow + 1;
      break;
    }
  }
  //  }
  return pivots; /** 0 or 1 only */
  // if one, need to update the current pivot list
}

 
/** @brief return permutation = decreasing probabilities (increasing LLR) */
mzp_t * sort_by_llr(mzp_t *perm, const qllr_t vLLR[], params_t const * const p){
  assert(perm->length == p->nvar);
  /** prepare array of ippairs */
  ippair_t * pairs = malloc(p->nvar * sizeof(ippair_t));
  if (!pairs)
    ERROR("memory allocation failed\n");
  for(int i=0; i<p->nvar; i++){
    pairs[i].index = i;
    //    pairs[i].prob = p->vP[i];
    pairs[i].llr = vLLR[i];
  }
  qsort(pairs, p->nvar, sizeof(ippair_t), cmp_ippairs);
  for(int i=0; i<p->nvar; i++)
    perm->values[i] = pairs[i].index;

  //  for(int i=0; i<p->nvar; i++) printf("i=%d p[perm[i]]=%g\n",i,p->vP[perm->values[i]]);
  free(pairs);
  return perm;
}

/** @brief prepare an ordered pivot-skip list of length `n-rank` */
mzp_t * do_skip_pivs(const size_t rank, const mzp_t * const pivs){
  const rci_t n=pivs->length;
  rci_t j1=rank; /** position to insert the next result */
  mzp_t * ans = mzp_copy(NULL,pivs);
  qsort(ans->values, rank, sizeof(pivs->values[0]), cmp_rci_t);

  for(rci_t j=0; j<n; j++){
    if(!bsearch(&j, ans->values, rank, sizeof(ans->values[0]), cmp_rci_t)){
      ans->values[j1++]=j;
    }
  }
  assert(j1==n);

  int j=rank;
  for(size_t i=0; j<n; j++)
    ans->values[i++] = ans->values[j];
  ans->length = n-rank;

  if(prm.debug & 8){/** in skip_pivs */
    printf("skip_pivs of len=%d: ",ans->length);
    for(int i=0; i< ans->length; i++)
      printf(" %d%s",ans->values[i],i+1 == ans->length ?"\n":"");
    printf("pivs of len=%d, rank=%zu: ",pivs->length, rank);
    for(size_t i=0; i< rank; i++)
      printf(" %d%s",pivs->values[i],i+1 == rank ?"\n":"");
  }
  return ans;
}

/** @brief Read in vectors from `NZLIST` file `fnam` to hash.
    @return number of entries read */
long long int nzlist_read(const char fnam[], params_t *p){
  long long int count=0, lineno;
  assert(fnam);
  FILE * f=nzlist_r_open(fnam, &lineno); /** no need to check if open */
  one_vec_t *entry=NULL;
  while((entry=nzlist_r_one(f,NULL, fnam, &lineno))){
    if((p->maxC) && (count >= p->maxC))
      ERROR("number of entries in file %s exceeds maxC=%lld\n", p->finC, p->maxC);

    qllr_t energ=0;
    for(int i=0; i < entry->weight; i++) 
      energ += p->vLLR[entry -> arr[i]];
    entry->energ=energ;
      /** `assume` unique vectors stored in the file */
    const size_t keylen = entry->weight * sizeof(rci_t);
    HASH_ADD(hh, p->codewords, arr, keylen, entry); /** store in the `hash` */
    count ++;
  }
  fclose(f);
  return count;
}

/** @brief Write vectors from hash to `NZLIST` file `fnam`.
    @return number of entries written */
long long int nzlist_write(const char fnam[], const char comment[], params_t *p){
  long long int count=0;
  assert(fnam);
  FILE * f=nzlist_w_new(fnam, comment); /** no need to check if open */
  one_vec_t *pvec;
  
  for(pvec = p->codewords; pvec != NULL; pvec = (one_vec_t *)(pvec->hh.next)){
    count ++;
    nzlist_w_append(f,pvec);
  }
  fclose(f);
  return count;
}

/** @brief update minE and minW values using vectors from hash */
int do_hash_min(params_t * const p){
  /** TODO: do we want to drop codewords not satisfying the W and E bounds? */
  for(one_vec_t *pvec = p->codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
    if(p->minE > pvec->energ)
      p->minE = pvec->energ;
    if(p->minW > pvec->weight)
      p->minW = pvec->weight;
  }
  return p->minW;
}

/** @brief Verify codewords `cz` in the hash.  Valid `c` satisfies `c*Ht=0` and `c*Lt!=0`.
 * @param mHt matrix `H=Hx` (transposed)
 * @param mLt matrix `L=Lx` (transposed) or `NULL` to skip the second part of the check.
 *  
 */
int do_hash_verify_CW(const csr_t * const mHt, const csr_t * const mLt, const params_t * const p){
  const int k= mLt!=NULL ? mLt->cols : 0;
  const int r=mHt->cols;
  if (k) assert(mHt->rows == mLt->rows);
  mzd_t *vHt = mzd_init(1, r);
  mzd_t *vLt = k>0 ? mzd_init(1, k) : NULL;
  one_vec_t *pvec;
  int not_cw = 0; 
  long long int count=0;

  for(pvec = p->codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
    mzd_set_ui(vHt,0);
    for(int i=0; i < pvec->weight; i++){
      const int pos = pvec->arr[i];
      for(int j=mHt->p[pos]; j < mHt->p[pos+1] ; j++)
	mzd_flip_bit(vHt,0,mHt->i[j]);
    }
    if(k){
      mzd_set_ui(vLt,0);
      for(int i=0; i < pvec->weight; i++){
	const int pos = pvec->arr[i];
	for(int j=mLt->p[pos]; j < mLt->p[pos+1] ; j++)
	  mzd_flip_bit(vLt,0,mLt->i[j]);
      }
      not_cw = mzd_is_zero(vLt); 
    }    
    if((!mzd_is_zero(vHt))||(not_cw)){
      printf("v=");
      print_one_vec(pvec);
      printf("v*Ht=");
      mzd_print(vHt);
      if(k){
	printf("v*Lt=");
	mzd_print(vLt);
      }
      ERROR("invalid hash vector[%lld]\n",count);
    }
    count++;
  }

  mzd_free(vHt);
  if (vLt) 
    mzd_free(vLt);

  return 0;
}

/** @brief using codewords in the hash, estimate fail probability */
double do_hash_fail_prob( params_t * const p){

  one_vec_t *pvec;
  
  if(p->debug & 1024) {/** `print` the list of cws found by energy */
    HASH_SORT(p->codewords, by_energy);
    for(pvec = p->codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next))
      print_one_vec(pvec);
  }
  /** finally calculate and output fail probability here */
  /** TODO: use the limit on `W` and `E` (just ignore codewords outside the limit) */
  double pfail=0, pmax=0;
  int minW=p->nvar + 1;
  for(pvec = p->codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
    double prob=do_prob_one_vec(pvec, p);
    pfail += prob;
    if(prob>pmax)
      pmax=prob;
    if(minW> pvec->weight)
      minW=pvec->weight;
  }
  /** todo: prefactor calculation */
  if(p->debug&1)
    printf("# sumP(fail) maxP(fail) min_weight num_found\n");
  printf("%g %g %d %lld\n", pfail, pmax, minW, p->num_cws);
  
  return pfail; 
}

void do_hash_clear(params_t *const p){
    /** prescribed way to clean the hashing table */
  one_vec_t *cw, *tmp;
  HASH_ITER(hh, p->codewords, cw, tmp) {
    HASH_DEL(p->codewords, cw);
    free(cw);
  }


}

/** @brief see if the codeword needs to be added to hash, return pointer */
one_vec_t * do_hash_check(const int ee[], int weight, params_t * const p){
  const size_t keylen = weight * sizeof(rci_t);
  one_vec_t *pvec=NULL;
  HASH_FIND(hh, p->codewords, ee, keylen, pvec);
  if(pvec){
    pvec->cnt++; /** just increment the counter */
    if(p->debug &16)
      printf("vector exists, cnt=%d\n",pvec->cnt);
  }
  else{ /** vector not found, inserting */
    qllr_t energ=0;
    for(int i=0; i<weight; i++) 
      energ += p->vLLR[ee[i]];
    if((p->dE < 0) || (energ <= p->minE + p->dE)){      
      pvec = (one_vec_t *) malloc(sizeof(one_vec_t)+keylen);
      if(!pvec)
	ERROR("memory allocation failed!\n");
      pvec->energ = energ; /** energy value */
      pvec->weight = weight;
      pvec->cnt = 1; /** encountered `1`st time */
      memcpy(pvec->arr, ee, keylen);
      HASH_ADD(hh, p->codewords, arr, keylen, pvec); /** store in the `hash` */
      ++(p->num_cws);     /** update the counter */
    }
    else{
      /** silently ignore */
    }
  }
  return pvec;
}


/** @brief Random Information Set search for small-E logical operators.
 *
 *  Uses hashing storage to identify unique vectors.  Only vectors of weight no
 *  more that `minW`+`dW` will be recorded, where `minW` is the current minimum
 *  weight.
 * 
 * @param dW weight increment from the minimum found
 * @param p pointer to global parameters structure
 * @return minimum `weight` of a CW found (or `-weigt` if early termination condition is reached). 
 */
int do_LLR_dist(params_t  * const p){
  const int dW=p->dW;
  //  qllr_t dE=p->dE;
  /** whether to verify logical ops as a vector or individually */
  const int use_vector = p->mLt->cols >= 16 ? 1 : 0;
  if(p->nvec == 16) /** default value */
    p->nvec=0;
  mzd_t * mH = mzd_from_csr(NULL, p->mH);
  mzd_t *mLt = NULL, *eemLt = NULL, *mL = NULL;
  if(use_vector){ /** all logical ops at once */
    mLt = mzd_from_csr(NULL, p->mLt);
    eemLt = mzd_init(1, p->mLt->cols);
  }
  else /** just one logical op */
    mL = mzd_from_csr(NULL, p->mL);
    
  rci_t *ee = malloc(p->mH->cols*sizeof(rci_t)); /** actual `vector` */
  
  if((!mH) || (!ee))
    ERROR("memory allocation failed!\n");
  //  if(p->debug & 16)  mzd_print(mH);
  /** 1. Construct random column permutation P */

  mzp_t * perm=mzp_init(p->nvar); /** identity column permutation */
  mzp_t * pivs=mzp_init(p->nvar); /** list of pivot columns */
  if((!pivs) || (!perm))
    ERROR("memory allocation failed!\n");

  int iwait=0, ichanged=0;
  perm = sort_by_llr(perm, p->vLLR, p);   /** order of decreasing `p` */
  for (int ii=0; ii< p->steps; ii++){
    if(ii!=0){
      pivs=mzp_rand(pivs); /** random pivots LAPAC-style */
      mzp_set_ui(perm,1);
      perm=perm_p_trans(perm,pivs,0); /**< corresponding permutation */
      /** todo: make it probability-dependent */
    }
    /** full row echelon form of `H` (gauss) using the order in `perm` */
    int rank=0;
    for(int i=0; i< p->nvar; i++){
      int col=perm->values[i];
      int ret=gauss_one(mH, col, rank);
      if(ret)
        pivs->values[rank++]=col;
    }
    /** construct skip-pivot permutation */
    mzp_t * skip_pivs = do_skip_pivs(rank, pivs);

    /** calculate sparse version of each vector (list of positions)
     *  `p`    `p``p`               # pivot columns marked with `p`   
     *  [1  a1        b1 ] ->  [a1  1  a2 a3 0 ]
     *  [   a2  1     b2 ]     [b1  0  b2 b3 1 ]
     *  [   a3     1  b3 ]
     */
    int k = p->nvar - rank;
    for (int ir=0; ir< k; ir++){ /** each row in the dual matrix */
      int cnt=0; /** how many non-zero elements */
      const int col = ee[cnt++] = skip_pivs->values[ir];
      for(int ix=0; ix<rank; ix++){
        if(mzd_read_bit(mH,ix,col))
          ee[cnt++] = pivs->values[ix];          
      }
      /** sort the column indices */
      qsort(ee, cnt, sizeof(rci_t), cmp_rci_t);
      
      /** verify logical operator */
      int nz;
      if(use_vector){ /** use vector */
        mzd_set_ui(eemLt,0);
        for(int i=0; i<cnt; i++) 
          mzd_combine_even_in_place(eemLt,0,0,mLt,ee[i],0);
        nz = mzd_is_zero(eemLt) ? 0 : 1;
      }
      else{ /** for each logical operator = row of `mL` */
	nz=0;
	for(int j=0; j < mL->nrows; j++){
	  nz=0;
	  for(int i=0; i<cnt; i++) /** bits in one row */
	    nz ^= mzd_read_bit(mL,j,ee[i]);
	  if(nz)
	    break;
	}
      }
      if(nz){ /** we got non-trivial codeword! */
        /** TODO: try local search to `lerr` (if 2 or larger) */
        /** calculate the energy and compare */
        /** at this point we have `cnt` codeword indices in `ee`, and its `energ` */
        if (cnt < p->minW){
          p->minW=cnt;
	  if (p->minW <= p->dmin){ /** early termination condition */
	    p->minW = - p->minW; /** this distance value is of little interest; */
	  }
	}
        if((dW<0) || (cnt <= abs(p->minW) + dW)){ /** try to add to hashing storage */
	  one_vec_t *ptr=do_hash_check(ee,cnt,p); 
	  if(ptr->cnt == 1){ /** new codeword just added to hash */
	    ichanged++; 
	    if (ptr->energ < p->minE){  /** legacy code */
	      if(p->debug&1)
		printf("nz=%d cnt=%d energ=%g\n",nz,cnt,dbl_from_llr(ptr->energ));
	      p->minE=ptr->energ;
	    }
	    if (p->minW<0)
	      goto alldone; 
          }
        }
      }
    } /** end of the dual matrix rows loop */
    if(p->debug & 16)
      printf(" round=%d of %d minE=%g minW=%d num_cws=%lld ichanged=%d iwait=%d\n",
             ii+1, p->steps, dbl_from_llr(p->minE), p->minW, p->num_cws, ichanged, iwait);
    
    mzp_free(skip_pivs);
    
    iwait = ichanged > 0 ? 0 : iwait+1 ;
    ichanged=0;
    if((p->swait > 0)&&(iwait > p->swait)){
      if(p->debug & 16)
        printf("  iwait=%d >swait=%d, terminating after %d steps\n", iwait, p->swait, ii+1);
      break;
    }
    if((p->maxC > 0) && (p->num_cws >= p->maxC)) 
      break; /** limit was set, not an error */

  }/** end of `steps` random window */

  /** TODO: prefactor calculation */

 alldone: /** early termination label */

  /** clean up */
  mzp_free(perm);
  mzp_free(pivs);
  free(ee);
  if(use_vector){
    mzd_free(eemLt);
    mzd_free(mLt);
  }
  else
    mzd_free(mL);
  mzd_free(mH);
  
  return p->minW;
}

int do_energ_verify(const qllr_t * const vE, const mzd_t * const mE, const params_t * const p){
  int nfail=0;
  for (int i=0;i < mE->nrows;i++){
# ifdef USE_QLLR
    if(vE[i] != mzd_row_energ(p->vLLR,mE,i)){
      nfail++;
#ifndef NDEBUG      
      printf("mismatch i=%d vE=%d row_energ=%d\n",i,vE[i],mzd_row_energ(p->vLLR,mE,i));
#endif       
    }
# else       
    if(fabs(vE[i] - mzd_row_energ(p->vLLR,mE,i))>1e-5){
      nfail++;
#ifndef NDEBUG      
      printf("mismatch i=%d vE=%g row_energ=%g\n",i,vE[i],mzd_row_energ(p->vLLR,mE,i));
#endif       
    }
# endif
  }
  return nfail;
}

/**
 * @brief do local search up to `p->lerr` inclusive recursively
 * @param vE0 best energy values for each `e` a row in `mE0`
 * @param mE0 best error vectors so far for each syndrome (rows)
 * @param jstart start local search from this column in `mH`
 * @param lev recusion level, must not exceed `p->lerr`
 * @param vE current energy values for each `e` a row in `mE`
 * @param mE input error vectors (rows)
 * @param mH check matrix in row echelon form (pivots in `pivs`)
 * @param skip_pivs `sorted` list of `n-rank` non-pivot positions in `mH`
 * @param pivs list of `rank` pivots returned by gauss (length = `n`)
 * @param p pointer to the remaining program parameters
 * @return number of updated error vectors
 * @todo: see if binary ops + transposition is faster
 */
int do_local_search(qllr_t *vE0, mzd_t * mE0, rci_t jstart, int lev,
		    const qllr_t * const vE, const mzd_t * const mE, const mzd_t * const mH,
		    const mzp_t * const skip_pivs, const mzp_t * const pivs,
		    const params_t * const p){
  assert(lev<=p->lerr);
  int last_lev = lev < p->lerr ? 0 : 1;
  qllr_t *vE1 = NULL;
  mzd_t *mE1 = NULL; //mzd_copy(NULL,mE); /** error vectors to update */
  if(p->debug&128)
    printf("entering %slev=%d / %d of recursion jstart=%d\n",last_lev? "last " :"",lev,p->lerr, jstart);  
#ifndef NDEBUG
  if(0 != do_energ_verify(vE0,mE0,p)){
    mzd_print(mE0);
    ERROR("energy value mismatch vE0, mE0 at lev=%d of %d \n",lev, p->lerr);
  }
  if(0 != do_energ_verify(vE,mE,p)){
    mzd_print(mE0);
    ERROR("energy value mismatch vE, mE at lev=%d of %d \n",lev, p->lerr);
  }
#endif   
  int ich_here=0, ich_below=0;
  rci_t knum = skip_pivs->length; /** number of non-pivot cols in `mH` to go over */
  rci_t rank = p->nvar - knum; /** number of valid pivot cols */
  rci_t rnum; /** number of non-zero entries in rlis (for each `j`) */
  int * rlis = malloc(rank * sizeof(int));
  if(!rlis) ERROR("memory allocation failed!");
  if(!last_lev){
    vE1 = malloc(sizeof(qllr_t) * mE0->nrows);
    if(!vE1) ERROR("memory allocation failed!");
  }

  for(rci_t j=jstart; j<knum; j++){ /** `j`th `non-pivot` entry */
    if (!last_lev){
      mE1 = mzd_copy(mE1,mE); /** fresh copy of error vectors to update */
      memcpy(vE1,vE,sizeof(qllr_t) * mE0->nrows );
    }
    rci_t jj=skip_pivs->values[j]; /** actual `non-pivot` column we are looking at */
    rnum=0; /** prepare list of positions to update for `jj`th col of `mH` */
    for(rci_t ir=0; ir<rank; ir++)
      if(mzd_read_bit(mH,ir,jj)) /** non-zero bit */
        rlis[rnum++] = pivs->values[ir]; /** column of `mH` to update */
#ifndef NDEBUG
    if (p->debug & 128){
      printf("jj=%d rlis: ",jj);
      for(int ir=0; ir< rnum; ir++)
        printf(" %d%s",rlis[ir],ir+1==rnum?"\n":"");
    }
#endif

    for(rci_t is=0; is < mE->nrows; is++){ /** syndrome rows */
#ifndef NDEBUG      
      if(mzd_read_bit(mE,is,jj)) /** sanity check */
        ERROR("bit found at is=%d jj=%d\n",is,jj);
#endif 
      if(vE0[is] >= 2 * p->LLRmin){/** TODO: add more reasonable cutoff here */
        qllr_t energ = vE[is];
	/** calculate updated energy only */
        energ += p->vLLR[jj];
        for(rci_t ir = 0 ;  ir < rnum; ++ir){
	  const int ii = rlis[ir]; /** position to update */
          if(mzd_read_bit(mE,is,ii))
            energ -= p->vLLR[ii]; /** `1->0` flip */
          else
            energ += p->vLLR[ii]; /** `0->1` flip */
        }
	if(!last_lev){ /* calculate updated error vector */
	  vE1[is]=energ;
	  mzd_flip_bit(mE1,is,jj);
	  for(rci_t ir = 0 ;  ir < rnum; ++ir)
	    mzd_flip_bit(mE1, is, rlis[ir]);
	}
        if(energ < vE0[is]-1e-10){
#ifndef NDEBUG
          if(p->debug & 128){/** inf set decoding */
            printf("lev=%d j=%d jj=%d is=%d E0=%g E=%g success\n", lev,j,jj,is,
		   dbl_from_llr(vE0[is]),dbl_from_llr(energ));
          }
#endif
          vE0[is]=energ;
	  if (!last_lev)
	    mzd_copy_row(mE0,is, mE1,is);
	  else{
	    mzd_copy_row(mE0,is, mE,is);
	    mzd_flip_bit(mE0,is,jj);
	    for(rci_t ir = 0 ;  ir < rnum; ++ir)
	      mzd_flip_bit(mE0, is, rlis[ir]);
	  }
          ich_here++;
        }
      }
    }

    if(!last_lev){ /** go up one recursion level */
      if(j+1<knum){
        ich_below+=do_local_search(vE0,mE0,j+1,lev+1,vE1,mE1, mH, skip_pivs, pivs,p);
      }
    }
  }
  if(p->debug & 128)
    if(ich_here + ich_below)
      printf("exiting lev=%d of recursion, here ch=%d below ch=%d\n",
             lev,ich_here, ich_below);
  free(rlis);
  if (mE1)
    mzd_free(mE1);
  if (vE1)
    free(vE1);
  return ich_below+ich_here;
}

/**
 * @brief actual syndrome-based decoding routine
 * @param mS the matrix with syndromes (each column)
 * @param p structure with error model information
 * @return binary matrix of min weight errors / each syndrome
 ***************** todo: reuse some matrices? ***************/
mzd_t *do_decode(mzd_t *mS, params_t const * const p){
  mzd_t * mH = mzd_from_csr(NULL, p->mH);

  mzd_t * mE = mzd_init(mH->ncols,mS->ncols); /**< error vectors by col */
  qllr_t *vE = calloc(mS->ncols,sizeof(qllr_t)); /**< best energies */
  if((!mE) || (!vE))
    ERROR("memory allocation failed!\n");

  mzp_t * perm=mzp_init(p->nvar); /** identity column permutation */
  mzp_t * pivs=mzp_init(p->nvar); /** list of pivot columns */
  if((!pivs) || (!perm))
    ERROR("memory allocation failed!\n");

  /** first pass ******************************************* */
  perm = sort_by_llr(perm, p->vLLR, p);   /** order of decreasing `p` */
  /** full row echelon form (gauss elimination) using the order of `p`,
   * on the block matrix `[H|S]` (in fact, two matrices).
   */
  int rank=0;
  for(int i=0; i< p->nvar; i++){
    int col=perm->values[i];
    int ret=twomat_gauss_one(mH,mS, col, rank);
    if(ret)
      pivs->values[rank++]=col;
  }
  if((p->debug &4)&&(p->debug &512)){ /** debug gauss */
    printf("rank=%d\n",rank);
    //    printf("perm: "); mzp_out(perm);
    //    printf("pivs: "); mzp_out(pivs);
    //    printf("mH:\n");
    //    mzd_print(mH);
    //    printf("mS:\n");
    //    mzd_print(mS);
  }

  // for each syndrome, calculate error vector and energy
  mzd_set_ui(mE,0); /** zero matrix to store `errors by column` */
  for(int i=0;i< rank; i++)
    mzd_copy_row(mE,pivs->values[i],mS,i);
  mzd_t *mEt0 = mzd_transpose(NULL,mE);
  for(int i=0; i< mS->ncols; i++)
    vE[i]=mzd_row_energ(p->vLLR,mEt0,i);

  if(p->debug & 512){
    printf("mEt0 after round 0:\n");
    mzd_print(mEt0);
  }

  if(p->lerr>0){  /** do information-set decoding `**********************` */
    mzp_t * skip_pivs = do_skip_pivs(rank, pivs);
    mzd_t * mEt1=NULL;
    qllr_t *vE1=NULL;
    if (p->lerr > 1){
      mEt1 = mzd_copy(NULL, mEt0);
      vE1 = malloc(sizeof(qllr_t) * mEt0->nrows);  if(!vE1) ERROR("memory allocation failed!");
      memcpy(vE1,vE,sizeof(qllr_t) * mEt0->nrows);
    }
    do_local_search(vE, mEt0, 0, 1,
		    p->lerr > 1 ? vE1 : vE,
		    p->lerr > 1 ? mEt1 : mEt0, mH, skip_pivs, pivs, p);
    if(p->debug & 512){
      printf("mEt0 after local search:\n");
      mzd_print(mEt0);
    }
    mzp_free(skip_pivs);
    if (p->lerr > 1){
      mzd_free(mEt1);
      free(vE1);
    }
  }

  int iwait=0, ichanged=0;
  /** main loop over permutations * `***************************` */
  for (int ii=1; ii< p->steps; ii++){
    pivs=mzp_rand(pivs); /** random pivots LAPAC-style */
    mzp_set_ui(perm,1); perm=perm_p_trans(perm,pivs,0); /**< corresponding permutation */
    /** todo: make it probability-dependent */
    rank=0;
    for(int i=0; i< p->nvar; i++){
      int col=perm->values[i];
      int ret=twomat_gauss_one(mH,mS, col, rank);
      if(ret)
        pivs->values[rank++]=col;
    }
    // for each syndrome, calculate error vector and energy; update minima
    mzd_set_ui(mE,0); /** zero matrix */
    for(int i=0;i< rank; i++)
      mzd_copy_row(mE,pivs->values[i],mS,i);
    ichanged=0;
    mzd_t * mEt = mzd_transpose(NULL, mE);
    qllr_t *vE1=NULL;
    if(p->lerr > 0){  /** do information-set decoding `**********************` */
      vE1 = malloc(sizeof(qllr_t) * mEt0->nrows);
      if(!vE1) ERROR("memory allocation failed!");
    }
    for(int i=0; i< mEt->nrows; i++){
      qllr_t energ=mzd_row_energ(p->vLLR,mEt,i);
      if(p->lerr > 0)
	vE1[i]=energ;
      if(energ < vE[i]){
	vE[i]=energ;
	mzd_copy_row(mEt0,i,mEt,i);
	ichanged++;
      }
    }
    if(ichanged){
      if(p->debug & 512){
        if(p->debug &4){ /** debug gauss */
          printf("after round %d rank=%d\n",ii,rank);
          //          printf("perm: "); mzp_out(perm);
          //          printf("pivs: "); mzp_out(pivs);
          //          printf("mH:\n");
          //          mzd_print(mH);
          //          printf("mS:\n");
          //          mzd_print(mS);
        }
        printf("mEt0:\n");
        mzd_print(mEt0);
      }
    }
    if(p->lerr > 0){  /** do information-set decoding `**********************` */
      mzp_t * skip_pivs = do_skip_pivs(rank, pivs);
      do_local_search(vE, mEt0, 0, 1, vE1, mEt, mH, skip_pivs, pivs, p);
      if(p->debug & 512){
        printf("mEt0 after local search:\n");
        mzd_print(mEt0);
      }
      mzp_free(skip_pivs);
      if (p->lerr > 1)
	free(vE1);
    }

    mzd_free(mEt);
    iwait = ichanged > 0 ? 0 : iwait+1 ;
    if(p->debug&8) /** convergence information */
      if(ichanged)
        printf(" ii=%d of %d changed=%d\n",ii,p->steps,ichanged);
    if((p->swait > 0)&&(iwait > p->swait))
      break;
  }
  /** clean-up */
  mzd_free(mH);
  mzd_free(mE);
  free(vE);
  mzp_free(perm);
  mzp_free(pivs);
  return mEt0;
}

/** @brief one more matrix initialization routine 
 */
void init_Ht(params_t *p){
  const int n = p->nvar;
  p->mHt = csr_transpose(p->mHt, p->mH);
  p->mLt = csr_transpose(p->mLt,p->mL);
  /** todo: fix reallocation logic to be able to reuse the pointers model */
  //  if(p->vP)    free(p->vP);
  //  p->vP=inP;
  p->vLLR = malloc(n*sizeof(qllr_t));
  assert(p->vLLR !=0);
  p->LLRmin=1e9;
  p->LLRmax=-1e9;

  for(int i=0;  i < n; i++){
    qllr_t val=llr_from_P(p->vP[i]);
    p->vLLR[i] = val;
    if(val < p->LLRmin)
      p->LLRmin = val;
    if(val > p->LLRmax)
      p->LLRmax = val;
  }
  if(p->LLRmin<=0)
    ERROR("LLR values should be positive!  LLRmin=%g LLRmax=%g",
	  dbl_from_llr(p->LLRmin),dbl_from_llr(p->LLRmax));
  if(p->debug & 2){/** `print` out the entire error model ******************** */
    printf("# error model read: r=%d k=%d n=%d LLR min=%g max=%g\n",
	   p->nchk, p->ncws, p->nvar, dbl_from_llr(p->LLRmin),dbl_from_llr(p->LLRmax));  
  }

#ifndef NDEBUG
  if((p->debug & 2)&&(p->debug &512)){ /** print resulting vectors and matrices */
    if(p->useP>0)
      printf("# uniform error probability useP=%g LLR=%g\n",p->useP,dbl_from_llr(p->vLLR[0]));
    else{
      printf("# P=[");
      for(int i=0; i< p->nvar; i++)
	printf("%g%s",p->vP[i],i+1<p->nvar?" ":"]\n");    
      printf("# LLR=[");
      for(int i=0; i< p->nvar; i++)
	printf("%g%s",dbl_from_llr(p->vLLR[i]),i+1<p->nvar?" ":"]\n");
    }
    //    if(init_mat){
    mzd_t *mdH = mzd_from_csr(NULL,p->mH);
    printf("mH:\n");
    //    csr_out(mH);
    mzd_print(mdH);
    mzd_free(mdH);

    printf("mL:\n");
    //    csr_out(mL);
    mzd_t *mdL = mzd_from_csr(NULL,p->mL);
    mzd_print(mdL);
    mzd_free(mdL);
    //    }
  }
#endif

}


/** @brief initialize parameters set at command line */
int var_init(int argc, char **argv, params_t *p){

  int dbg=0;
  long long int lldbg=0;
  double val;

  if(argc<=1)
    ERROR("try \"%s -h\" for help",argv[0]);
#ifdef USE_QLLR
  p->d1=12; p->d2=300; p->d3=7; /** recommended values */
#endif   

  for(int i=1; i<argc; i++){  
    int pos=0;
    if(sscanf(argv[i],"debug=%d",& dbg)==1){/** `debug` */
      if(dbg==0)
	p->debug = 0;
      else{
        if(i==1)
          p->debug = dbg; /** just assign if in the `1st position` */
        else
          p->debug ^= dbg; /** otherwise `XOR` */
        if(p->debug &1)
	  printf("# read %s, debug=%d octal=%o\n",argv[i],p->debug,p->debug);
      }
    }
    else if(sscanf(argv[i],"mode=%d%n",& dbg, &pos)==1){ /** `mode.submode` */
      p->mode = dbg;
      if((int) strlen(argv[i]) > pos){
	 if(sscanf(argv[i]+pos,".%d",& dbg)==1)
	   p->submode = dbg;
	 else
	   ERROR("invalid format argv[%d]: %s\n",i,argv[i]);
      }
      if(p->debug&1)
	printf("# read %s, mode=%d submode=%d [octal=%o]\n",argv[i],
	       p->mode,p->submode,p->submode);      
    }
    else if (sscanf(argv[i],"qllr1=%d",&dbg)==1){ /** QLLR `d1` param */
      p -> d1 = dbg;
      if (p->debug&1)
	printf("# read %s, QLLR parameter d1=%d\n",argv[i],p-> d1);
    }
    else if (sscanf(argv[i],"qllr2=%d",&dbg)==1){ /** QLLR `d2` param */
      p -> d2 = dbg;
      if (p->debug&1)
	printf("# read %s, QLLR parameter d2=%d\n",argv[i],p-> d2);
    }
    else if (sscanf(argv[i],"qllr3=%d",&dbg)==1){ /** QLLD `d3` param */
      p -> d3 = dbg;
      if (p->debug&1)
	printf("# read %s, QLLR parameter d3=%d\n",argv[i],p-> d3);
    }
    else if (sscanf(argv[i],"nvec=%d",&dbg)==1){ /** `nvec` */
      p -> nvec = dbg;
      if (p->debug&1)
	printf("# read %s, nvec=%d\n",argv[i],p-> nvec);
    }
    else if (sscanf(argv[i],"dmin=%d",&dbg)==1){ /** `dmin` */
      p -> dmin = dbg;
      if (dbg<0)
	ERROR("command line argument %d dmin=%d must be non-negative\n",i,dbg);
      if (p->debug&1)
	printf("# read %s, dmin=%d\n",argv[i],p-> dmin);
    }
    else if (sscanf(argv[i],"steps=%d",&dbg)==1){ /** `steps` */
      p -> steps = dbg;
      if (p->debug&1)
	printf("# read %s, steps=%d\n",argv[i],p-> steps);
    }
    else if (sscanf(argv[i],"swait=%d",&dbg)==1){ /** `swait` */
      p -> swait = dbg;
      if (p->debug&1)
	printf("# read %s, swait=%d\n",argv[i],p-> swait);
    }
    else if (sscanf(argv[i],"lerr=%d",&dbg)==1){ /** `lerr` */
      p -> lerr = dbg;
      if (p->debug&1)
	printf("# read %s, lerr=%d\n",argv[i],p-> lerr);
    }
    else if (sscanf(argv[i],"maxosd=%d",&dbg)==1){ /** `maxosd` */
      p -> maxosd = dbg;
      if (p->debug&1)
	printf("# read %s, maxosd=%d\n",argv[i],p-> maxosd);
    }
    else if (sscanf(argv[i],"ntot=%lld",&lldbg)==1){ /** `ntot` */
      p -> ntot = lldbg;
      if (p->debug&1)
	printf("# read %s, ntot=%lld\n",argv[i],p-> ntot);
    }
    else if (sscanf(argv[i],"maxC=%lld",&lldbg)==1){ /** `maxC` */
      p -> maxC = lldbg;
      if (p->debug&1)
	printf("# read %s, maxC=%lld\n",argv[i],p-> maxC);
    }
    else if (sscanf(argv[i],"nfail=%lld",&lldbg)==1){ /** `nfail` */
      p -> nfail = lldbg;
      if (p->debug&1)
	printf("# read %s, nfail=%lld\n",argv[i],p-> nfail);
    }
    else if (sscanf(argv[i],"seed=%lld",&lldbg)==1){ /** `seed` */
      p->seed=lldbg;
      if (p->debug&1)
	printf("# read %s, seed=%lld\n",argv[i],p->seed);
    }
    else if (sscanf(argv[i],"bpalpha=%lg",&val)==1){ /** `bpalpha` */
      p -> bpalpha = val; /** WARNING: no value verification!!! */
      if (p->debug&1)
	printf("# read %s, bpalpha=%g\n",argv[i],p-> bpalpha);
    }
    else if (sscanf(argv[i],"bpretry=%d",&dbg)==1){ /** `bpretry` */
      p -> bpretry = dbg;
      if(dbg<1)
	ERROR("read arg[%d]=%s, bpretry=%d should be at least 1",i,argv[i],dbg);
      if (p->debug&1)
	printf("# read %s, bpretry=%d\n",argv[i],p-> bpretry);
    }
    else if (sscanf(argv[i],"useP=%lg",&val)==1){ /** `useP` */
      p->useP = val;
      if (p->debug&1)
	printf("# read %s, useP=%g\n",argv[i],p-> useP);
    }
    else if (sscanf(argv[i],"dE=%lg",&val)==1){ /** `dE` */
      p -> dEdbl = val;
      if (p->debug&1){
	printf("# read %s, dE=%g\n", argv[i], p->dEdbl);
	if (p->dEdbl < 0)
	  printf("# no limit on error/codeword energies to store\n");
	else
	  printf("# setting upper limit on error/codeword energies to store\n");
      }

    }
    else if (sscanf(argv[i],"dW=%d",&dbg)==1){ /** `dW` */
      p -> dW = dbg;
      if (p->debug&1){
	printf("# read %s, dW=%d\n",argv[i],p-> dW);
	if (p->dW < 0)
	  printf("# no limit on error/codeword weight to store\n");
	else
	  printf("# setting upper limit 'minW+%d' on error/codeword weight to store\n",p->dW);
      }	
    }
    else if (0==strncmp(argv[i],"fout=",5)){ /** `fout` */
      if(strlen(argv[i])>5){
        p->fout = argv[i]+5;
	if (p->debug&1)
	  printf("# read %s, fout=%s\n",argv[i],p->fout);
      }
      else
	ERROR("Please specify argument for 'fout=[string]' w/o space\n");
    }
    else if (0==strncmp(argv[i],"f=",2)){/** back compatibility */
      if(strlen(argv[i])>2)
        p->fdem = argv[i]+2;
      else
        p->fdem = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, (fdem) f=%s\n",argv[i],p->fdem);
    }
    else if (0==strncmp(argv[i],"fdem=",5)){ /** fdem */
      if(strlen(argv[i])>5)
        p->fdem = argv[i]+5;
      else
        p->fdem = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fdem=%s\n",argv[i],p->fdem);
    }
    else if (0==strncmp(argv[i],"finH=",5)){ /** `finH` */
      if(strlen(argv[i])>5)
        p->finH = argv[i]+5;
      else
        p->finH = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finH=%s\n",argv[i],p->finH);
    }
    else if (0==strncmp(argv[i],"finP=",5)){ /** `finP` probabilities */
      if(strlen(argv[i])>5)
        p->finP = argv[i]+5;
      else
        p->finP = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finP=%s\n",argv[i],p->finP);
    }
    else if (0==strncmp(argv[i],"finL=",5)){ /** `finL` logical */
      if(strlen(argv[i])>5)
        p->finL = argv[i]+5;
      else
        p->finL = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finL=%s\n",argv[i],p->finL);
    }
    else if (0==strncmp(argv[i],"finK=",5)){ /** `finK` logical */
      if(strlen(argv[i])>5)
        p->finK = argv[i]+5;
      else
        p->finK = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finK=%s\n",argv[i],p->finK);
    }
    else if (0==strncmp(argv[i],"finC=",5)){ /** `finC` in low-weight codeword list */
      if(strlen(argv[i])>5)
        p->finC = argv[i]+5;
      else
        p->finC = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finC=%s\n",argv[i],p->finC);
    }
    else if (0==strncmp(argv[i],"outC=",5)){ /** `outC` out low-weight codeword list */
      if(strlen(argv[i])>5)
        p->outC = argv[i]+5;
      else
        p->outC = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, outC=%s\n",argv[i],p->outC);
    }
    else if (0==strncmp(argv[i],"finG=",5)){/** `finG` degeneracy generator matrix */
      if(strlen(argv[i])>5)
        p->finG = argv[i]+5;
      else
        p->finG = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finG=%s\n",argv[i],p->finG);
    }
    else if (0==strncmp(argv[i],"fdet=",5)){ /** detector events / 01 file */
      if(strlen(argv[i])>5)
        p->fdet = argv[i]+5;
      else
        p->fdet = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fdet=%s\n",argv[i],p->fdet);
    }
    else if (0==strncmp(argv[i],"fobs=",5)){/** observable events / 01 file */
      if(strlen(argv[i])>5)
        p->fobs = argv[i]+5;
      else
        p->fobs = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fobs=%s\n",argv[i],p->fobs);
    }
    else if (0==strncmp(argv[i],"ferr=",5)){ /** errors / 01 file */
      if(strlen(argv[i])>5)
        p->ferr = argv[i]+5;
      else
        p->ferr = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, ferr=%s\n",argv[i],p->ferr);
    }
    else if((strcmp(argv[i],"--help")==0)
            ||(strcmp(argv[i],"-h")==0)
            ||(strcmp(argv[i],"help")==0)){
      printf( USAGE , argv[0],argv[0]);
      exit (-1);
    }
    else if((strcmp(argv[i],"--morehelp")==0)
            ||(strcmp(argv[i],"-mh")==0)
            ||(strcmp(argv[i],"morehelp")==0)){
      printf( USAGE "\n", argv[0],argv[0]);
      printf( MORE_HELP,argv[0]);
      exit (-1);
    }    
    else{ /* unrecognized option */
      printf("# unrecognized parameter \"%s\" at position %d\n",argv[i],i);
      ERROR("try \"%s -h\" for help",argv[0]);
    }

  }

  
  if (p->seed == 0){
    p->seed=time(NULL)+1000000ul*getpid(); /* ensure a different seed even if started at the same time */
    if((p->debug)&&(p->mode!=3))
      printf("# initializing seed=%lld from time(NULL)+1000000ul*getpid()\n",p->seed);
    /** use `tinymt64_generate_double(&pp.tinymt)` for double [0,1] */
  }
  
  tinymt64_init(&tinymt,p->seed);

  if((! p->fdem) && (! p->finH))
    ERROR("Please specify the DEM file 'fdem' or check matrix file 'finH'\n");
  if((p->fdem) && (p->finH))
    ERROR("Please specify fdem=%s OR finH=%s but not both\n",p->fdem, p->finH);

  if(p->fdem){
    if((p->finL)||(p->finP)||(p->finG)||(p->finH))
      ERROR("DEM file fdem=%s can not be specified together with \n"
	    " finH=%s or finL=%s or finP=%s\n", p->fdem,  
	    p->finH ? p->finH : "",  p->finL ? p->finL : "",  p->finP ? p->finP : "");

    /** read in the DEM file, initialize sparse matrices */
    void * ptrs[]= {p->mH,p->mL,p->vP};
    read_dem_file(p->fdem, ptrs, p->debug);
    p->mH=ptrs[0];
    p->mL=ptrs[1];
    p->vP=ptrs[2];
    //    csr_out(p->mH);
    p->nvar = p->mH->cols;
    p->nchk = p->mH->rows;
    p->ncws = p->mL->rows;
  }

  if(p->finH){
    if(! p->finL){
      if((p->fobs) || (p->fdet))
	ERROR("Without L matrix, cannot specify fdet=%s or fobs=%s\n",
	      p->fdet ? p->fdet : "",  p->fobs ? p->fobs : "");      
    }
    if((! p->finL) && (! p->finG)){ /** only `H` matrix specified */
      p->classical=1;
      p->internal=1;
    }
    p->mH=csr_mm_read(p->finH, p->mH, 0, p->debug);
    p->nvar = p->mH->cols;
    p->nchk = p->mH->rows;
  }
  /** at this point we should have `H` matrix for sure */
  assert(p->mH);

  if(p->finL){
    p->mL=csr_mm_read(p->finL, p->mL, 0, p->debug);
    p->ncws = p->mL->rows;
    if(p->mL->cols != p->nvar)
      ERROR("column number mismatch L[%d,%d] and H[%d,%d]\n",
	    p->mL->rows, p->mL->cols, p->mH->rows, p->mH->cols);
  }  

  if(p->finK){
    p->mK=csr_mm_read(p->finK, p->mK, 0, p->debug);
    int ncws = p->mK->rows;
    if(p->ncws){
      if (ncws != p->ncws)
	ERROR("number of codewords mismatch: K[%d,%d] and L[%d,%d]\n",
	      p->mK->rows, p->mK->cols, p->mL->rows, p->mL->cols);
    }
    else
      p->ncws = ncws;  
    if(p->mK->cols != p->nvar)
      ERROR("column number mismatch K[%d,%d] and H[%d,%d]\n",
	    p->mK->rows, p->mK->cols, p->mH->rows, p->mH->cols);
  }  

  if(p->finG){
    p->mG=csr_mm_read(p->finG, p->mG, 0, p->debug);
    if(p->mG->cols != p->nvar)
      ERROR("column number mismatch G[%d,%d] and H[%d,%d]\n",
	    p->mG->rows, p->mG->cols, p->mH->rows, p->mH->cols);

    /** verify row orthogonality of `G` and `H` */
    csr_t *Gt = csr_transpose(NULL, p->mG);
    mzd_t *MGt = mzd_from_csr(NULL,Gt); /* convert to dense matrix */
    csr_free(Gt);
    if((p->mH)&&(product_weight_csr_mzd(p->mH,MGt,0)))
      ERROR("rows of H=Hx and G=Hz should be orthogonal \n");

    if(!p->mL) /** create `Lx` */
      p->mL=Lx_for_CSS_code(p->mH,p->mG);
    p->ncws = p->mL->rows;

    /** verify row orthogonality of `G` and `L` */
    if((p->mL)&&(product_weight_csr_mzd(p->mL,MGt,0)))
      ERROR("rows of L=Lx and G=Hz should be orthogonal \n");
    mzd_free(MGt);        
  }

  if(p->useP > 0){/** override probability values */
    if(!p->vP){
      p->vP=malloc(p->nvar * sizeof(double));
      if(!p->vP)
	ERROR("memory allocation failed");
    }
    double *ptr=p->vP;
    for(int i=0;i<p->nvar;i++, ptr++)
      *ptr = p->useP;    
  }
  if(!p->vP)
    ERROR("probabilities missing, specify 'fdem', 'finP', or 'useP'");
    
  LLR_table = init_LLR_tables(p->d1,p->d2,p->d3);

  p->dE = llr_from_dbl(p->dEdbl);
  
  if(!p->mL){
    p->mL=csr_identity(p->nvar, p->nvar);
    p->ncws = p->nvar;
  }

  switch(p->mode){
  case 1:  /** currently only needed for BP */
    if(p->debug&1){
      out_LLR_params(LLR_table);
      printf("# submode=%d, %s BP using %s LLR\n",
	     p->submode,
	     p->submode & 4 ? (p->submode & 8 ? "serial-V" : "serial-C") : "parallel",
	     (((p->submode & 3)==3) || ((p->submode & 3)==0)) ? "both regular and average" :
	     p->submode & 2 ? "average" : "instantaneous");
      if(((p->submode & 4)==0) && (p->submode & 8))
	ERROR("mode=%d submode=%d invalid (add 4 to submode ->%d for serial-V)",p->mode,p->submode,p->submode | 4);
      if(((p->submode & 2)) || ((p->submode & 3)==0))
	printf("# use average LLR: aLLR = %g * aLLR + %g * LLR\n",p->bpalpha,1-p->bpalpha);
      if((p->submode & 4) && (p->submode & 16))
	printf("# randomize node order at each step");
    }
    if ((p->submode>=32))
      ERROR(" mode=%d BP : submode='%d' currently unsupported\n", p->mode, p->submode);    
    /* fall through */
  case 0:
    if((p->ferr) && (p->fobs))
      ERROR("With ferr=%s cannot specify fobs=%s (will calculate)\n",p->fdem, p->fobs);
    if((p->ferr) && (p->fdet))
      ERROR("With ferr=%s cannot specify fdet=%s (will calculate)\n",p->fdem, p->fdet);    
    if((p->fdet==NULL)&&(p->fobs!=NULL))
      ERROR(" mode=%d fobs='%s' need detection events file 'fdet'\n",
	    p->mode, p->fobs);
    if ((p->fdet!=NULL)&&(p->fobs==NULL))
      ERROR(" mode=%d fdet='%s' need observables file 'fobs'\n",
	    p->mode, p->fdet);

    if(p->fdet){/** expect both `fdet` and `fobs` to be defined */
      p->file_det=fopen(p->fdet, "r");
      if(p->file_det==NULL)
	ERROR("can't open the (det) file %s for reading\n",p->fdet);
      p->file_obs=fopen(p->fobs, "r");
      if(p->file_obs==NULL)
	ERROR("can't open the (obs) file %s for reading\n",p->fobs);
    }
    else if (p->ferr){
      p->file_err=fopen(p->ferr, "r");
      if(p->file_err==NULL)
	ERROR("can't open the (err) file %s for reading\n",p->ferr);
      p->internal=2;
    }
    else{
      p->internal=1;
    }
    if(p->ferr)
      p->internal=2; /* generate observables and detector events from errors provided */
    else if ((! p->fobs) && (! p->fdet))
      p->internal=1; /* generate observables and detector events internally */

    if ((p->mode==0)&&(p->submode!=0))
      ERROR(" mode=%d : non-zero submode='%d' unsupported\n",
	    p->mode, p->submode);
    
    break;
  case 2: /** estimate success probability */
    if((p->fdet!=NULL)||(p->fobs!=NULL) ||(p->ferr!=NULL))
      ERROR(" mode=%d, must not specify 'ferr' or 'fobs' or 'fdet' files\n",
	    p->mode);
    if (p->submode!=0)
      ERROR(" mode=%d : non-zero submode='%d' unsupported\n",
	    p->mode, p->submode);
    break;
    
  case 3: /** read in DEM file and output the H, L, G matrices and P vector */
    if(strcmp(p->fout,"stdout")==0)
      p->use_stdout=1;
    if (p->submode!=0)
      ERROR(" mode=%d : non-zero submode='%d' unsupported\n",
	    p->mode, p->submode);    
    break;
    
  default:
    ERROR(" mode=%d is currently not supported\n",p->mode);
    break;
  }

  if(p->debug & 64){ /** print matrices */
    assert(p->nvar != 0);
    mzd_t *mH0 = mzd_from_csr(NULL,p->mH);
    printf("matrix mH0:\n");  mzd_print(mH0);
    mzd_free(mH0);
    
    mzd_t *mL0 = mzd_from_csr(NULL,p->mL);
    printf("matrix mL0:\n");  mzd_print(mL0);
    mzd_free(mL0);

    //    for(int i=0; i < p->nvar; i++)
    //      printf(" P[%d]=%g \n",i,p->vP[i]);
  }
  
  
  if(p->fdet){/** expect both `fdet` and `fobs` to be defined */
    p->internal=0;
    p->file_det=fopen(p->fdet, "r");
    if(p->file_det==NULL)
      ERROR("can't open the (det) file %s for reading\n",p->fdet);
    p->file_obs=fopen(p->fobs, "r");
    if(p->file_obs==NULL)
      ERROR("can't open the (obs) file %s for reading\n",p->fobs);
  }

  if(p->mode<=1){ /* vecdec RIS or BP decoder */
    p->mHe = mzd_init(p->nchk, p->nvec); /** each column a syndrome vector `H*e` */
    p->mLe = mzd_init(p->ncws,  p->nvec); /** each column `L*e` vector */
    if(p->ferr)
      p->mE = mzd_init(p->nvar, p->nvec);
    if(p->debug &1){
      printf("# mode=%d, running %s decoder ",  p->mode, p->mode==0 ? "vecdec RIS" : "BP");
      switch(p->internal){
      case 0: printf("use DET and OBS files\n"); break;
      case 1: printf("internal error generator\n"); break;
      case 2: printf("reading error vectors from ERR file\n"); break;
      default: ERROR("this should not happen");
	break;
      }
    }

    /** initialize the counters */
    for(extr_t i=0; i< EXTR_MAX; i++){
      cnt[i]=0;
      iter1[i]=0;
      iter2[i]=0;
    }

  }

  return 0;
};

/** @brief clean up variables and open files */
void var_kill(params_t *p){
  if(p->file_det)  fclose(p->file_det);
  if(p->file_obs)  fclose(p->file_obs);  
  if(p->file_err)  fclose(p->file_err);  
  if(p->vP){        free(p->vP);    p->vP = NULL;  }
  if(p->vLLR){      free(p->vLLR);  p->vLLR = NULL;}
  if(LLR_table){ free(LLR_table);  LLR_table = NULL;}

  p->mH =  csr_free(p->mH);
  p->mHt = csr_free(p->mHt);
  p->mL =  csr_free(p->mL);
  p->mLt = csr_free(p->mLt);
  p->mG = csr_free(p->mG);
  if(p->mE) mzd_free(p->mE);
  if(p->mHe) mzd_free(p->mHe);
  if(p->mLe) mzd_free(p->mLe);
  if(p->mHeT) mzd_free(p->mHeT);
  if(p->mLeT) mzd_free(p->mLeT);

}

int do_err_vecs(params_t * const p){

  int il1=p->nvec, il2;
  /** prepare error vectors ************************************************************/
  switch(p->internal){
  case 0: /** read `det` and `obs` files */
    il1=read_01(p->mHe,p->file_det, &p->line_det, p->fdet, p->debug);
    /** TODO: enable external processing of observables */
    il2=read_01(p->mLe,p->file_obs, &p->line_obs, p->fobs, p->debug);
    if(il1!=il2)
      ERROR("mismatched DET %s (line %lld) and OBS %s (line %lld) files!",
	    p->fdet,p->line_det,p->fobs,p->line_obs);
    if(p->debug&1)
      printf("read %d det/obs pairs\n",il1);		 
    break;
  case 1: /** generate errors internally */
    do_errors(p->mHe,p->mLe,p->mHt, p->mLt, p->vP);
    if(p->debug&1)
      printf("# generated %d det/obs pairs\n",p->mHe->ncols);
    if((p->debug&128)&&(p->nvar <= 256)&&(p->nvec <= 256)&&(p->debug &512)){
      printf("He:\n");
      mzd_print(p->mHe);
      printf("Le:\n");
      mzd_print(p->mLe);
    }

    break;
  case 2: /** read errors from file and generate corresponding `obs` and `det` matrices */
    il1=read_01(p->mE,p->file_err, &p->line_err, p->ferr, p->debug);
    if(p->debug&1)
      printf("# read %d errors from file %s\n",il1,p->ferr);
    csr_mzd_mul(p->mHe,p->mH,p->mE,1);
    csr_mzd_mul(p->mLe,p->mL,p->mE,1);
    if((p->debug&128)&&(p->nvar <= 256)&&(p->nvec <= 256)&&(p->debug &512)){
      printf("error columns read:\n");
      mzd_print(p->mE);
      //      printf("He:\n");
      //      mzd_print(mHe);
      //      printf("Le:\n");
      //      mzd_print(mLe);
    }
    break;
  default:
    ERROR("internal=%d, this should not happen",p->internal);
  }
  return il1;
}


int main(int argc, char **argv){
  params_t * const p=&prm;
  
  /** initialize variables, read in the DEM file, initialize sparse matrices */
  var_init(argc,argv,  p);
  assert(p->mL);
    //    p->mL=csr_identity(p->nvar, p->nvar);
  init_Ht(p);
  //  mat_init(p);
  if(p->debug & 1)
    printf("# mode=%d submode=%d debug=%d\n",p->mode,p->submode,p->debug);

  long long int ierr_tot=0, rounds=(long long int )ceil((double) p->ntot / (double) p->nvec);
  if(((p->mode == 0) || (p->mode == 1)) && (p->debug & 2))
    printf("# ntot=%lld nvec=%d will do calculation in %lld rounds\n",p->ntot,p->nvec,rounds);

  switch(p->mode){
    long long int synd_tot, synd_fail; /** case 0 */
    qllr_t *ans;                  /** case 1 */
    size_t size;                  /** case 3 */
    char * comment;
  case 0: /** `mode=0` internal `vecdec` decoder */
    /** at least one round always */
    synd_tot=0, synd_fail=0;
    for(long long int iround=0; iround < rounds; iround++){
      if(p->debug &1)
	printf("# starting round %lld of %lld\n", iround, rounds);
    
      if( !(ierr_tot = do_err_vecs(p)))
	break; /** no more rounds */

      // actually decode and generate error vectors
      mzd_t *mE0=NULL;
#ifndef NDEBUG  /** need `mHe` later */
      mzd_t *mS=mzd_copy(NULL,p->mHe);
      mE0=do_decode(mS, p); /** each row a decoded error vector */
      mzd_free(mS); mS=NULL;
#else
      mE0=do_decode(p->mHe, p); /** each row a decoded error vector */
#endif /* NDEBUG */
      mzd_t *mE0t = mzd_transpose(NULL, mE0);
      mzd_free(mE0); mE0=NULL;
        
#ifndef NDEBUG
      mzd_t *prodHe = csr_mzd_mul(NULL,p->mH,mE0t,1);
      mzd_add(prodHe, prodHe, p->mHe);
      if(!mzd_is_zero(prodHe)){
	if((p->debug&512)||(p->nvec <=64)){
	  printf("syndromes difference:\n");
	  mzd_print(prodHe);
	}
	ERROR("some syndromes are not matched!\n");
      }
      mzd_free(prodHe); prodHe = NULL;
      //      mzd_free(p->mHe);    mHe    = NULL;
#endif

      mzd_t *prodLe = csr_mzd_mul(NULL,p->mL,mE0t,1);

      if(p->debug & 512){ /** print matrices */
	printf("prodLe:\n");
	mzd_print(prodLe);
	printf("mLe:\n");
	mzd_print(p->mLe);
      }

      mzd_add(prodLe, prodLe, p->mLe);
      //      mzd_free(mLe); mLe=NULL;

      int fails=0;
      for(rci_t ic=0; ic< ierr_tot; ic++){
	rci_t ir=0;
	if(mzd_find_pivot(prodLe, ir, ic, &ir, &ic)){
	  fails++;
	}
	else /** no more pivots */
	  break;
      }
      /** update the global counts */
      synd_tot  += ierr_tot; /** was: prodLe->ncols */
      synd_fail += fails;
      mzd_free(prodLe); prodLe=NULL;
      if((p->nfail > 0) && (synd_fail >= p->nfail))
	break;
    }
    if(p->debug&1)
      printf("# fail_fraction total_cnt sycces_cnt\n");
    printf(" %g %lld %lld # %s\n",(double) synd_fail/synd_tot, synd_tot, synd_tot-synd_fail, p->fdem ? p->fdem : p->finH );
      
    break;

  case 1: /** `mode=1` various BP flavors */    
    ans = calloc(p->nvar, sizeof(qllr_t));
      if(!ans) ERROR("memory allocation failed!"); 

    for(long long int iround=0; iround < rounds; iround++){
      if(p->debug &1)
	printf("# starting round %lld of %lld\n", iround, rounds);
    
      if( !(ierr_tot = do_err_vecs(p)))
	break; /** no more rounds */
      p->mHeT = mzd_transpose(p->mHeT,p->mHe);
      p->mLeT = mzd_transpose(p->mLeT,p->mLe);
      for(long long int ierr = 0; ierr < ierr_tot; ierr++){ /** cycle over errors */
	cnt[TOTAL]++;
	int succ_BP = 0, succ_OSD = 0;
	if(mzd_row_is_zero(p->mHeT,ierr)){
	  //	  printf("ierr=%d of tot=%d\n",ierr,ierr_tot);

	  cnt_update(CONV_TRIVIAL,0); /** trivial convergence after `0` steps */
	  if(mzd_row_is_zero(p->mLeT,ierr)){
	    cnt[SUCC_TRIVIAL]++;
	    cnt[SUCC_TOT]++;
	  }
	  if(p->debug & 128)
	    printf("error %lld of %lld is trivial\n",ierr+1,ierr_tot);
	  continue ; /** next error / syndrome vector pair */     
	}
	else{ /** non-trivial syndrome */
#ifndef NDEBUG	  
	  if((p->debug&8)&&(p->nvar <= 256)&&(p->debug &512)){
	    printf("non-trivial error %lld of %lld:\n",ierr+1,ierr_tot);
	    if(p->mE) /** print column as row */	      
	      for(int i=0; i<p->nvar; i++)
		printf("%s%d%s",i==0?"[":" ",mzd_read_bit(p->mE,i,ierr),i+1<p->nvar?"":"]\n");
	    mzd_print_row(p->mHeT,ierr);
	    mzd_print_row(p->mLeT,ierr);
	    out_llr("i",p->nvar,p->vLLR);
	  }
#endif 	  
	  mzd_t * const srow = mzd_init_window(p->mHeT, ierr,0, ierr+1,p->nchk); /* syndrome row */
	  if(p->submode&4){ /** bit 2 is set, use serial schedule */
	    if(p->submode&8)
	      succ_BP = do_serialV_BP(ans, srow, p->mH, p->mHt, p->vLLR, p);
	    else
	      succ_BP = do_serialC_BP(ans, srow, p->mH, p->mHt, p->vLLR, p);
	  }
	  else{             /** bit 2 is not set, parallel schedule */
	    succ_BP = do_parallel_BP(ans, srow, p->mH, p->mHt, p->vLLR, p);	   
	  }
	  if((!succ_BP) && (p->lerr >=0)){
	      if(p->debug&128)
		printf("ierr=%lld starting OSD lerr=%d maxosd=%d\n",ierr,p->lerr, p->maxosd);
	      do_osd_start(ans,srow,p->mH,p);
	      succ_OSD=1;
	    }

	  if((succ_BP)||(succ_OSD)){              /** `convergence success`  */
	    mzd_t * const obsrow = mzd_init_window(p->mLeT, ierr,0, ierr+1,p->mLeT->ncols);
	    if(syndrome_check(ans, obsrow, p->mL, p)){
	      cnt[SUCC_TOT]++;
	      if(succ_BP)
		cnt[SUCC_BP]++;
	      else
		cnt[SUCC_OSD]++;
	    }
	    mzd_free(obsrow);
	  }
	  mzd_free(srow);
	  
	  if(p->debug&16)
	    printf("i=%lld of %lld succ=%d\n",ierr,ierr_tot,succ_BP);
	}
	if((p->nfail) && cnt[TOTAL]-cnt[SUCC_TOT] >= p->nfail)
	  break;
      }
    }
    cnt_out(p->debug&1);
    free(ans);
    break;
  case 2: /** `mode=2` */
    if(p->debug&1)
      printf("# mode=%d, estimating fail probability in %d steps\n",p->mode, p->steps);
    if(p->finC){
      p->num_cws = nzlist_read(p->finC,p);
      if(p->debug&1)
	printf("# %lld codewords read from %s ...",p->num_cws, p->finC);
      do_hash_verify_CW(p->mHt, p->mLt, p);
      if(p->debug&1)
	printf("all verified\n");
      do_hash_min(p); /** set values `minE` and `minW` from stored CWs */
    }
    do_LLR_dist(p);
    do_hash_fail_prob(p); /** output results */
    if(p->outC){
      char * name;
      if(p->fdem)
	name=p->fdem;
      else if (p->finH)
	name=p->finH;
      else
	name="(unknown source)";
      size_t size = 1 + snprintf(NULL, 0, "codewords computed from '%s'", name);
      char *buf = malloc(size * sizeof(char));
      if(!buf)
	ERROR("memory allocation");
      snprintf(buf, size, "# codewords from '%s'", name);
      
      long long cnt=nzlist_write(p->outC, buf, p);
      if(p->debug & 1)
	printf("wrote %lld computed codewords to file %s\n",cnt,p->outC);
    }
    do_hash_clear(p);
    break;
    
  case 3: /** read in DEM file and output the H, L, G matrices and P vector */
    size = snprintf(NULL, 0, "H matrix from DEM file %s", p->fdem);
    if(!(comment = malloc(size + 1)))
      ERROR("memory allocation");
    sprintf(comment, "H matrix from DEM file %s", p->fdem);
    if(p->debug&1)
      printf("# writing H=Hx matrix [ %d x %d ] to \t%s%s\n",
	     p->mH->rows, p->mH->cols, p->fout, p->use_stdout ? "\n" :"H.mmx");
    csr_mm_write(p->fout,"H.mmx",p->mH,comment);
    if(p->debug&1)
      printf("# writing L=Lx matrix [ %d x %d ] to \t%s%s\n",
	     p->mL->rows, p->mL->cols, p->fout, p->use_stdout ? "\n" :"L.mmx");
    comment[0]='L';
    csr_mm_write(p->fout,"L.mmx",p->mL,comment);

    if(p->debug&1)
      printf("# writing P vector [ %d ] to      \t%s%s\n",
	     p->nvar, p->fout, p->use_stdout ? "\n" :"P.mmx");
    comment[0]='P';
    //    printf("%% %s\n", comment);
    dbl_mm_write(p->fout,"P.mmx",1,p->nvar,p->vP,comment);

    if(!p->mG){
      if(p->debug&1)
	printf("# creating G=Hz matrix and writing to\t%s%s\n",
	       p->fout, p->use_stdout ? "\n" :"G.mmx");
      p->mG = do_G_matrix(p->mHt,p->mLt,p->vLLR, p->debug);
      comment[0]='G';
    }
    else{
      size_t size2 = snprintf(NULL, 0, "G matrix from file %s", p->finG);
      if(size2>size){
	comment = realloc(comment, size2 + 1);
	size=size2;
      }    
      assert(comment);
      sprintf(comment, "G matrix from file %s", p->finG);
    }
    csr_mm_write(p->fout,"G.mmx",p->mG,comment);

    if(!p->mK){ /** create `Lz` */
      if(p->debug&1)
	printf("# creating K=Lz matrix and writing to\t%s%s\n",
	       p->fout, p->use_stdout ? "\n" :"K.mmx");
      p->mK=Lx_for_CSS_code(p->mG,p->mH);
      comment[0]='K';
    }
    else{
      size_t size2 = snprintf(NULL, 0, "K matrix from file %s", p->finK);
      if(size2>size)
	if(!(comment = realloc(comment, size2 + 1)))
	  ERROR("memory allocation");
      sprintf(comment, "K matrix from file %s", p->finK);
    }
    csr_mm_write(p->fout,"K.mmx",p->mK,comment);
    
    free(comment);
    break;
    
  default:
    ERROR("mode=%d not supported\n",p->mode);
    break;
  }

  var_kill(p);
  return 0;
}
