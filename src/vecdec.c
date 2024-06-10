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
#include <stdbool.h>

params_t prm={ .nchk=-1, .nvar=-1, .ncws=-1, .steps=50, .pads=0,
  .rankH=0, .rankG=-1, .rankL=-1, 
  .lerr=-1, .maxosd=100, .bpalpha=0.5, .bpretry=1, .swait=0, .maxC=0,
  .dW=0, .minW=INT_MAX, .maxW=0, .dE=-1, .dEdbl=-1, .minE=INT_MAX,
  .nvec=1024, .ntot=1, .nfail=0, .seed=0, .epsilon=1e-8, .useP=0, .dmin=0,
  .debug=1, .fdem=NULL, .fdet=NULL, .fobs=NULL, .fout="tmp", .ferr=NULL,
  .mode=-1, .submode=0, .use_stdout=0, 
  .LLRmin=0, .LLRmax=0, .codewords=NULL, .num_cws=0,
  .finH=NULL, .finL=NULL, .finG=NULL, .finK=NULL, .finP=NULL,
  .finC=NULL, .outC=NULL, 
  .vP=NULL, .vLLR=NULL, .mH=NULL, .mHt=NULL,
  .mL=NULL, .mLt=NULL, .internal=0, 
  .file_det=NULL, .file_obs=NULL, .line_det=0, .line_obs=0,
  .mE=NULL, .mHe=NULL, .mLe=NULL, .mHeT=NULL, .mLeT=NULL,
  .nzH=0, .nzL=0
};

params_t prm_default={  .steps=50, .pads=0, 
  .lerr=-1, .maxosd=100, .bpalpha=0.5, .bpretry=1, .swait=0, .maxC=0,
  .dW=0, .minW=INT_MAX, .maxW=0, .dE=-1, .dEdbl=-1, .minE=INT_MAX,
  .nvec=1024, .ntot=1, .nfail=0, .seed=0, .epsilon=1e-8, .useP=0, .dmin=0,
  .debug=1, .fout="tmp", .ferr=NULL,
  .mode=-1, .submode=0, .use_stdout=0, 
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

/** @brief construct `L` matrix for a classical code orthogonal to rows of `H`
 * 
 * The matrix should have `k` rows of weight `1` each.
 */
csr_t * do_L_classical(const csr_t * const H, const params_t * const p){
  assert((H!=NULL) && (H->nz == -1));
  mzd_t * mH = mzd_from_csr(NULL, H);

  /** construct a list of non-pivot columns in `H` */
  mzp_t * skip_pivs=mzp_init(p->nvar); 
  if(!skip_pivs)
    ERROR("memory allocation failed!\n");
  int rank=0, k=0;
  for(int col=0; col< p->nvar; col++){
    int ret=gauss_one(mH, col, rank);
    if(ret)
      rank++;
    else
      skip_pivs->values[k++]=col;
  }
  mzd_free(mH);

  /** construct the matrix to return */
  csr_t *ans = csr_init(NULL,k,p->nvar,k);
  for(int i=0; i<k; i++){
    ans->p[i] = i;
    ans->i[i] = skip_pivs->values[i];
  }
  ans->p[k]=k;
  ans->nz = -1;

  mzp_free(skip_pivs);
  return ans;
}

/** @brief Read in vectors from `NZLIST` file `fnam` to hash.
    @return number of entries read */
long long int nzlist_read(const char fnam[], params_t *p){
  long long int count = 0, lineno;
  assert(fnam);
  FILE * f=nzlist_r_open(fnam, &lineno); /** no need to check if open */
  one_vec_t *entry=NULL;
  while((entry=nzlist_r_one(f,NULL, fnam, &lineno))){
    if((p->maxC) && (count >= p->maxC))
      ERROR("number of entries in file %s exceeds maxC=%lld\n", p->finC, p->maxC);
    if((p->maxW==0) ||((p->maxW) && (entry->weight <= p->maxW))){
      const size_t keylen = entry->weight * sizeof(rci_t);
      /** `DO NOT assume` unique vectors stored in the file */
      one_vec_t *pvec=NULL;
      HASH_FIND(hh, p->codewords, entry->arr, keylen, pvec);
      if(!pvec){ /** vector not found, inserting */
	
	qllr_t energ=0;
	for(int i=0; i < entry->weight; i++) 
	  energ += p->vLLR[entry -> arr[i]];
	entry->energ=energ;
	
	HASH_ADD(hh, p->codewords, arr, keylen, entry); /** store in the `hash` */
	count ++;
	if(p->minE > entry->energ)
	  p->minE = entry->energ;

	if(p->minW > entry->weight)
	  p->minW = entry->weight;
      }
    }
  }
  p->num_cws += count; 
  fclose(f);
  return count; 
}

/** @brief Write vectors from hash to `NZLIST` file `fnam`.
    @return number of entries written */
long long int nzlist_write(const char fnam[], const char comment[], params_t *p){
  long long int count=0;
  assert(fnam);
  FILE * f = nzlist_w_new(fnam, comment); /** no need to check if open */
  one_vec_t *pvec;
  
  for(pvec = p->codewords; pvec != NULL; pvec = (one_vec_t *)(pvec->hh.next)){
    if((p->maxW==0) ||((p->maxW) && (pvec->weight <= p->maxW))){
      count ++;
      nzlist_w_append(f,pvec);
    }
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
  if (k)
    assert(mHt->rows == mLt->rows);
  //  else    assert(p->classical);
  
  mzd_t *vHt = mzd_init(1, r);
  mzd_t *vLt = k>0 ? mzd_init(1, k) : NULL;
  one_vec_t *pvec;
  long long int count=0;

  for(pvec = p->codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
    int not_cw = 0; 
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
      ERROR("invalid hash vector [%lld]\n",count);
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
      if ((p->maxW==0) || ((p->maxW != 0) && (pvec->weight <= p->maxW)))
	print_one_vec(pvec);
  }
  /** finally calculate and output fail probability here */
  /** TODO: use the limit on `W` and `E` (just ignore codewords outside the limit) */
  double pfail=0, pmax=0;
  int minW=p->nvar + 1;
  for(pvec = p->codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
    if ((p->maxW==0) || ((p->maxW != 0) && (pvec->weight <= p->maxW))){
      double prob=do_prob_one_vec(pvec, p);
      pfail += prob;
      if(prob>pmax)
	pmax=prob;
      if(minW> pvec->weight)
	minW=pvec->weight;
    }
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
 * @param classical set to `1` for classical code (do not use `L` matrix), `0` otherwise  
 * @return minimum `weight` of a CW found (or `-weigt` if early termination condition is reached). 
 */
int do_LLR_dist(params_t  * const p, const int classical){
  const int dW=p->dW;
  //  qllr_t dE=p->dE;
  /** whether to verify logical ops as a vector or individually */
  int use_vector=0;
  if (!classical)
    use_vector = p->mLt->cols >= 16 ? 1 : 0;
  if(p->nvec == 16) /** default value */
    p->nvec=0;
  mzd_t * mH = mzd_from_csr(NULL, p->mH);
  mzd_t *mLt = NULL, *eemLt = NULL, *mL = NULL;
  if(!classical){
    if(use_vector){ /** all logical ops at once */
      mLt = mzd_from_csr(NULL, p->mLt);
      eemLt = mzd_init(1, p->mLt->cols);
    }
    else /** just one logical op */
      mL = mzd_from_csr(NULL, p->mL);
  }
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
      if (classical)
	nz=1; /** no need to verify */
      else{
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
    if(eemLt)
      mzd_free(eemLt);
    if (mLt)
      mzd_free(mLt);
  }
  else
    if (mL)
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
    //    if (p->lerr > 1){
    if (mEt1){
      mzd_free(mEt1);
      mEt1=NULL;
    }
    if (vE1){
      free(vE1);
      vE1=NULL;
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
      //      if (p->lerr > 1)
      if(vE1){
	free(vE1);
	vE1=NULL;
      }
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

/**
 * @brief Calculate the energy change when a row from mH is combined with a row from mEt0.
 * @param coeff Coefficients array for energy calculation.
 * @param A The matrix containing current error vectors.
 * @param H The parity-check matrix.
 * @param i The index of the row in A to be updated.
 * @param row_to_add The index of the row in H to add to the i-th row of A.
 * @return The change in energy.
 */
qllr_t calculate_energy_change(qllr_t *coeff, const mzd_t *A, const mzd_t *B, const int i, const int row_to_add) {
    qllr_t energy_change = 0;
    // Loop over all columns to calculate the change in energy due to the row addition.
    for (rci_t j = 0; j < A->ncols; ++j) {
        // Check if there is a change in the bit after addition.
        // If the bits in A and B are different, the result is 1 (change), otherwise it's 0 (no change).
        if (mzd_read_bit(A, i, j) != mzd_read_bit(B, row_to_add, j)) {
            // If the original bit in A is 1, the energy will decrease (because the bit flips to 0).
            // If the original bit in A is 0, the energy will increase (because the bit flips to 1).
            energy_change += (mzd_read_bit(A, i, j) ? -coeff[j] : coeff[j]);
        }
    }
    return energy_change;
}

/**
 * Decompress a CSR matrix if it is in compressed form.
 * @param mat The CSR matrix to decompress.
 */
void csr_decompress(csr_t *mat) {
    if (mat == NULL || mat->nz != -1) {
        return;  // Matrix is already decompressed or invalid
    }

    int nz = mat->p[mat->rows];
    int *new_i = malloc(nz * sizeof(int));
    if (new_i == NULL) {
        ERROR("Memory allocation failed in csr_decompress.\n");
        return;
    }

    int idx = 0;
    for (int i = 0; i < mat->rows; i++) {
        for (int j = mat->p[i]; j < mat->p[i + 1]; j++) {
            new_i[idx++] = mat->i[j];
        }
    }

    free(mat->i);
    mat->i = new_i;
    mat->nz = nz;
}


/**
 * @brief Copy a CSR matrix
 * @param src The source CSR matrix to copy from.
 * @return A new CSR matrix that is a copy of the source.
 */
csr_t *csr_copy(const csr_t *src) {
    if (src == NULL) {
        return NULL;
    }

    csr_t *dst = malloc(sizeof(csr_t));
    if (dst == NULL) {
        ERROR("Memory allocation failed for CSR matrix copy.\n");
        return NULL;
    }

    dst->rows = src->rows;
    dst->cols = src->cols;
    dst->nzmax = src->nzmax;

    dst->p = malloc((src->rows + 1) * sizeof(int));
    if (dst->p == NULL) {
        free(dst);
        ERROR("Memory allocation failed for CSR matrix data.\n");
        return NULL;
    }

    dst->i = malloc(src->nzmax * sizeof(int));
    if (dst->i == NULL) {
        free(dst->p);
        free(dst);
        ERROR("Memory allocation failed for CSR matrix data.\n");
        return NULL;
    }

    memcpy(dst->p, src->p, (src->rows + 1) * sizeof(int));
    memcpy(dst->i, src->i, src->nzmax * sizeof(int));

    dst->nz = src->nz;

    return dst;
}



/**
 * @brief Add row j of matrix B to row i of matrix A, and store the result in row i of matrix A.
 * @param A The CSR matrix to be updated.
 * @param i The row index in matrix A to be updated.
 * @param B The CSR matrix providing the row to be added.
 * @param j The row index in matrix B to be added to row i of matrix A.
 */
void csr_add_row_to_row(csr_t *A, int i, const csr_t *B, int j) {
    if (A == NULL || B == NULL || i >= A->rows || j >= B->rows) {
        ERROR("Invalid input parameters in csr_add_row_to_row.\n");
        return;
    }

    if (A->nz != -1 || B->nz != -1) {
        ERROR("Matrices must be in compressed form for csr_add_row_to_row.\n");
        return;
    }

    int start_A = A->p[i];
    int end_A = A->p[i + 1];
    int start_B = B->p[j];
    int end_B = B->p[j + 1];

    int *toggle_columns = calloc(A->cols, sizeof(int));
    if (toggle_columns == NULL) {
        ERROR("Memory allocation failed in csr_add_row_to_row.\n");
        return;
    }

    for (int idx_A = start_A; idx_A < end_A; ++idx_A) {
        toggle_columns[A->i[idx_A]] ^= 1;
    }
    for (int idx_B = start_B; idx_B < end_B; ++idx_B) {
        toggle_columns[B->i[idx_B]] ^= 1;
    }

    int nz_required = 0;
    for (int col = 0; col < A->cols; col++) {
        if (toggle_columns[col]) {
            nz_required++;
        }
    }

    int old_nz = end_A - start_A;
    int total_nonzeros = A->p[A->rows];
    int new_nzmax = total_nonzeros - old_nz + nz_required;

    if (new_nzmax > A->nzmax) {
        int *new_i = realloc(A->i, new_nzmax * sizeof(int));
        if (new_i == NULL) {
            free(toggle_columns);
            ERROR("Memory reallocation failed in csr_add_row_to_row.\n");
            return;
        }
        A->i = new_i;
        A->nzmax = new_nzmax;
    }

    if (nz_required != old_nz) {
        memmove(&A->i[start_A + nz_required], &A->i[end_A], (total_nonzeros - end_A) * sizeof(int));
    }

    int idx = 0;
    for (int col = 0; col < A->cols; col++) {
        if (toggle_columns[col]) {
            A->i[start_A + idx] = col;
            idx++;
        }
    }

    int new_nz = idx;
    //A->p[i + 1] = A->p[i] + new_nz;
    for (int k = i + 1; k <= A->rows; k++) {
        A->p[k] += new_nz - old_nz;
    }

    free(toggle_columns);
}






/**
 * @brief Calculate the energy of row `i` in CSR matrix.
 * @param coeff Coefficients array for energy calculation.
 * @param csr The CSR matrix.
 * @param i The row index.
 * @return The energy of the row.
 */
qllr_t csr_row_energ(qllr_t *coeff, const csr_t *csr, const int i) {
    if (coeff == NULL || csr == NULL || i >= csr->rows || i < 0) {
        ERROR("Invalid input parameters in csr_row_energ: coeff=%p, csr=%p, i=%d, rows=%d\n",
              coeff, csr, i, (csr ? csr->rows : -1));
        return 0;
    }

    qllr_t ans = 0;

    if (csr->nz != -1) {
        ERROR("Matrix must be in compressed form for csr_row_energ.\n");
        return 0;
    }

    int start = csr->p[i];
    int end = csr->p[i + 1];

    if (start < 0 || end > csr->p[csr->rows] || start > end) {
        ERROR("Invalid row pointer range in csr_row_energ: start=%d, end=%d, nz=%d\n", start, end, csr->nz);
        return 0;
    }

    for (int j = start; j < end; ++j) {
        int col = csr->i[j];
        if (col < 0 || col >= csr->cols) {
            ERROR("Invalid column index in csr_row_energ: col=%d, cols=%d\n", col, csr->cols);
            return 0;
        }
        ans += coeff[col];
    }

    return ans;
}




/**
 * @brief Calculate the energy change when a row from csr B is combined with a row from csr A.
 * @param coeff Coefficients array for energy calculation.
 * @param A The CSR matrix containing current error vectors.
 * @param B The CSR matrix containing the rows to add.
 * @param i The index of the row in A to be updated.
 * @param row_to_add The index of the row in B to add to the i-th row of A.
 * @return The change in energy.
 */
qllr_t csr_calculate_energy_change(qllr_t *coeff, const csr_t *A, const csr_t *B, const int i, const int row_to_add) {
    if (coeff == NULL || A == NULL || B == NULL || i >= A->rows || row_to_add >= B->rows || i < 0 || row_to_add < 0) {
        ERROR("Invalid input parameters in csr_calculate_energy_change: coeff=%p, A=%p, B=%p, i=%d, row_to_add=%d\n",
              coeff, A, B, i, row_to_add);
        return 0;
    }

    qllr_t energy_change = 0;

    if (A->nz != -1 || B->nz != -1) {
        ERROR("Matrices must be in compressed form for csr_calculate_energy_change.\n");
        return 0;
    }

    int start_A = A->p[i];
    int end_A = A->p[i + 1];
    int start_B = B->p[row_to_add];
    int end_B = B->p[row_to_add + 1];

    int *toggle_columns = calloc(A->cols, sizeof(int));
    if (toggle_columns == NULL) {
        ERROR("Memory allocation failed in csr_calculate_energy_change.\n");
        return 0;
    }

    for (int j = start_A; j < end_A; ++j) {
        toggle_columns[A->i[j]] ^= 1;
    }

    for (int j = start_B; j < end_B; ++j) {
        int col = B->i[j];
        if (toggle_columns[col]) {
            energy_change -= coeff[col];
        } else {
            energy_change += coeff[col];
        }
        toggle_columns[col] ^= 1;
    }

    free(toggle_columns);
    return energy_change;
}










/**
 * @brief Replace the row i of destination CSR matrix with the row i of source CSR matrix.
 * @param dest The destination CSR matrix.
 * @param src The source CSR matrix.
 * @param i The row index to replace.
 */
void csr_replace_row(csr_t *dest, const csr_t *src, int i) {
    if (dest == NULL || src == NULL || i >= dest->rows || i >= src->rows) {
        ERROR("Invalid input parameters in csr_replace_row.\n");
        return;
    }

    if (dest->nz != -1) {
        csr_compress(dest);
    }
    if (src->nz != -1) {
        csr_compress((csr_t *)src);
    }

    int start_src = src->p[i];
    int end_src = src->p[i + 1];
    int start_dest = dest->p[i];
    int end_dest = dest->p[i + 1];
    int nz_src = end_src - start_src;
    int nz_dest = end_dest - start_dest;

    int new_nz = dest->p[dest->rows] - nz_dest + nz_src;
    if (dest->nzmax < new_nz) {
        dest->nzmax = new_nz * 2;
        dest->i = realloc(dest->i, dest->nzmax * sizeof(int));
        if (dest->i == NULL) {
            ERROR("Memory reallocation failed in csr_replace_row.\n");
            return;
        }
    }

    if (nz_src != nz_dest) {
        memmove(&dest->i[start_dest + nz_src], &dest->i[end_dest], (dest->p[dest->rows] - end_dest) * sizeof(int));
    }

    memcpy(&dest->i[start_dest], &src->i[start_src], nz_src * sizeof(int));

    for (int j = i + 1; j <= dest->rows; j++) {
        dest->p[j] += nz_src - nz_dest;
    }

    dest->nz = new_nz;
}






/**
 * @brief Executes simulated annealing to adjust error matrix based on energy changes.
 * @param mEt0_csr The initial error matrix in CSR format to be updated.
 * @param mHz CSR matrix representing the syndrome structure.
 * @param mEt_csr CSR matrix where the results are stored.
 * @param p Structure containing model parameters and error probabilities.
 */
void simulate_annealing_mcmc(csr_t *mEt0_csr, params_t const * const p,int i) {
    int max_iterations = 5000;
    double final_temperature = 1.0;
    double cooling_rate = 1.0;
    double reheating_rate = 1.0;
    int check_interval = 100;
    double min_acceptance_rate = 0.2;

    //double temperature = csr_row_energ(p->vLLR,mEt0_csr,i);
    double temperature = 1.0;
    int iterations = 0;
    int accepted_updates = 0;

    //while (iterations < max_iterations && temperature >= final_temperature) {
    while (iterations < max_iterations && temperature >= final_temperature) {
            iterations++;
            unsigned int seed = time(NULL);  // Use the current time as seed for randomness
            int row_to_add = rand_r(&seed) % (p->mG->rows); //here no fault
            
            
            double delta_energy = csr_calculate_energy_change(p->vLLR, mEt0_csr, p->mG, i, row_to_add);

            if (delta_energy < 0 || exp(-delta_energy / temperature) > (double)rand_r(&seed) / RAND_MAX) {
                csr_add_row_to_row(mEt0_csr, i, p->mG, row_to_add);
                accepted_updates++;
            }

            if (iterations % check_interval == 0) {
                double acceptance_rate = (double)accepted_updates / iterations;
                if (acceptance_rate < min_acceptance_rate) {
                    temperature *= reheating_rate;
                }
            }
            temperature *= cooling_rate;
        }
    //double energ=csr_row_energ(p->vLLR,mEt0_csr,i);
    //printf("energy=%f",energ);
}

/**
 * @brief Solves the equation to determine which ensemble has a lower free energy.
 * @param avg_exp_U_A_minus_U_B_A The average of <exp(U_A - U_B)>_A.
 * @param avg_exp_U_A_minus_U_B_B The average of <exp(U_A - U_B)>_B.
 * @param N_A The number of MCMC steps in ensemble A.
 * @param N_B The number of MCMC steps in ensemble B.
 * @return 1 if ensemble A is preferred, -1 if ensemble B is preferred.
 */
int solve_free_energy(double avg_exp_U_A_minus_U_B_A, double avg_exp_U_A_minus_U_B_B, int N_A, int N_B) {
    double x = 1.0;
    double tolerance = 1e-10;
    int max_iterations = 1000;

    for (int iter = 0; iter < max_iterations; ++iter) {
        double f_x = x - (1 + avg_exp_U_A_minus_U_B_A * N_B / N_A * x) /
                           ((avg_exp_U_A_minus_U_B_A * x) * (1 + avg_exp_U_A_minus_U_B_B * N_B / N_A * x));
        double f_x_prime = 1 - ((avg_exp_U_A_minus_U_B_A * N_B / N_A * (1 + avg_exp_U_A_minus_U_B_B * N_B / N_A * x)) -
                               (1 + avg_exp_U_A_minus_U_B_A * N_B / N_A * x) * avg_exp_U_A_minus_U_B_A * 
                               (1 + avg_exp_U_A_minus_U_B_B * N_B / N_A * x)) /
                               pow((avg_exp_U_A_minus_U_B_A * x * (1 + avg_exp_U_A_minus_U_B_B * N_B / N_A * x)), 2);
        double x_new = x - f_x / f_x_prime;

        if (fabs(x_new - x) < tolerance) {
            x = x_new;
            break;
        }
        x = x_new;
    }

    if (x >= 1.0) {
        return 1;
    } else {
        return -1;
    }
}


void BAR_MCMC(csr_t *mEt0_csr, csr_t *mEt_csr, params_t const * const p, int i) {
    int max_iterations = 5000;
    double temperature = 1.0;
    int num_ensembles = p->mK->rows + 1;

    // Create arrays to store flipped error vectors for each Lz row and the original error vector
    csr_t **error_vectors = malloc(num_ensembles * sizeof(csr_t *));
    // if (error_vectors == NULL) {
    //     ERROR("Memory allocation failed for error_vectors array.\n");
    //     return;
    // }
    //printf("Malloc done\n");

    // Initialize error vectors
    error_vectors[0] = csr_copy(mEt0_csr);

    // if (error_vectors[0] == NULL) {
    //     ERROR("Memory allocation failed for error_vectors[0].\n");
    //     //free(error_vectors);
    //     return;
    // }
    for (int k = 1; k < num_ensembles; k++) {
        error_vectors[k] = csr_copy(mEt0_csr);
        
        // if (error_vectors[k] == NULL) {
        //     ERROR("Memory allocation failed for error_vectors[%d].\n", k);
        //     for (int j = 0; j < k; j++) {
        //         csr_free(error_vectors[j]);
        //     }
        //     //free(error_vectors);
        //     return;
        // }
        csr_add_row_to_row(error_vectors[k], i, p->mK, k - 1);
    }
    //printf("Initialize ensembles done\n");

    // Variables for accumulating statistics
    double exp_U_A_minus_U_B_A = 0.0;
    double exp_U_A_minus_U_B_B = 0.0;
    int N_A = 0;
    int N_B = 0;
    double U_A_A = 0;
    double U_B_A = 0;
    double U_A_B = 0;
    double U_B_B = 0;
    csr_t *current_ensemble = csr_copy(error_vectors[0]);
    // if (current_ensemble == NULL) {
    //     ERROR("Memory allocation failed for current_ensemble.\n");
    //     for (int k = 0; k < num_ensembles; k++) {
    //         csr_free(error_vectors[k]);
    //     }
    //     //free(error_vectors);
    //     return;
    // }
    //printf("Current ensemble init done\n");

    // Iterate through each ensemble starting from the second one
    for (int k = 1; k < num_ensembles; k++) {
        exp_U_A_minus_U_B_A = 0.0;
        exp_U_A_minus_U_B_B = 0.0;
        N_A = 0;
        N_B = 0;

        // Perform MCMC in the current ensemble
        for (int iterations = 0; iterations < max_iterations; iterations++) {
            unsigned int seed = time(NULL);
            int row_to_add = rand_r(&seed) % (p->mG->rows);
            // if (error_vectors[k] == NULL || current_ensemble == NULL) {
            //     ERROR("NULL pointer in BAR_MCMC at iteration %d, ensemble %d\n", iterations, k);
            //     return;
            // }

            // Perform MCMC step for the current ensemble
            double delta_energy_A = csr_calculate_energy_change(p->vLLR, current_ensemble, p->mG, i, row_to_add);
            //printf("delta_A\n");
            if (delta_energy_A < 0 || exp(-delta_energy_A / temperature) > (double)rand_r(&seed) / RAND_MAX) {
                csr_add_row_to_row(current_ensemble, i, p->mG, row_to_add);
                csr_add_row_to_row(error_vectors[k], i, p->mG, row_to_add);
            }
            //printf("Add_done_A\n");
            U_A_A = csr_row_energ(p->vLLR, current_ensemble, i);// csr_row_energ has memory leakage
            //printf("Energ_done_A1\n");
            U_B_A = csr_row_energ(p->vLLR, error_vectors[k], i);// csr_row_energ has memory leakage
            //printf("Energ_done_A2\n");
            exp_U_A_minus_U_B_A += exp(U_A_A - U_B_A);
            //printf("exp_A\n");
            N_A++;
        }
        //printf("Current ensemble done\n");

        // Perform MCMC in the k-th ensemble
        for (int iterations = 0; iterations < max_iterations; iterations++) {
            unsigned int seed = time(NULL);
            int row_to_add = rand_r(&seed) % (p->mG->rows);
            // if (error_vectors[k] == NULL || current_ensemble == NULL) {
            //     ERROR("NULL pointer in BAR_MCMC at iteration %d, ensemble %d\n", iterations, k);
            //     return;
            // }

            // Perform MCMC step for the k-th ensemble
            double delta_energy_B = csr_calculate_energy_change(p->vLLR, error_vectors[k], p->mG, i, row_to_add);
            //printf("delta_B\n");
            if (delta_energy_B < 0 || exp(-delta_energy_B / temperature) > (double)rand_r(&seed) / RAND_MAX) {
                csr_add_row_to_row(error_vectors[k], i, p->mG, row_to_add);
                csr_add_row_to_row(current_ensemble, i, p->mG, row_to_add);
            }
            //printf("Add_done_B\n");
            U_A_B = csr_row_energ(p->vLLR, current_ensemble, i);
            //printf("Energ_done_B1\n");
            U_B_B = csr_row_energ(p->vLLR, error_vectors[k], i);
            //printf("Energ_done_B2\n");
            exp_U_A_minus_U_B_B += exp(U_A_B - U_B_B);
            //printf("exp_B\n");
            N_B++;
        }
        //printf("k-th ensemble done\n");

        double average_exp_U_A_minus_U_B_A = exp_U_A_minus_U_B_A / N_A;
        double average_exp_U_A_minus_U_B_B = exp_U_A_minus_U_B_B / N_B;
        int result = solve_free_energy(average_exp_U_A_minus_U_B_A, average_exp_U_A_minus_U_B_B, N_A, N_B);

        if (result == -1) {
            csr_free(current_ensemble);
            current_ensemble = csr_copy(error_vectors[k]);
            // if (current_ensemble == NULL) {
            //     ERROR("Memory allocation failed for current_ensemble during update.\n");
            //     for (int k = 0; k < num_ensembles; k++) {
            //         csr_free(error_vectors[k]);
            //     }
            //     free(error_vectors);
            //     return;
            // }
        }
        //printf("Update current ensemble\n");
    }
    //printf("BAR done\n");

    csr_replace_row(mEt_csr, current_ensemble, i);
    //printf("Replace done\n");

    for (int k = 0; k < num_ensembles; k++) {
        //printf("Freeing error_vector[%d]\n", k);
        csr_free(error_vectors[k]);
    }
    free(error_vectors);
    csr_free(current_ensemble);
    //printf("Free done\n");
}


/**
 * @brief Simulated Annealing-based syndrome decoding routine with improved stopping mechanism. (CSR form)
 * @param mS Matrix with syndromes (each column).
 * @param p Structure with error model information.
 * @return Binary matrix of min weight errors for each syndrome.
 */
mzd_t *do_decode_BAR(mzd_t *mS, params_t const * const p) {
    mzd_t *mHx_dense = mzd_from_csr(NULL, p->mH);
    mzd_t *mE_dense = mzd_init(mHx_dense->ncols, mS->ncols);

    mzp_t *perm = mzp_init(p->nvar);
    mzp_t *pivs = mzp_init(p->nvar);

    perm = sort_by_llr(perm, p->vLLR, p);
    int rank = 0;
    for (int i = 0; i < p->nvar; i++) {
        int col = perm->values[i];
        int ret = twomat_gauss_one(mHx_dense, mS, col, rank);
        if (ret)
            pivs->values[rank++] = col;
    }

    mzd_set_ui(mE_dense, 0);
    for (int i = 0; i < rank; i++) {
        mzd_copy_row(mE_dense, pivs->values[i], mS, i);
    }
    mzp_free(perm);
    mzp_free(pivs);
    //printf("rows of G matrix:%d\n",p->mG->rows);
    csr_t *mE_csr = csr_from_mzd(NULL,mE_dense);
    csr_t *mEt0_csr = csr_transpose(NULL, mE_csr);

    //printf("nz of mG=%d\n", p->mG->nz);

    int nrows_E = mHx_dense->ncols;
    int ncols_E = mS->ncols;
    
    int nzmax_E = nrows_E * ncols_E;
    csr_t *mEt_csr = csr_init(NULL, ncols_E, nrows_E, nzmax_E);
    //csr_t *mEt_csr=csr_copy(mEt0_csr);
    if (mEt_csr == NULL) {
        ERROR("Memory allocation failed for CSR structure.\n");
        return NULL;
    }
    csr_compress(mEt_csr);
    //#pragma omp parallel for
    //int tracker = 1;
    for (int i = 0; i < mEt0_csr->rows; i++) {
        //simulate_annealing_mcmc(mEt0_csr, p, i);
        BAR_MCMC(mEt0_csr,mEt_csr,p,i);
        //tracker++;
        //printf("iter=%d\n", tracker);
    }

    //printf("finished calculation");

    mzd_t *mEt = mzd_from_csr(NULL, mEt_csr);
    return mEt;
}



/** @brief one more matrix initialization routine 
 */
void init_Ht(params_t *p){
  const int n = p->nvar;
  p->mHt = csr_transpose(p->mHt, p->mH);
  if(p->mL)
    p->mLt = csr_transpose(p->mLt,p->mL);
  /** todo: fix reallocation logic to be able to reuse the pointers model */
  //  if(p->vP)    free(p->vP);
  //  p->vP=inP;
  if(!(p->vLLR = malloc(n*sizeof(qllr_t))))
    ERROR("memory allocation!");
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

    if(p->mL){
      printf("mL:\n");
      //    csr_out(mL);
      mzd_t *mdL = mzd_from_csr(NULL,p->mL);
      mzd_print(mdL);
      mzd_free(mdL);
      //    }
    }
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
    else if (sscanf(argv[i],"pads=%d",&dbg)==1){ /** `pads` */
      p -> pads = dbg;
      if (p->debug&1)
	printf("# read %s, pads=%d\n",argv[i],p-> pads);
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
    else if (sscanf(argv[i],"epsilon=%lg",&val)==1){ /** `epsilon` */
      p->epsilon = val;
      if (p->debug&1)
	printf("# read %s, epsilon=%g\n",argv[i],p-> epsilon);
      if((val<0) || (val>1))
	ERROR("arg[%d]=%s invalid probability cut-off, read %g",i,argv[i],val);
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
    else if (sscanf(argv[i],"maxW=%d",&dbg)==1){ /** `maxW` */
      p -> maxW = dbg;
      if (p->debug&1){
	printf("# read %s, maxW=%d\n",argv[i],p-> maxW);
	if (p->maxW <= 0)
	  printf("# no hard limit on error/codeword weight\n");
	else
	  printf("# hard upper limit w<=%d on codeword weight\n",p->maxW);
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
      switch(p->mode){
      case -1:  /** todo: insert help messages specific to each mode */
      default:
	printf( USAGE , argv[0],argv[0]);
	break;
      }    
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
  
  if(p->mode==-1)
    p->mode=0; /** the default value if `mode` is not set on command line */
  
  if (p->seed <= 0){
    long long int seed_old= - p->seed; 
    p->seed=-(seed_old) + time(NULL)+1000000ul*getpid(); /* ensure a different seed even if started at the same time */
    if((p->debug)&&(p->mode!=3))
      printf("# initializing seed=%lld from time(NULL)+1000000ul*getpid()+%lld\n",p->seed, seed_old);
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
    p->rankH = rank_csr(p->mH);
    p->ncws = p->nvar - p->rankH;
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
      /** WARNING: this does not necessarily have minimal row weights */
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
  else if (p->finP){/** read probabilities */
    int rows, cols, siz;
    p->vP = dbl_mm_read(p->finP, &rows, &cols, &siz, NULL);
    assert(rows == 1);
    assert(cols == p->nvar);
    if(p->debug&1)
      printf("# read %d probability values from %s\n", cols, p->finP);
  }
  if(!p->vP)
    ERROR("probabilities missing, specify 'fdem', 'finP', or 'useP'");
    
  LLR_table = init_LLR_tables(p->d1,p->d2,p->d3);

  p->dE = llr_from_dbl(p->dEdbl);
  
  switch(p->mode){
  case 1: /** both `mode=1` (BP) and `mode=0` */
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
    if(p->classical){
      p->mL=do_L_classical(p->mH, p);
      p->ncws = p->mL->rows;
    }
    
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

    if ((p->mode==0)&&(p->submode!=0)&&(p->submode!=1))
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
    
  case 3: /** read in DEM file and output the H, L, G, K matrices and P vector */
    if(strcmp(p->fout,"stdout")==0)
      p->use_stdout=1;
    if (p->submode>=64)
      ERROR(" mode=%d : submode='%d' unsupported\n",
	    p->mode, p->submode);
    if(p->submode & 32){ // matrix modification
      if(!(p->submode &31))
	p->submode ^= (4+8+16);
      else 
	ERROR("mode=%d : submode=%d should equal 32 ('mode=3.32' for matrix modification)\n",p->mode,p->submode);
    }
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

    if(p->mL){
      mzd_t *mL0 = mzd_from_csr(NULL,p->mL);
      printf("matrix mL0:\n");  mzd_print(mL0);
      mzd_free(mL0);
    }
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
  if(p->mL)
    p->mL =  csr_free(p->mL);
  if(p->mLt)
    p->mLt = csr_free(p->mLt);
  if(p->mG)
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
    il1=read_01(p->mHe,p->file_det, &p->line_det, p->fdet, p->pads, p->debug);
    /** TODO: enable external processing of observables */
    il2=read_01(p->mLe,p->file_obs, &p->line_obs, p->fobs, 0, p->debug);
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
    il1=read_01(p->mE,p->file_err, &p->line_err, p->ferr, 0, p->debug);
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
  init_Ht(p);

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
    synd_tot = 0;
    synd_fail = 0;
    for (long long int iround = 1; iround <= rounds; iround++) {
        if (p->debug & 1) {
            printf("# starting round %lld of %lld fail=%lld total=%lld\n", iround, rounds, synd_fail, synd_tot);
            fflush(stdout);
        }

        if (!(ierr_tot = do_err_vecs(p)))
            break; /** no more rounds */

        mzd_t *mE0 = NULL;
        mzd_t *mE0t = NULL;

        if (p->submode & 1) { /** submode 1 */
            // In submode1, we use Maximum Likelihood decoder
            if (!p->mG || !p->mK) {
              for (int i = 0; i < p->nchk; i++) {
                  p->vP[i] = P_from_llr(p->vLLR[i]);
              }
              p->rankH = rank_csr(p->mH);
              p->rankL = rank_csr(p->mL);
              p->rankG = p->nvar - p->rankH - p->rankL;
              if (p->debug & 1){
                    printf("rankG=%d", p->rankG);
              }
              p->mG = do_G_matrix(p->mHt, p->mLt, p->vLLR, p->rankG, p->debug);
              int rankG = rank_csr(p->mG);
              if (p->debug & 1){
                        printf("G matrix created with: rankG=%d \n", rankG);
              }
              p->mK = Lx_for_CSS_code(p->mG, p->mH);
              int rankK = rank_csr(p->mK);
              if (p->debug & 1){
                        printf("K matrix created with: rankK=%d \n", rankK);
              }
            }
            #ifndef NDEBUG  /** need `mHe` later */
            mzd_t *mS = mzd_copy(NULL, p->mHe);
            mE0 = do_decode_BAR(mS, p); /** each row a decoded error vector */
            mzd_free(mS);
            mS = NULL;
            #else
            mE0 = do_decode_BAR(p->mHe, p); /** each row a decoded error vector */
            #endif /* NDEBUG */
        } else {
            #ifndef NDEBUG  /** need `mHe` later */
            mzd_t *mS = mzd_copy(NULL, p->mHe);
            mE0 = do_decode(mS, p); /** each row a decoded error vector */
            mzd_free(mS);
            mS = NULL;
            #else
            mE0 = do_decode(p->mHe, p); /** each row a decoded error vector */
            #endif /* NDEBUG */
        }

        mE0t = mzd_transpose(NULL, mE0);
        mzd_free(mE0);
        mE0 = NULL;

        #ifndef NDEBUG
        mzd_t *prodHe = csr_mzd_mul(NULL, p->mH, mE0t, 1);
        mzd_add(prodHe, prodHe, p->mHe);
        if (!mzd_is_zero(prodHe)) {
            if ((p->debug & 512) || (p->nvec <= 64)) {
                printf("syndromes difference:\n");
                mzd_print(prodHe);
            }
            ERROR("some syndromes are not matched!\n");
        }
        mzd_free(prodHe);
        prodHe = NULL;
        #endif

        mzd_t *prodLe = csr_mzd_mul(NULL, p->mL, mE0t, 1);

        if (p->debug & 512) { /** print matrices */
            printf("prodLe:\n");
            mzd_print(prodLe);
            printf("mLe:\n");
            mzd_print(p->mLe);
        }

        mzd_add(prodLe, prodLe, p->mLe);

        int fails = 0;
        for (rci_t ic = 0; ic < ierr_tot; ic++) {
            rci_t ir = 0;
            if (mzd_find_pivot(prodLe, ir, ic, &ir, &ic)) {
                fails++;
            } else /** no more pivots */
                break;
        }
        /** update the global counts */
        synd_tot += ierr_tot; /** was: prodLe->ncols */
        synd_fail += fails;
        mzd_free(prodLe);
        prodLe = NULL;
        if ((p->nfail > 0) && (synd_fail >= p->nfail))
            break;
    }
    if (p->debug & 1)
        printf("# fail_fraction total_cnt succes_cnt\n");
    printf(" %g %lld %lld # %s\n", (double) synd_fail / synd_tot, synd_tot, synd_tot - synd_fail, p->fdem ? p->fdem : p->finH);

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
    do_LLR_dist(p, p->classical);
    do_hash_fail_prob(p); /** output results */
    if(p->outC){
      char * name;
      if(p->fdem)
	name=p->fdem;
      else if (p->finH)
	name=p->finH;
      else
	name="(unknown source)";
      size_t size = 1 + snprintf(NULL, 0, "codewords computed from '%s', maxW=%d ", name, p->maxW);
      char *buf = malloc(size * sizeof(char));
      if(!buf)
	ERROR("memory allocation");
      snprintf(buf, size, "# codewords computed from '%s', maxW=%d ", name, p->maxW);
      
      long long cnt=nzlist_write(p->outC, buf, p);
      if(p->debug & 1)
	printf("# wrote %lld computed codewords to file %s\n",cnt,p->outC);
    }
    do_hash_clear(p);
    break;
    
  case 3: /** read in DEM file and output the H, L, G matrices and P vector */
    if(p->submode&32){ /** matrix transformation */
      if(!p->finC)
	ERROR("mode=%d submode=%d (bit 5 set) must set 'finC' for matrix transformation",
	      p->mode,p->submode);
      if((p->classical)||(!(p->mL)))
	ERROR("mode=%d submode=32 (bit 5 set) must be a quantum code for matrix transformation",p->mode);
      p->num_cws = nzlist_read(p->finC, p);
      if(p->debug&1)
	printf("# %lld codewords read from %s ...",p->num_cws, p->finC);
      do_hash_verify_CW(p->mHt, NULL, p);
      if(p->debug&1)
	printf("all verified\n");
      do_hash_min(p); /** set values `minE` and `minW` from stored CWs */
      if(p->minW > 3)
	ERROR("minW=%d, codewords of weight <= 3 required for matrix transformation\n",p->minW);
      star_triangle(p->mHt, p->mLt, p->vLLR, p->codewords, p->debug);

      /** update the three matrices and the error vector */
      csr_free(p->mH);
      p->mH = csr_transpose(NULL, p->mHt);
      csr_free(p->mL);
      p->mL = csr_transpose(NULL,p->mLt);
      p->nchk = p->mH->rows;
      p->nvar = p->mH->cols;
      for(int i=0; i < p->nchk; i++)
	p->vP[i] = P_from_llr(p->vLLR[i]);
      
      assert((p->submode & (4+8+16)) == 28); /** only write Hx, Lx, and P */
    }
    size = snprintf(NULL, 0, "H matrix from DEM file %s", p->fdem);
    if(!(comment = malloc(size + 1)))
      ERROR("memory allocation");
    sprintf(comment, "H matrix from DEM file %s", p->fdem);
    
    if(!(p->submode & 31) || (p->submode & 4)){ /** see USAGE in `vecdec.h` for codes */
      if(p->debug&1)
	printf("# writing H=Hx matrix [ %d x %d ] to \t%s%s\n",
	       p->mH->rows, p->mH->cols, p->fout, p->use_stdout ? "\n" :"H.mmx");
      csr_mm_write(p->fout,"H.mmx",p->mH,comment);
    }
    if(!(p->submode & 31) || (p->submode & 8)){
      if(p->mL){
	if(p->debug&1)
	  printf("# writing L=Lx matrix [ %d x %d ] to \t%s%s\n",
		 p->mL->rows, p->mL->cols, p->fout, p->use_stdout ? "\n" :"L.mmx");
	comment[0]='L';
	csr_mm_write(p->fout,"L.mmx",p->mL,comment);
      }
      else {
	if(p->debug&1)
	  printf("# skippling L=Lx matrix for a classical code\n");	
      }
    }

    if(!(p->submode & 31) || (p->submode & 16)){      
      if(p->debug&1)
	printf("# writing P vector [ %d ] to      \t%s%s\n",
	       p->nvar, p->fout, p->use_stdout ? "\n" :"P.mmx");
      comment[0]='P';
      dbl_mm_write(p->fout,"P.mmx",1,p->nvar,p->vP,comment);
    }    

    if(!(p->submode & 31) || (p->submode & 1) || (p->submode & 2)){/** need `G=Hz` or `L=Lz` matrices */
      if(p->finC){/** use the list of codewords from file */
	nzlist_read(p->finC, p);
      }
      if(! p->rankH)
	p->rankH=rank_csr(p->mH);
      if (p->mL){
	p->rankL=rank_csr(p->mL);
	p->rankG=p->nvar - p->rankH - p->rankL;
	if(p->ncws != p->rankL)
	  ERROR("code parameters mismatch: rkH=%d rkL=%d Num(codewords)=%d\n",p->rankH, p->rankL, p->ncws);
	if(p->debug &1)
	  printf("quantum code parameters: n=%d k=%d rkH=%d rkG=%d\n",p->nvar,p->ncws,p->rankH, p->rankG);
      }
      else{
	p->rankL=p->nvar - p->rankH;
	p->rankG=0;
	if(p->debug &1)
	  printf("classical code parameters: n=%d k=%d rkH=%d \n",p->nvar,p->rankL,p->rankH);       
      }
      
      if(!(p->submode & 31) || (p->submode & 1)){
	if(p->classical)
	  ERROR("mode=%d submode=%d : create G matrix not supported for a classical code\n"
		"use 'mode=3.2' to create K matrix instead (dual to H)\n",p->mode,p->submode);
	if(!p->mG){
	  if(p->debug&1)
	    printf("# creating G=Hz matrix and writing to\t%s%s\n",
		   p->fout, p->use_stdout ? "\n" :"G.mmx");
	  if(p->finC)
	    p->mG = do_G_from_C(p->mLt,p->codewords,p->rankG,p->minW, p->minW + p->dW, p->debug);
	  else 
	    p->mG = do_G_matrix(p->mHt,p->mLt,p->vLLR, p->rankG, p->debug);
	  comment[0]='G';
	}
	else{
	  if(p->debug&1)
	    printf("# writing G=Hz matrix to\t%s%s\n",
		   p->fout, p->use_stdout ? "\n" :"G.mmx");
	  
	  size_t size2 = snprintf(NULL, 0, "G matrix from file %s", p->finG);
	  if(size2>size){
	    comment = realloc(comment, size2 + 1);
	    size=size2;
	  }    
	  assert(comment);
	  sprintf(comment, "G matrix from file %s", p->finG);
#ifndef NDEBUG	  
	  int rankG=rank_csr(p->mG);
	  if(rankG != p->rankG)
	    ERROR("Incorrect rank=%d of provided matrix G=Hz, expected %d",rankG, p->rankG);
#endif 	  	  
	}
	csr_mm_write(p->fout,"G.mmx",p->mG,comment);
      }

      if(!(p->submode & 31) || (p->submode & 2)){      
	if(!p->mK){ /** create `Lz` */
	  if(p->debug&1)
	    printf("# creating K=Lz matrix and writing to\t%s%s\n",
		   p->fout, p->use_stdout ? "\n" :"K.mmx");
	  if(p->finC){/** use the list of codewords from file */
	    //	    nzlist_read(p->finC, p);
	    p->mK = do_K_from_C(p->mLt, p->codewords, p->ncws, p->nvar, p->minW, p->minW+p->dW, p->debug);
	    //	  csr_out(p->mK);
	  }
	  else{
	    if ((p->mG)&&(p->mH)){
	      p->mK=Lx_for_CSS_code(p->mG,p->mH);
	      comment[0]='K';
	    }
	    else
	      ERROR("use mode=3.3 or specify the codewords file `finC=...` or generator matrix `finG=...`");
	  }
	}
	else{
	  size_t size2 = snprintf(NULL, 0, "K matrix from file %s", p->finK);
	  if(size2>size)
	    if(!(comment = realloc(comment, size2 + 1)))
	      ERROR("memory allocation");
	  sprintf(comment, "K matrix from file %s", p->finK);
#ifndef NDEBUG	  
	  int rankK=rank_csr(p->mK);
	  if(rankK != p->rankL)
	    ERROR("Incorrect rank=%d of provided matrix K=Lz, expected %d",rankK, p->rankL);
#endif 	  	  	
	}
	csr_mm_write(p->fout,"K.mmx",p->mK,comment);
      }
    }
    
    free(comment);
    break;
    
  default:
    ERROR("mode=%d not supported\n",p->mode);
    break;
  }

  var_kill(p);
  return 0;
}
