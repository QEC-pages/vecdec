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
#include <limits.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <m4ri/m4ri.h>
#include <m4ri/mzd.h>
#include "utils.h"
#include "util_m4ri.h"
#include "vecdec.h"
#include "qllr.h"
#include <stdbool.h>
#include <float.h>

params_t prm={ .nchk=-1, .nvar=-1, .ncws=-1, .steps=50, .pads=0,
  .rankH=0, .rankG=-1, .rankL=-1, 
  .lerr=-1, .maxosd=100, .swait=0, .maxC=0,
  .dW=0, .minW=INT_MAX, .maxW=0, .dE=-1, .dEdbl=-1, .minE=INT_MAX,
  .bpalpha=1, .bpbeta=1, .bpgamma=0.5, .bpretry=1, 
  .uW=2, .uR=4, //.uEdbl=-1, .uE=-1,
  .numU=0, .numE=0, .maxU=0,
  .hashU_error=NULL, .hashU_syndr=NULL, .permHe=NULL,
  .nvec=1024, .ntot=1, .nfail=0, .seed=0, .epsilon=1e-8, .useP=0, .mulP=0, .dmin=0,
  .debug=1, .fdem=NULL, .fout="tmp",
  .fdet=NULL, .fobs=NULL,  .ferr=NULL,
  .gdet=NULL, .gobs=NULL, // .gerr=NULL,
  .pdet=NULL, .pobs=NULL,  .perr=NULL,  
  .mode=-1, .submode=0, .use_stdout=0, 
  .LLRmin=0, .LLRmax=0, .codewords=NULL, .num_cws=0,
  .finH=NULL, .finL=NULL, .finG=NULL, .finK=NULL, .finP=NULL,
  .finC=NULL, .outC=NULL, 
  .finU=NULL, .outU=NULL, 
  .vP=NULL, .vLLR=NULL, .mH=NULL, .mHt=NULL,
  .mL=NULL, .mLt=NULL,  .mA=NULL, .mAt=NULL,
  .internal=0, 
  .file_err=NULL,  .file_det=NULL, .file_obs=NULL,
  //  .file_err_g=NULL,
  .file_gdet=NULL, .file_gobs=NULL,
  .file_perr=NULL,  .file_pdet=NULL, .file_pobs=NULL,
  .line_err=0,   .line_er0=0,  .line_det=0, .line_obs=0,
  .mE0=NULL,
  .mE=NULL, .mHe=NULL, .mLe=NULL, .mHeT=NULL, .mLeT=NULL,
  .nzH=0, .nzL=0,
  .buffer=NULL, .buffer_size = 0, .v0=NULL, .v1=NULL,
  .err=NULL, .svec=NULL, .obs=NULL,
  .ufl=NULL
};

params_t prm_default={  .steps=50, .pads=0, 
  .lerr=-1, .maxosd=100, .bpgamma=0.5, .bpretry=1, .swait=0, .maxC=0,
  .dW=0, .minW=INT_MAX, .maxW=0, .dE=-1, .dEdbl=-1, .minE=INT_MAX,
  .uW=2, .uR=4, //.uEdbl=-1, .uE=-1,
  .maxU=0, .bpalpha=1, .bpbeta=1,
  .nvec=1024, .ntot=1, .nfail=0, .seed=0, .epsilon=1e-8, .useP=0, .mulP=0, .dmin=0,
  .debug=1, .fout="tmp", .ferr=NULL,
  .mode=-1, .submode=0, .use_stdout=0, 
};

/** various success counters */
long long int cnt[EXTR_MAX];
long long int iter1[EXTR_MAX]; /** sums of BP iteration numbers */
long long int iter2[EXTR_MAX]; /** sums of BP iteration numbers squared */


/** @brief return `1` if `one` is a subvector of `two`.  
 *  @details input vectors  `one` and `two` should have  sorted coordinates */
static inline int is_subvec(const one_vec_t *const one, const one_vec_t *const two){
  assert(one->weight < two->weight); /* sanity check */
  for(int i=0, j=0; i < one->weight; i++){
    if(++j == two->weight)
      return 0; 
    while(one->arr[i] > two->arr[j]){
      if((++j) == two->weight)
	return 0; 
    }
    if(one->arr[i] < two->arr[j])
      return 0;
  }
  return 1;
}

/** @brief remove all entries from hash which contain smaller weight vector as subvectors
 * 
 */


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

/** @brief calculate exact fail probability for one vector */
static inline double do_prob_one_vec_exact(const one_vec_t * const pvec, const params_t * const p){
  _maybe_unused const int max_len = sizeof(unsigned long)*8 -1; 
  double ans=0;
  const int w = pvec->weight;
  assert(w<=max_len);
  const unsigned long min = 1;
  const unsigned long max = (1<<w)-1;
  /** TODO: accurately estimate limits using `min_llr` and `max_llar` */
  for(unsigned long bitmap=min; bitmap<max; bitmap++){
    double prod0=1, prod1=1;
    unsigned long bit = 1<<(w-1);
    const int *idx = pvec->arr;
    /** TODO: speed-up this loop to calculate faster.  Use array of
     * previously computed values and only start with the last bit changed. */
    for(int i=w-1; i>=0; i--, bit>>=1, idx++){
      const double prob = p->vP[*idx];
      if(bitmap&bit){
	prod1 *= prob;
	prod0 *= (1-prob);
      }
      else{ 
	prod0 *= prob;
	prod1 *= (1-prob);
      }
    }
    if (prod1 > 1.000001*prod0)
      ans += prod0;
    else if ((prod1 < 1.000001*prod0) && (prod1 > 0.999999*prod0))
      ans += 0.5*prod0; /** almost the same; added for numerical stability */
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
 * The resulting matrix should have `k` rows of weight `1` each.
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
  if(!f){
    if ((p->outC ==NULL) || (strcmp(fnam,p->outC)!=0)){      
      printf("codeword input file I/O ERROR: %s, outC=%s\n", strerror(errno),p->outC);
      ERROR("can't open file %s for reading",fnam);
    }
    else
      return 0; /** not an error: `finC = outC` */
  }
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
  if(p->debug&1)
    printf("read %lld codewords from %s, total %lld\n",count, fnam, p->num_cws);
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

/** @brief using codewords in the hash, estimate fail probability and min weigth */
int do_hash_fail_prob( params_t * const p){

  one_vec_t *pvec;
  
  if(p->debug & 1024) {/** `print` the list of cws found by energy */
    HASH_SORT(p->codewords, by_energy);
    for(pvec = p->codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next))
      if ((p->maxW==0) || ((p->maxW != 0) && (pvec->weight <= p->maxW)))
	print_one_vec(pvec);
  }
  /** finally calculate and output fail probability here */
  /** TODO: use the limit on `W` and `E` (just ignore codewords outside the limit) */
  double pfail=0, pfail_two=0, pmax=0, pmax_two=0;
  int count = 0, count_min = 0, count_tot = 0;
  int minW=p->nvar + 1;
  for(pvec = p->codewords; pvec != NULL; pvec=(one_vec_t *)(pvec->hh.next)){
    if ((p->maxW==0) || ((p->maxW != 0) && (pvec->weight <= p->maxW))){
      if(p->submode &1){ /* original simple estimate */
	double prob=do_prob_one_vec(pvec, p);
	pfail += prob;
	if(prob>pmax)
	  pmax=prob;
      }
      if(p->submode &2){ /* exact fail prob */      
	double prob_two=do_prob_one_vec_exact(pvec, p);
	pfail_two += prob_two;
	if(prob_two>pmax_two)
	  pmax_two=prob_two;
      }
      if(minW > pvec->weight){
	minW=pvec->weight;
	count_min=1;
      }
      else if (minW == pvec->weight)
	count_min++;
      count ++;
    }
    count_tot++;
  }  
  /** todo: prefactor calculation */
  if(p->debug&1)
    printf("# %s%smin_weight N_min N_use N_tot\n",p->submode&1 ?"sumP(fail) maxP(fail) ":"", p->submode&2 ?"sumP_exact(fail) maxP_exact(fail) ":"" );
  if (p->submode & 1)
    printf("%g %g ",pfail, pmax);
  if (p->submode & 2)
    printf("%g %g ",pfail_two, pmax_two);
  printf("%d %d %d %d\n", minW, count_min, count, count_tot);
  
  return minW; 
}

void do_hash_clear(params_t *const p){
    /** prescribed way to clean the hashing table */
  one_vec_t *cw, *tmp;
  HASH_ITER(hh, p->codewords, cw, tmp) {
    HASH_DEL(p->codewords, cw);
    free(cw);
  }
}

/** @brief remove reducible codewords from hash
 * @input min_dW the minimum weight of a `trivial` codeword (weight increment) */
void do_hash_remove_reduc(const int min_dW, params_t *const p){
  int max=0, siz = p->maxW > 0 ? p->maxW+1 : 30 ;
 
  one_vec_t *** by_w = calloc(siz, sizeof(one_vec_t **));
  if(!by_w) ERROR("memory allocation");
  one_vec_t *cw, *tmp;
  
  HASH_SORT(p->codewords, by_weight_pos);

  /** distribute `irreducible` vectors by weight and 1st entry position */
  int w;
  for(w = p->minW; p->codewords != NULL ; w++){
    if (w >= siz){
      siz=2*siz;
      by_w = realloc(by_w, sizeof(one_vec_t **) * siz);
      if(!by_w) ERROR("memory allocation");
    }
    const size_t keylen = w * sizeof(int);
    by_w[w]=calloc(p->nvar, sizeof(one_vec_t **));
    HASH_ITER(hh, p->codewords, cw, tmp) {
      if(cw->weight == w){
	int good=1;  /** verify irreducibility */
	for(int w1=p->minW; w1 + min_dW <= w; w1++){
	  const int b1max = cw->arr[w-1] - w + 2 > p->nvar ? cw->arr[w-1] - w + 2 : p->nvar - 1;
	  for(int beg1 = cw->arr[0]; beg1 < b1max; beg1++){
	    one_vec_t *cw1, *tmp1;
	    HASH_ITER(hh, by_w[w1][beg1], cw1, tmp1){
	      if(is_subvec(cw1, cw)){
		if(p->debug & 1){
		  printf("reducible pair w1=%d w=%d:\n",w1,w);
		  print_one_vec(cw1);
		  print_one_vec(cw);
		}
		good=0;
		break;	       
	      }
	    }
	    if(!good)
	      break;
	  }
	  if(!good)
	    break;
	}
	HASH_DEL(p->codewords,cw);
	if(good){
	  //	  	  printf("adding to w=%d b=%d cw: ",w,cw->arr[0]); print_one_vec(cw);
	  HASH_ADD(hh, by_w[w][cw->arr[0]], arr, keylen, cw);
	}
	else{
	  //	  printf("skipping cw: "); print_one_vec(cw);
	  free(cw);
	}
      }
      else
	break;
    }
  }
  max=w;
  if(p->codewords !=NULL)
    ERROR("unexpected");
  long long int count=0;
  /** and now move the entries back to `codewords` */
  for(w=p->minW; w < max; w++){
    const int keylen = w*sizeof(int);
    for(int beg = 0; beg < p->nvar; beg++){
      if(by_w[w][beg]){
	HASH_ITER(hh, by_w[w][beg], cw, tmp){
	  HASH_DEL(by_w[w][beg],cw);
	  HASH_ADD(hh, p->codewords, arr, keylen, cw);
	  count++;
	}
      }
    }
    free(by_w[w]);
  }
  free(by_w);
  if(p->debug & 1){
    if (p->num_cws > count)
      printf("# removed %lld reducible codewords count=%lld -> %lld\n",
	     p->num_cws - count, p->num_cws, count);
    else
      printf("# scanned for reducible codewords, found none of %lld\n",
	     p->num_cws);
  }
  p->num_cws = count;
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
    dst->nz = src->nz;  // Copy the nz value, which might be -1 for compressed matrices

    // For compressed matrices, we need to calculate the actual number of non-zero elements
    int actual_nz = (src->nz == -1) ? src->p[src->rows] : src->nz;

    dst->p = malloc((src->rows + 1) * sizeof(int));
    if (dst->p == NULL) {
        free(dst);
        ERROR("Memory allocation failed for CSR matrix data.\n");
        return NULL;
    }

    dst->i = malloc(dst->nzmax * sizeof(int));
    if (dst->i == NULL) {
        free(dst->p);
        free(dst);
        ERROR("Memory allocation failed for CSR matrix data.\n");
        return NULL;
    }

    memcpy(dst->p, src->p, (src->rows + 1) * sizeof(int));
    memcpy(dst->i, src->i, actual_nz * sizeof(int));

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
    int new_total_nonzeros = total_nonzeros - old_nz + nz_required;

    if (new_total_nonzeros > A->nzmax) {
        int new_nzmax = new_total_nonzeros * 2;  // Double the size for future growth
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

    for (int k = i + 1; k <= A->rows; k++) {
        A->p[k] += nz_required - old_nz;
    }

    free(toggle_columns);
}






/**
 * @brief Calculate the energy of row `i` in CSR matrix.
 * @param coeff Coefficients array for energy calculation.
 * @param csr The CSR matrix.
 * @param i The row index.
 * @return The energy of the row, or 0 if an error occurs.
 */
qllr_t csr_row_energ(qllr_t *coeff, const csr_t *csr, const int i) {
    // Initial parameter checks
    if (coeff == NULL || csr == NULL || csr->p == NULL || csr->i == NULL) {
        ERROR("Invalid input parameters in csr_row_energ: coeff=%p, csr=%p\n", coeff, csr);
        return 0;
    }

    if (i < 0 || i >= csr->rows) {
        ERROR("Invalid row index in csr_row_energ: i=%d, rows=%d\n", i, csr->rows);
        return 0;
    }

    // Check if the matrix is in compressed form
    if (csr->nz != -1) {
        ERROR("Matrix must be in compressed form for csr_row_energ.\n");
        return 0;
    }

    qllr_t ans = 0;
    int start = csr->p[i];
    int end = csr->p[i + 1];

    // Check row pointer range
    if (start < 0 || end > csr->nzmax || start > end) {
        ERROR("Invalid row pointer range in csr_row_energ: start=%d, end=%d, nzmax=%d\n", start, end, csr->nzmax);
        return 0;
    }

    for (int j = start; j < end; ++j) {
        // Check if j is within valid range
        if (j < 0 || j >= csr->nzmax) {
            ERROR("Index j (%d) out of bounds in csr_row_energ. nzmax: %d\n", j, csr->nzmax);
            return 0;
        }

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
 * @return The change in energy, or 0 if an error occurs.
 */
qllr_t csr_calculate_energy_change(qllr_t *coeff, const csr_t *A, const csr_t *B, const int i, const int row_to_add) {
    // Initial parameter checks
    if (coeff == NULL || A == NULL || B == NULL || A->p == NULL || A->i == NULL || B->p == NULL || B->i == NULL) {
        ERROR("Invalid input parameters in csr_calculate_energy_change: coeff=%p, A=%p, B=%p\n", coeff, A, B);
        return 0;
    }

    if (i < 0 || i >= A->rows || row_to_add < 0 || row_to_add >= B->rows) {
        ERROR("Invalid row indices in csr_calculate_energy_change: i=%d, A->rows=%d, row_to_add=%d, B->rows=%d\n",
              i, A->rows, row_to_add, B->rows);
        return 0;
    }

    // Check if matrices are in compressed form
    if (A->nz != -1 || B->nz != -1) {
        ERROR("Matrices must be in compressed form for csr_calculate_energy_change.\n");
        return 0;
    }

    qllr_t energy_change = 0;

    int start_A = A->p[i];
    int end_A = A->p[i + 1];
    int start_B = B->p[row_to_add];
    int end_B = B->p[row_to_add + 1];

    // Check row pointer ranges
    if (start_A < 0 || end_A > A->nzmax || start_A > end_A ||
        start_B < 0 || end_B > B->nzmax || start_B > end_B) {
        ERROR("Invalid row pointer ranges: A(start=%d, end=%d, nzmax=%d), B(start=%d, end=%d, nzmax=%d)\n",
              start_A, end_A, A->nzmax, start_B, end_B, B->nzmax);
        return 0;
    }

    int *toggle_columns = calloc(A->cols, sizeof(int));
    if (toggle_columns == NULL) {
        ERROR("Memory allocation failed in csr_calculate_energy_change.\n");
        return 0;
    }

    // Process matrix A
    for (int j = start_A; j < end_A; ++j) {
        if (j < 0 || j >= A->nzmax) {
            ERROR("Index j (%d) out of bounds for matrix A. nzmax: %d\n", j, A->nzmax);
            free(toggle_columns);
            return 0;
        }
        int col = A->i[j];
        if (col < 0 || col >= A->cols) {
            ERROR("Invalid column index in matrix A: col=%d, cols=%d\n", col, A->cols);
            free(toggle_columns);
            return 0;
        }
        toggle_columns[col] ^= 1;
    }

    // Process matrix B
    for (int j = start_B; j < end_B; ++j) {
        if (j < 0 || j >= B->nzmax) {
            ERROR("Index j (%d) out of bounds for matrix B. nzmax: %d\n", j, B->nzmax);
            free(toggle_columns);
            return 0;
        }
        int col = B->i[j];
        if (col < 0 || col >= B->cols) {
            ERROR("Invalid column index in matrix B: col=%d, cols=%d\n", col, B->cols);
            free(toggle_columns);
            return 0;
        }
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




/* Decompress a CSR matrix that was previously compressed */
void csr_decompress(csr_t *mat) {
    if (mat->nz != -1) {
        fprintf(stderr, "Error: Matrix is not in a compressed state.\n");
        return;
    }

    int nz = 0;
    for (int r = 0; r < mat->rows; r++) {
        nz += mat->p[r + 1] - mat->p[r];
    }

    mat->nz = nz;
}





/**
 * @brief Replace the row i of destination CSR matrix with the row i of source CSR matrix.
 * @param dest The destination CSR matrix.
 * @param src The source CSR matrix.
 * @param i The row index to replace.
 */
void csr_replace_row(csr_t *dst, const csr_t *src, int row) {
    // Check if the matrix is compressed, and decompress if necessary
    if (dst->nz == -1) {
        csr_decompress(dst);
    }
    if (src->nz == -1) {
        csr_decompress((csr_t *)src); // src is const, so cast it temporarily
    }

    // Validate input matrices and row index
    if (!dst || !src) {
        fprintf(stderr, "Error: Null pointer passed to csr_replace_row.\n");
        return;
    }

    if (row < 0 || row >= dst->rows || row >= src->rows) {
        fprintf(stderr, "Error: Row index %d out of bounds. dst rows: %d, src rows: %d\n", row, dst->rows, src->rows);
        return;
    }

    // Determine the start and end pointers for the row in both matrices
    int dst_row_start = dst->p[row];
    int dst_row_end = dst->p[row + 1];
    int dst_row_nnz = dst_row_end - dst_row_start;

    int src_row_start = src->p[row];
    int src_row_end = src->p[row + 1];
    int src_row_nnz = src_row_end - src_row_start;

    // Ensure the source row is non-empty
    if (src_row_nnz == 0) {
        fprintf(stderr, "Warning: Source row %d is empty, no changes made.\n", row);
        return;
    }

    int nnz_diff = src_row_nnz - dst_row_nnz;

    // Calculate the new total non-zero entries required
    int new_nz = dst->nz + nnz_diff;

    // Check if reallocation is needed
    if (new_nz > dst->nzmax) {
        int new_nzmax = new_nz + 1024; // Add extra space to reduce frequent reallocations

        // Reallocate column indices array
        int *new_i = realloc(dst->i, new_nzmax * sizeof(int));
        if (!new_i) {
            fprintf(stderr, "Error: Failed to reallocate memory for column indices.\n");
            return;
        }
        dst->i = new_i;

        dst->nzmax = new_nzmax;
    }

    // Move data only if necessary
    if (nnz_diff != 0) {
        if (dst_row_end + nnz_diff < dst->nz && dst_row_end + nnz_diff >= 0) {
            memmove(dst->i + dst_row_end + nnz_diff, dst->i + dst_row_end, (dst->nz - dst_row_end) * sizeof(int));
        } else {
            //fprintf(stderr, "Error: Invalid memory move operation detected.\n");
            return;
        }

        for (int r = row + 1; r <= dst->rows; r++) {
            dst->p[r] += nnz_diff;
        }

        dst->nz = new_nz;
    }

    // Copy column indices from source to destination
    memcpy(dst->i + dst->p[row], src->i + src->p[row], src_row_nnz * sizeof(int));

    // printf("Replaced row %d: src nnz = %d, dst nnz = %d, nnz_diff = %d, new dst nz = %d\n",
    //        row, src_row_nnz, dst_row_nnz, nnz_diff, dst->nz);
}







/**
 * @brief Executes simulated annealing to adjust error matrix based on energy changes.
 * @param mEt0_csr The initial error matrix in CSR format to be updated.
 * @param p Structure containing model parameters and error probabilities.
 */
void simulate_annealing_mcmc(csr_t *mEt0_csr, params_t const * const p,int i) {
    int max_iterations = 400000;
    double final_temperature = 1.0;
    double cooling_rate = 0.6;
    double reheating_rate = 1.5;
    int check_interval = 5000;
    double min_acceptance_rate = 0.2;

    double temperature = 1e5;
    //printf("initial temp is %f\n",temperature);
  
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

// Helper function for Fermi function
double fermi_function(double x) {
    return 1.0 / (1.0 + exp(x));
}

double BAR_equation(double ***U, int **N, int i, int j, double Q_ratio) {
    // Check input validity
    if (!U || !N || i < 0 || j < 0 || !isfinite(Q_ratio)) {
        fprintf(stderr, "Error: Invalid input in BAR_equation\n");
        return NAN;
    }

    // Check if N[i][j] and N[j][i] are valid
    if (N[i][j] <= 0 || N[j][i] <= 0) {
        fprintf(stderr, "Error: Invalid sample count N[%d][%d] = %d, N[%d][%d] = %d\n", 
                i, j, N[i][j], j, i, N[j][i]);
        return NAN;
    }

    double C = log(Q_ratio * N[j][i] / N[i][j]);
    if (!isfinite(C)) {
        fprintf(stderr, "Error: Invalid C value in BAR_equation\n");
        return NAN;
    }

    // Check if U[j][i] and U[i][j] are valid pointers
    if (!U[j][i] || !U[i][j]) {
        fprintf(stderr, "Error: Null pointer encountered for U[%d][%d] or U[%d][%d]\n", j, i, i, j);
        return NAN;
    }

    double sum_i = 0.0, sum_j = 0.0;

    // Calculate sum for state i
    for (int n = 0; n < N[i][j]; n++) {
        double arg = U[j][i][n] - U[i][j][n] - C;
        if (!isfinite(arg)) {
            fprintf(stderr, "Warning: Invalid argument for fermi_function: %f\n", arg);
            continue;
        }
        sum_i += fermi_function(arg);
    }
    
    // Calculate sum for state j
    for (int n = 0; n < N[j][i]; n++) {
        double arg = U[i][j][n] - U[j][i][n] + C;
        if (!isfinite(arg)) {
            fprintf(stderr, "Warning: Invalid argument for fermi_function: %f\n", arg);
            continue;
        }
        sum_j += fermi_function(arg);
    }
    
    // Check if sum_i and sum_j are valid
    if (sum_i == 0.0 || sum_j == 0.0 || !isfinite(sum_i) || !isfinite(sum_j)) {
        //fprintf(stderr, "Error: Invalid sum_i (%f) or sum_j (%f)\n", sum_i, sum_j);
        return NAN;
    }

    // Calculate the final result
    double result = Q_ratio - (sum_j / N[j][i]) / (sum_i / N[i][j]) * exp(C);
    
    // Check if the result is valid
    if (!isfinite(result)) {
        fprintf(stderr, "Error: Invalid result in BAR_equation: %f\n", result);
        return NAN;
    }

    return result;
}



void solve_pairwise_BAR_ratio(double ***U, int **N, int K, double **Q_ratios) {
    long double tolerance = 1e-12;
    int max_iterations = 10000;

    for (int i = 0; i < K; i++) {
        for (int j = i + 1; j < K; j++) {
            long double Q_low = 0.1;  // Lower bound for Q_ratio
            long double Q_high = 10;  // Upper bound for Q_ratio
            long double Q_mid;

            for (int iter = 0; iter < max_iterations; iter++) {
                Q_mid = (Q_low + Q_high) / 2;  // Arithmetic mean
                double f_mid = BAR_equation(U, N, i, j, (double)Q_mid);

                if (fabs((double)f_mid) < tolerance) {
                    break;
                }

                if (f_mid > 0) {
                    Q_high = Q_mid;
                } else {
                    Q_low = Q_mid;
                }

                if (fabsl(Q_high - Q_low) < tolerance * fabsl(Q_mid)) {
                    break;
                }
            }

            Q_ratios[i][j] = (double)Q_mid;
            Q_ratios[j][i] = 1.0 / (double)Q_mid;
            //printf("Q_ratio(%d,%d) = %Le\n", i, j, Q_mid);
        }
    }
}

#define MAX_K 32
#define MAX_ITERATIONS 100000

typedef struct {
    double ***data;
    int **sizes;
} EnergyDifferences;

typedef struct {
    double **data;
} QRatios;

typedef struct {
    csr_t **vectors;
} ErrorVectors;

// Helper function to find the state with the largest partition function
int find_max_partition_function(QRatios *Q_ratios, int K) {
    int max_Q_index = 0;
    double max_Q = -INFINITY;
    for (int k = 0; k < K; k++) {
        double Q_k = 0.0;
        for (int j = 0; j < K; j++) {
            if (j != k) {
                Q_k += Q_ratios->data[k][j];
            }
        }
        if (Q_k > max_Q) {
            max_Q = Q_k;
            max_Q_index = k;
        }
    }
    return max_Q_index;
}

void BAR_MCMC(csr_t *mEt0_csr, csr_t *mEt_csr, params_t const * const p, int i) {
    int mK_rows = p->mK->rows;
    int K = 1 << mK_rows;  // 2^mK_rows

    if (K > MAX_K) {
        fprintf(stderr, "Error: K exceeds MAX_K in BAR_MCMC\n");
        return;
    }

    EnergyDifferences U = {0};
    QRatios Q_ratios = {0};
    ErrorVectors error_vectors = {0};

    // Allocate memory
    U.data = malloc(K * sizeof(double **));
    U.sizes = malloc(K * sizeof(int *));
    Q_ratios.data = malloc(K * sizeof(double *));
    error_vectors.vectors = malloc(K * sizeof(csr_t *));

    for (int k = 0; k < K; k++) {
        U.data[k] = malloc(K * sizeof(double *));
        U.sizes[k] = calloc(K, sizeof(int));
        Q_ratios.data[k] = malloc(K * sizeof(double));
        for (int j = 0; j < K; j++) {
            U.data[k][j] = NULL;
            Q_ratios.data[k][j] = 1.0;
        }
    }

    // Initialize error_vectors
    for (int k = 0; k < K; k++) {
        error_vectors.vectors[k] = csr_copy(mEt0_csr);
        if (!error_vectors.vectors[k]) {
            fprintf(stderr, "Failed to copy error vector for k=%d\n", k);
            goto cleanup;
        }
        // Apply mK transformations
        for (int row = 0; row < mK_rows; row++) {
            if (k & (1 << row)) {
                csr_add_row_to_row(error_vectors.vectors[k], i, p->mK, row);
            }
        }
    }

    // Monte Carlo sampling
    for (int k = 0; k < K; k++) {
        for (int j = k + 1; j < K; j++) {
            csr_t *state_A = csr_copy(error_vectors.vectors[k]);
            csr_t *state_B = csr_copy(error_vectors.vectors[j]);
            if (!state_A || !state_B) {
                fprintf(stderr, "Failed to copy states for k=%d, j=%d\n", k, j);
                csr_free(state_A);
                csr_free(state_B);
                goto cleanup;
            }

            // Sample state A and B
            sample_state(state_A, k, j, &U, &error_vectors, p, i);
            sample_state(state_B, j, k, &U, &error_vectors, p, i);

            csr_free(state_A);
            csr_free(state_B);
        }
    }

    // Solve BAR equations
    solve_pairwise_BAR_ratio(U.data, U.sizes, K, Q_ratios.data);

    // Find the state with the largest partition function
    int max_Q_index = find_max_partition_function(&Q_ratios, K);

    // Update mEt_csr with the chosen ensemble
    csr_replace_row(mEt_csr, error_vectors.vectors[max_Q_index], i);

cleanup:
    // Clean up
    for (int k = 0; k < K; k++) {
        for (int j = 0; j < K; j++) {
            free(U.data[k][j]);
        }
        free(U.data[k]);
        free(U.sizes[k]);
        free(Q_ratios.data[k]);
        csr_free(error_vectors.vectors[k]);
    }
    free(U.data);
    free(U.sizes);
    free(Q_ratios.data);
    free(error_vectors.vectors);
}

// Helper function for sampling states
void sample_state(csr_t *state, int k, int j, EnergyDifferences *U, ErrorVectors *error_vectors, params_t const * const p, int i) {
    int accepted = 0;
    int check_interval = 1000;
    double min_acceptance_rate = 0.004;
    double temperature = 1.0 * pow(10, 4);

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        unsigned int seed = (unsigned int)time(NULL) + iter;
        int row_to_add = rand_r(&seed) % (p->mG->rows);
        
        double delta_energy = csr_calculate_energy_change(p->vLLR, state, p->mG, i, row_to_add);
        
        if (delta_energy < 0 || exp(-delta_energy / temperature) > (double)rand_r(&seed) / RAND_MAX) {
            csr_add_row_to_row(state, i, p->mG, row_to_add);
            accepted++;
        }
        
        double U_state = csr_row_energ(p->vLLR, state, i);
        double U_other = csr_row_energ(p->vLLR, error_vectors->vectors[j], i);
        double energy_diff = U_other - U_state;

        // Increase array length and add new energy difference
        U->sizes[k][j]++;
        U->data[k][j] = realloc(U->data[k][j], U->sizes[k][j] * sizeof(double));
        if (!U->data[k][j]) {
            fprintf(stderr, "Memory reallocation failed for U->data[%d][%d]\n", k, j);
            return;
        }
        U->data[k][j][U->sizes[k][j] - 1] = energy_diff;

        if (iter % check_interval == check_interval - 1) {
            double acceptance_rate = (double)accepted / check_interval;
            if (acceptance_rate < min_acceptance_rate) {
                break;
            }
            accepted = 0;
        }
    }
}

/**
 * @brief Simulated Annealing-based syndrome decoding routine with improved stopping mechanism. (CSR form)
 * @param mS Matrix with syndromes (each column).
 * @param p Structure with error model information.
 * @return Binary matrix of min weight errors for each syndrome.
 */
mzd_t *do_decode_BAR(mzd_t *mS, params_t const * const p) {
    mzd_t *mHx_dense = mzd_from_csr(NULL, p->mH);
    if (!mHx_dense) {
        fprintf(stderr, "Failed to create dense matrix from CSR\n");
        return NULL;
    }

    mzd_t *mE_dense = mzd_init(mHx_dense->ncols, mS->ncols);
    if (!mE_dense) {
        fprintf(stderr, "Failed to initialize mE_dense\n");
        mzd_free(mHx_dense);
        return NULL;
    }

    mzp_t *perm = mzp_init(p->nvar);
    mzp_t *pivs = mzp_init(p->nvar);
    if (!perm || !pivs) {
        fprintf(stderr, "Failed to initialize permutations\n");
        mzd_free(mHx_dense);
        mzd_free(mE_dense);
        mzp_free(perm);
        mzp_free(pivs);
        return NULL;
    }

    perm = sort_by_llr(perm, p->vLLR, p);
    mzd_set_ui(mE_dense, 0);

    int rank = 0;
    if (p->steps > 0) {
        // Gaussian elimination
        for (int i = 0; i < p->nvar; i++) {
            int col = perm->values[i];
            int ret = twomat_gauss_one(mHx_dense, mS, col, rank);
            if (ret)
                pivs->values[rank++] = col;
        }

        // Calculate error vector and energy for each syndrome
        for (int i = 0; i < rank; i++)
            mzd_copy_row(mE_dense, pivs->values[i], mS, i);
    }

    mzp_free(perm);
    mzp_free(pivs);

    csr_t *mE0_csr = csr_from_mzd(NULL, mE_dense);
    if (!mE0_csr) {
        fprintf(stderr, "Failed to create CSR matrix from dense matrix\n");
        mzd_free(mHx_dense);
        mzd_free(mE_dense);
        return NULL;
    }

    csr_t *mEt0_csr = csr_init(NULL, mS->ncols, mHx_dense->ncols, mS->ncols * mHx_dense->ncols / 4);
    if (!mEt0_csr) {
        fprintf(stderr, "Failed to initialize mEt0_csr\n");
        mzd_free(mHx_dense);
        mzd_free(mE_dense);
        csr_free(mE0_csr);
        return NULL;
    }
    csr_transpose(mEt0_csr, mE0_csr);

    mzd_t *mEt_dense = mzd_init(mS->ncols, mHx_dense->ncols);
    if (!mEt_dense) {
        fprintf(stderr, "Failed to initialize mEt_dense\n");
        mzd_free(mHx_dense);
        mzd_free(mE_dense);
        csr_free(mE0_csr);
        csr_free(mEt0_csr);
        return NULL;
    }
    mzd_set_ui(mEt_dense, 0);

    csr_t *mEt_csr = csr_from_mzd(NULL, mEt_dense);
    if (!mEt_csr) {
        fprintf(stderr, "Failed to create CSR matrix from mEt_dense\n");
        mzd_free(mHx_dense);
        mzd_free(mE_dense);
        mzd_free(mEt_dense);
        csr_free(mE0_csr);
        csr_free(mEt0_csr);
        return NULL;
    }
    mEt_csr->nzmax = mS->ncols * mHx_dense->ncols / 4;

    // BAR_MCMC and simulated annealing
    for (int i = 0; i < mEt0_csr->rows; i++) {
        BAR_MCMC(mEt0_csr, mEt_csr, p, i);
    }

    csr_compress(mEt_csr);
    
    for (int i = 0; i < mEt0_csr->rows; i++) {
        simulate_annealing_mcmc(mEt_csr, p, i);
    }

    mzd_t *mEt = mzd_from_csr(NULL, mEt_csr);
    if (!mEt) {
        fprintf(stderr, "Failed to create dense matrix from final CSR matrix\n");
        mzd_free(mHx_dense);
        mzd_free(mE_dense);
        mzd_free(mEt_dense);
        csr_free(mE0_csr);
        csr_free(mEt0_csr);
        csr_free(mEt_csr);
        return NULL;
    }

    return mEt;
}

/** @brief one more matrix initialization routine 
 */
void init_Ht(params_t *p){
  const int n = p->nvar;
  p->mHt = csr_transpose(p->mHt, p->mH);
  //! construct v-v graph 
  //  csr_t *vv_gr = do_vv_graph(p->mH, p->mHt, p);
  if(p->uW>0){
    if(p->debug&1){
      if(p->uR==0)
	printf("# generating errors of weight up to %d for syndrome hash\n",p->uW);
      else 
	printf("# generating error clusters of w <= %d and radius <= %d for syndrome hash\n",p->uW, p->uR);
    }
    p->v0 = vec_init(p->nchk);
    p->v1 = vec_init(p->nchk);
    if(p->uW >=0)
      do_clusters(p);
  }
    if(p->debug&1)
  printf("# initialize matrices\n");
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

  /** scan for `debug` and `mode` first */
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
        if(p->debug &4)
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
      if(p->debug&4)
	printf("# read %s, mode=%d submode=%d [octal=%o]\n",argv[i],
	       p->mode,p->submode,p->submode);      
    }
  }
  
  for(int i=1; i<argc; i++){  
    if (0==strncmp(argv[i],"debug=",6)){ /** `debug` */
      /** do nothing */
    }
    else if (0==strncmp(argv[i],"mode=",5)){ /** `debug` */
      /** do nothing */
    }
    else if (sscanf(argv[i],"qllr1=%d",&dbg)==1){ /** QLLR `d1` param */
      p -> d1 = dbg;
      if (p->debug&4)
	printf("# read %s, QLLR parameter d1=%d\n",argv[i],p-> d1);
#ifndef USE_QLLR
      ERROR("program was compiled without \"USE_QLLR\" define");
#endif       
    }
    else if (sscanf(argv[i],"qllr2=%d",&dbg)==1){ /** QLLR `d2` param */
      p -> d2 = dbg;
      if (p->debug&4)
	printf("# read %s, QLLR parameter d2=%d\n",argv[i],p-> d2);
#ifndef USE_QLLR
      ERROR("program was compiled without \"USE_QLLR\" define");
#endif       
    }
    else if (sscanf(argv[i],"qllr3=%d",&dbg)==1){ /** QLLD `d3` param */
      p -> d3 = dbg;
      if (p->debug&4)
	printf("# read %s, QLLR parameter d3=%d\n",argv[i],p-> d3);
#ifndef USE_QLLR
      ERROR("program was compiled without \"USE_QLLR\" define");
#endif       
    }
    else if (sscanf(argv[i],"nvec=%d",&dbg)==1){ /** `nvec` */
      p -> nvec = dbg;
      if (p->debug&4)
	printf("# read %s, nvec=%d\n",argv[i],p-> nvec);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"pads=%d",&dbg)==1){ /** `pads` */
      p -> pads = dbg;
      if (p->debug&4)
	printf("# read %s, pads=%d\n",argv[i],p-> pads);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"dmin=%d",&dbg)==1){ /** `dmin` */
      p -> dmin = dbg;
      if (dbg<0)
	ERROR("command line argument %d dmin=%d must be non-negative\n",i,dbg);
      if (p->debug&4)
	printf("# read %s, dmin=%d\n",argv[i],p-> dmin);
      if(p->mode != 2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"steps=%d",&dbg)==1){ /** `steps` */
      p -> steps = dbg;
      if (p->debug&4)
	printf("# read %s, steps=%d\n",argv[i],p-> steps);
      if(p->mode>2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"swait=%d",&dbg)==1){ /** `swait` */
      p -> swait = dbg;
      if (p->debug&4)
	printf("# read %s, swait=%d\n",argv[i],p-> swait);
      if((p->mode >0)|| (p->mode!=2))
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"lerr=%d",&dbg)==1){ /** `lerr` */
      p -> lerr = dbg;
      if (p->debug&4)
	printf("# read %s, lerr=%d\n",argv[i],p-> lerr);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"maxosd=%d",&dbg)==1){ /** `maxosd` */
      p -> maxosd = dbg;
      if (p->debug&4)
	printf("# read %s, maxosd=%d\n",argv[i],p-> maxosd);
      if(p->mode!=1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"ntot=%lld",&lldbg)==1){ /** `ntot` */
      p -> ntot = lldbg;
      if (p->debug&4)
	printf("# read %s, ntot=%lld\n",argv[i],p-> ntot);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"maxC=%lld",&lldbg)==1){ /** `maxC` */
      p -> maxC = lldbg;
      if (p->debug&4)
	printf("# read %s, maxC=%lld\n",argv[i],p-> maxC);
      if(p->mode<2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"maxU=%lld",&lldbg)==1){ /** `maxU` */
      p -> maxU = lldbg;
      if (p->debug&4)
	printf("# read %s, maxU=%lld\n",argv[i],p-> maxU);
      if(p->mode>=2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
#if 0    
    else if (sscanf(argv[i],"uE=%lg",&val)==1){ /** `uE`  */
      p -> uEdbl = val;
      if (p->debug&4){
	printf("# read %s, uE=%g\n", argv[i], p->uEdbl);
	if (p->uEdbl < 0)
	  printf("# no limit on error/codeword energies to store\n");
	else
	  printf("# setting upper limit on error/codeword energies to store\n");
      }
      if(p->mode>=2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
#endif     
    else if (sscanf(argv[i],"uW=%d",&dbg)==1){ /** `uW` */
      p -> uW = dbg;
      if (p->debug&4){
	printf("# read %s, uW=%d\n",argv[i],p-> uW);
	if (p->uW < 0)
	  printf("# will not attempt pre-decoding\n");
	else if (p->uW == 0)
	  printf("# will only skip zero syndrome vectors\n");
	else
	  printf("# will pre-compute syndromes for clusters of weight up to %d\n",p->uW);
      }	
      if(p->mode>=2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"uR=%d",&dbg)==1){ /** `uR` */
      p -> uR = dbg;
      if (p->debug&4){
	printf("# read %s, uR=%d\n",argv[i],p-> uR);
	if (p->uR < 0)
	  ERROR("invalid uR value");
	else if (p->uR == 0)
	  printf("# will use any errors of weight up to uW=%d for syndrome hash\n",p->uW);
	else
	  printf("# will use r-local v-v error clusters for syndrome hash, r<=uR\n");
      }	
      if(p->mode>=2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"nfail=%lld",&lldbg)==1){ /** `nfail` */
      p -> nfail = lldbg;
      if (p->debug&4)
	printf("# read %s, nfail=%lld\n",argv[i],p-> nfail);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"seed=%lld",&lldbg)==1){ /** `seed` */
      p->seed=lldbg;
      if (p->debug&4)
	printf("# read %s, seed=%lld\n",argv[i],p->seed);
    }
    else if (sscanf(argv[i],"bpgamma=%lg",&val)==1){ /** `bpgamma` */
      p -> bpgamma = val; /** WARNING: no value verification!!! */
      if (p->debug&4)
	printf("# read %s, bpgamma=%g\n",argv[i],p-> bpgamma);
      if(p->mode!=1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"bpretry=%d",&dbg)==1){ /** `bpretry` */
      p -> bpretry = dbg;
      if(dbg<1)
	ERROR("read arg[%d]=%s, bpretry=%d should be at least 1",i,argv[i],dbg);
      if (p->debug&4)
	printf("# read %s, bpretry=%d\n",argv[i],p-> bpretry);
      if(p->mode!=1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"epsilon=%lg",&val)==1){ /** `epsilon` */
      p->epsilon = val;
      if (p->debug&4)
	printf("# read %s, epsilon=%g\n",argv[i],p-> epsilon);
      if((val<0) || (val>1))
	ERROR("arg[%d]=%s invalid probability cut-off, read %g",i,argv[i],val);
      if(p->mode>=0)
	ERROR("mode=%d, this parameter %s is irrelevant (for any mode)\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"useP=%lg",&val)==1){ /** `useP` */
      p->useP = val;
      if (p->debug&4)
	printf("# read %s, useP=%g\n",argv[i],p-> useP);
    }
    else if (sscanf(argv[i],"mulP=%lg",&val)==1){ /** `mulP` */
      p->mulP = val;
      if (p->debug&4)
	printf("# read %s, mulP=%g\n",argv[i],p-> mulP);
    }
    else if (sscanf(argv[i],"dE=%lg",&val)==1){ /** `dE` */
      p -> dEdbl = val;
      if (p->debug&4){
	printf("# read %s, dE=%g\n", argv[i], p->dEdbl);
	if (p->dEdbl < 0)
	  printf("# no limit on error/codeword energies to store\n");
	else
	  printf("# setting upper limit on error/codeword energies to store\n");
      }
      if(p->mode!=2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"dW=%d",&dbg)==1){ /** `dW` */
      p -> dW = dbg;
      if (p->debug&4){
	printf("# read %s, dW=%d\n",argv[i],p-> dW);
	if (p->dW < 0)
	  printf("# no limit on error/codeword weight to store\n");
	else
	  printf("# setting upper limit 'minW+%d' on error/codeword weight to store\n",p->dW);
      }	
      if((p->mode!=2) && (p->mode!=3))
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (sscanf(argv[i],"maxW=%d",&dbg)==1){ /** `maxW` */
      p -> maxW = dbg;
      if (p->debug&4){
	printf("# read %s, maxW=%d\n",argv[i],p-> maxW);
	if (p->maxW <= 0)
	  printf("# no hard limit on error/codeword weight\n");
	else
	  printf("# hard upper limit w<=%d on codeword weight\n",p->maxW);
      }	
      if(p->mode<2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"fout=",5)){ /** `fout` */
      if(strlen(argv[i])>5){
        p->fout = argv[i]+5;
	if (p->debug&4)
	  printf("# read %s, fout=%s\n",argv[i],p->fout);
      }
      else
	ERROR("Please specify argument for 'fout=[string]' w/o space\n");
      if(p->mode!=3)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"fdem=",5)){ /** fdem */
      if(strlen(argv[i])>5)
        p->fdem = argv[i]+5;
      else
        p->fdem = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, fdem=%s\n",argv[i],p->fdem);
    }
    else if (0==strncmp(argv[i],"finH=",5)){ /** `finH` */
      if(strlen(argv[i])>5)
        p->finH = argv[i]+5;
      else
        p->finH = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finH=%s\n",argv[i],p->finH);
    }
    else if (0==strncmp(argv[i],"finA=",5)){ /** `finA` */
      if(strlen(argv[i])>5)
        p->finA = argv[i]+5;
      else
        p->finA = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finA=%s\n",argv[i],p->finA);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
      p->mA=csr_mm_read(p->finA, p->mA, 0, p->debug);
      if(p->mAt != NULL)
	ERROR("must specify finA=%s only once\n",p->finA);
      p->mAt=csr_transpose(NULL,p->mA);
    }
    else if (0==strncmp(argv[i],"finP=",5)){ /** `finP` probabilities */
      if(strlen(argv[i])>5)
        p->finP = argv[i]+5;
      else
        p->finP = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finP=%s\n",argv[i],p->finP);
    }
    else if (0==strncmp(argv[i],"finL=",5)){ /** `finL` logical */
      if(strlen(argv[i])>5)
        p->finL = argv[i]+5;
      else
        p->finL = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finL=%s\n",argv[i],p->finL);
    }
    else if (0==strncmp(argv[i],"finK=",5)){ /** `finK` logical */
      if(strlen(argv[i])>5)
        p->finK = argv[i]+5;
      else
        p->finK = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finK=%s\n",argv[i],p->finK);
    }
    else if (0==strncmp(argv[i],"finC=",5)){ /** `finC` in low-weight codeword list */
      if(strlen(argv[i])>5)
        p->finC = argv[i]+5;
      else
        p->finC = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finC=%s\n",argv[i],p->finC);
    }
    else if (0==strncmp(argv[i],"outC=",5)){ /** `outC` out low-weight codeword list */
      if(strlen(argv[i])>5)
        p->outC = argv[i]+5;
      else
        p->outC = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, outC=%s\n",argv[i],p->outC);
    }
    else if (0==strncmp(argv[i],"finU=",5)){ /** `finU` in low-weight errors list */
      if(strlen(argv[i])>5)
        p->finU = argv[i]+5;
      else
        p->finU = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finU=%s\n",argv[i],p->finU);
      if(p->mode>=2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);      
    }
    else if (0==strncmp(argv[i],"outU=",5)){ /** `outU` out low-weight errors list */
      if(strlen(argv[i])>5)
        p->outU = argv[i]+5;
      else
        p->outU = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, outU=%s\n",argv[i],p->outU);
      if(p->mode>=2)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"finG=",5)){/** `finG` degeneracy generator matrix */
      if(strlen(argv[i])>5)
        p->finG = argv[i]+5;
      else
        p->finG = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, finG=%s\n",argv[i],p->finG);
    }
    else if (0==strncmp(argv[i],"fdet=",5)){ /** detector events / 01 file */
      if(strlen(argv[i])>5)
        p->fdet = argv[i]+5;
      else
        p->fdet = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, fdet=%s\n",argv[i],p->fdet);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"fobs=",5)){/** observable events / 01 file */
      if(strlen(argv[i])>5)
        p->fobs = argv[i]+5;
      else
        p->fobs = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, fobs=%s\n",argv[i],p->fobs);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"ferr=",5)){ /** errors / 01 file */
      if(strlen(argv[i])>5)
        p->ferr = argv[i]+5;
      else
        p->ferr = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, ferr=%s\n",argv[i],p->ferr);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"fer0=",5)){ /** errors / 01 file */
      if(strlen(argv[i])>5)
        p->fer0 = argv[i]+5;
      else
        p->fer0 = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, fer0=%s\n",argv[i],p->fer0);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"gdet=",5)){ /** generated detector events / 01 file */
      if(strlen(argv[i])>5)
        p->gdet = argv[i]+5;
      else
        p->gdet = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, gdet=%s\n",argv[i],p->gdet);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"gobs=",5)){/** generated observable events / 01 file */
      if(strlen(argv[i])>5)
        p->gobs = argv[i]+5;
      else
        p->gobs = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, gobs=%s\n",argv[i],p->gobs);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"pdet=",5)){ /** detector events / 01 file */
      if(strlen(argv[i])>5)
        p->pdet = argv[i]+5;
      else
        p->pdet = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, pdet=%s\n",argv[i],p->pdet);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"pobs=",5)){/** observable events / 01 file */
      if(strlen(argv[i])>5)
        p->pobs = argv[i]+5;
      else
        p->pobs = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, pobs=%s\n",argv[i],p->pobs);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if (0==strncmp(argv[i],"perr=",5)){ /** errors / 01 file */
      if(strlen(argv[i])>5)
        p->perr = argv[i]+5;
      else
        p->perr = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, perr=%s\n",argv[i],p->perr);
      if(p->mode>1)
	ERROR("mode=%d, this parameter %s is irrelevant\n",p->mode,argv[i]);
    }
    else if((strcmp(argv[i],"--help")==0)
            ||(strcmp(argv[i],"-h")==0)
            ||(strcmp(argv[i],"help")==0)){
      switch(p->mode){
      case 0:
	printf( HELP0 "\n" HELPU);
	break;	
      case 1:
	printf( HELP1 "\n" HELPU);
	break;	
      case 2:
	printf( HELP2);
	break;	
      case 3:
	printf( HELP3);
	break;	
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
    p->seed=(seed_old) + time(NULL)+1000000ul*getpid(); /* ensure a different seed even if started at the same time */
    if((p->debug)&&(p->mode!=3))
      printf("# initializing seed=%lld from time(NULL)+1000000ul*getpid()+%lld\n",p->seed, seed_old);
    /** use `tinymt64_generate_double(&pp.tinymt)` for double `[0,1)` */
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
      if(((p->fobs) || (p->fdet)) && ((p->perr == NULL) && (p->fer0 == NULL)))
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
  
  if(p->mA){
    if(p->mA->rows != p->nchk)
      ERROR("DET dimension mismatch nvar=%d and A %d by %d\n",p->nvar, p->mA->rows,p->mA->cols);
  }

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

  if(p->mulP > 0){/** scale probability values */
    double *ptr=p->vP;
    for(int i=0;i<p->nvar;i++, ptr++)
      *ptr *= p->mulP;    
  }
    
  LLR_table = init_LLR_tables(p->d1,p->d2,p->d3);

  p->dE = llr_from_dbl(p->dEdbl);
  //  p->uE = llr_from_dbl(p->uEdbl);
  
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
	printf("# use average LLR: aLLR = %g * aLLR + %g * LLR\n",p->bpgamma,1-p->bpgamma);
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

    if(p->perr){/** output predicted errors */
      p->file_perr = fopen(p->perr, "w");
      if(p->file_perr==NULL)
	ERROR("can't open the (predicted err) file %s for writing\n",p->perr);
    }
    if(p->pdet){/** output predicted syndrome (only really `makes sense for BP`) */
      p->file_pdet = fopen(p->pdet, "w");
      if(p->file_pdet==NULL)
	ERROR("can't open the (predicted det) file %s for writing\n",p->pdet);
    }
    if(p->pobs){/** output predicted observables */
      p->file_pobs = fopen(p->pobs, "w");
      if(p->file_pobs==NULL)
	ERROR("can't open the (predicted obs) file %s for writing\n",p->pobs);
    }
    if(p->gdet){/** output generated syndrome */
      p->file_gdet = fopen(p->gdet, "w");
      if(p->file_gdet==NULL)
	ERROR("can't open the (generated det) file %s for writing\n",p->gdet);
    }
    if(p->gobs){/** output generated observables */
      p->file_gobs = fopen(p->gobs, "w");
      if(p->file_gobs==NULL)
	ERROR("can't open the (generated obs) file %s for writing\n",p->gobs);
    }

    if((p->fer0)||(p->finA)){
      if((p->fer0==0)||(p->finA==NULL)||((p->fdet==0)&&(p->ferr==0))){
	printf("The parameters fer0=%s [e0] and finA=%s [A] must be set together, \n"
	       "also with either 'fdet' [s] or 'ferr' [e]; \n"
	       "to modify det evens as s'=s+A*e0 or s'=H*e+A*e0\n",p->fer0,p->finA);      
	ERROR("invalid input\n");
      }
    }
    
    if((p->ferr) && (p->fobs))
      ERROR("With ferr=%s cannot specify fobs=%s (will calculate)\n",p->fdem, p->fobs);
    if((p->ferr) && (p->fdet))
      ERROR("With ferr=%s cannot specify fdet=%s (will calculate)\n",p->fdem, p->fdet);    
    if((p->fdet==NULL)&&(p->fobs!=NULL))
      ERROR(" mode=%d fobs='%s' need detection events file 'fdet'\n",
	    p->mode, p->fobs);
    if ((p->fdet!=NULL)&&(p->fobs==NULL)&&(p->perr == NULL)&&(p->fer0 == NULL))
      ERROR(" mode=%d fdet='%s' need observables file 'fobs'\n",
	    p->mode, p->fdet);

    if(p->fer0){
      p->file_er0=fopen(p->fer0, "r");
      if(p->file_er0==NULL)
	ERROR("can't open the (er0) file %s for reading\n",p->fdet);
    }
    if(p->fdet){/** expect both `fdet` and `fobs` to be defined */
      p->file_det=fopen(p->fdet, "r");
      if(p->file_det==NULL)
	ERROR("can't open the (det) file %s for reading\n",p->fdet);
      if(p->fobs){
	p->file_obs=fopen(p->fobs, "r");
	if(p->file_obs==NULL)
	  ERROR("can't open the (obs) file %s for reading\n",p->fobs);
      }
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
    if (p->submode>3)
      ERROR(" mode=%d : submode='%d' unsupported\n",
	    p->mode, p->submode);
    break;
    
  case 3: /** read in DEM file and output the H, L, G, K matrices and P vector */
    if(strcmp(p->fout,"stdout")==0)
      p->use_stdout=1;
    if (p->submode>64)
      ERROR(" mode=%d : submode='%d' unsupported\n",
	    p->mode, p->submode);
    if(p->submode & 32){ // matrix modification
      if(!(p->submode & (31 | 64)))
	p->submode ^= (4+8+16);
      else 
	ERROR("mode=%d : submode=%d should equal 32 ('mode=3.32' for matrix modification)\n",p->mode,p->submode);
    }
    else if(p->submode & 64){ // write DEM model 
      /* no checks needed */
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
    if(p->fobs){
      p->file_obs=fopen(p->fobs, "r");
      if(p->file_obs==NULL)
	ERROR("can't open the (obs) file %s for reading\n",p->fobs);
    }
  }

  if(p->mode<=1){ /* vecdec RIS or BP decoder */
    p->mHe = mzd_init(p->nchk, p->nvec); /** each column a syndrome vector `H*e` */
    p->mLe = mzd_init(p->ncws,  p->nvec); /** each column `L*e` vector */
    if(p->ferr)
      p->mE = mzd_init(p->nvar, p->nvec);
    if(p->fer0)
      p->mE0 = mzd_init(p->mA->cols, p->nvec);
    if(p->debug &1){
      printf("# mode=%d, running %s decoder ",  p->mode, p->mode==0 ? "vecdec RIS" : "BP");
      switch(p->internal){
      case 0: printf(": use DET %s\n", p->fobs ? "and OBS files" : "file"); break;
      case 1: printf(": generating errors internally\n"); break;
      case 2: printf(": reading error vectors from ERR file\n"); break;
      default: ERROR("this should not happen");
	break;
      }
      if(p->fer0)
	printf("modifying DET events +A*e0 with A %s and e0 %s files\n",p->finA, p->fer0);
    }

    /** initialize the counters */
    for(extr_t i=0; i< EXTR_MAX; i++){
      cnt[i]=0;
      iter1[i]=0;
      iter2[i]=0;
    }

  }
  if(p->debug &1){
    printf("# using variables: mode=%d submode=%d\n",p->mode, p->submode);
    switch(p->mode){
    case 1:
      printf("# BP decoder parameters: bpgamma=%g bpretry=%d OSD level lerr=%d maxosd=%d\n",
	     p->bpgamma, p->bpretry, p->lerr,p->maxosd);
      /* fall through */      
    case 0:
      if(p->mode==0)
	printf("# RIS decoder parameters: swait=%d\n",p->swait);
      printf("# ntot=%lld nvec=%d steps=%d nfail=%lld lerr=%d\n",p->ntot,p->nvec,p->steps,p->nfail,p->lerr);
      break;
    case 2:
      printf("# Analyze small-weight codewords in %d RIS steps; swait=%d\n",p->steps,p->swait);
      printf("# maxC=%lld dE=%g dW=%d maxW=%d\n",p->maxC, dbl_from_llr(p->dE), p->dW, p->maxW);
      printf("# calculating %s%sminimum weight\n",p->submode&1 ? "est fail prob, ":"", p->submode&1 ? "exact fail prob, ":"");
      break;
    case 3:
      if(p->submode<32){
	printf("# Output %s%s%s%smatrices %s\n",
	       p->submode&1 ? "G=Hz " : "",
	       p->submode&2 ? "K=Lz " : "",
	       p->submode&4 ? "H=Hx " : "",
	       p->submode&8 ? "L=Lx " : "",
	       p->submode&16 ? "and probabilities vector" : "");
      }
      else if (p->submode &32)
	printf("# Do code modification\n");
      else if (p->submode &64)
	printf("# Write DEM file\n");
      break;
    default:
      ERROR("mode=%d not defined\n",p->mode);
    }
  }
  /** reset the counters */
  for (int i=0; i < EXTR_MAX; i++)
    cnt[i]=0;
  
  return 0;
};

/** @brief clean up variables and open files */
void var_kill(params_t *p){
  if(p->file_gdet)  fclose(p->file_gdet);
  if(p->file_gobs)  fclose(p->file_gobs);  
  if(p->file_pdet)  fclose(p->file_pdet);
  if(p->file_pobs)  fclose(p->file_pobs);  
  if(p->file_perr)  fclose(p->file_perr);  
  if(p->file_det)  fclose(p->file_det);
  if(p->file_obs)  fclose(p->file_obs);  
  if(p->file_err)  fclose(p->file_err);  
  if(p->file_er0)  fclose(p->file_er0);  
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
  if(p->v0) free(p->v0);
  if(p->v1) free(p->v1);
  if(p->svec) free(p->svec);
  if(p->err) free(p->err);
  if(p->ufl) ufl_free(p->ufl);
  if(p->hashU_syndr)
    kill_clusters(p);
}

int do_err_vecs(params_t * const p){

  int il1=p->nvec, il2;
  /** prepare error vectors ************************************************************/
  switch(p->internal){
  case 0: /** read `det` and `obs` files (each line a column) */
    il1=read_01(p->mHe,p->file_det, &p->line_det, p->fdet, 1, p->pads, p->debug);
    if(p->fobs){
      il2=read_01(p->mLe,p->file_obs, &p->line_obs, p->fobs, 1, 0, p->debug);
      if(il1!=il2)
	ERROR("mismatched DET %s (line %lld) and OBS %s (line %lld) files!",
	      p->fdet,p->line_det,p->fobs,p->line_obs);
      if(p->debug&1)
	printf("# read %d det/obs pairs\n",il1);
    }
    else
      if(p->debug&1)
	printf("# read %d det events\n",il1);
    if(p->fer0){
      il2=read_01(p->mE0,p->file_er0, &p->line_er0, p->fer0, 1, 0, p->debug);
      if(il1!=il2)
	ERROR("mismatched DET %s (line %lld) and ER0 %s (line %lld) files!",
	      p->fdet,p->line_det,p->fer0,p->line_er0);
      csr_mzd_mul(p->mHe,p->mA,p->mE0,0);
      if(p->debug&1)
	printf("updated %d det += A*e0 events\n",il1);
    }
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
    il1=read_01(p->mE,p->file_err, &p->line_err, p->ferr, 1, 0, p->debug);
    if(p->debug&1)
      printf("# read %d errors from file %s\n",il1,p->ferr);
    csr_mzd_mul(p->mHe,p->mH,p->mE,1);
    if(p->mL)
      csr_mzd_mul(p->mLe,p->mL,p->mE,1);
    if(p->fer0){
      il2=read_01(p->mE0,p->file_er0, &p->line_er0, p->fer0, 1, 0, p->debug);
      if(il1!=il2)
	ERROR("mismatched ERR %s (line %lld) and ER0 %s (line %lld) files!",
	      p->ferr,p->line_err,p->fer0,p->line_er0);
      csr_mzd_mul(p->mHe,p->mA,p->mE0,0);
      if(p->debug&1)
	printf("computed %d det events using H*e+A*e0\n",il1);
    }
    else if(p->debug&1)
      printf("computed %d det events H*e\n",il1);
    break;
  default:
    ERROR("internal=%d, this should not happen",p->internal);
  }
  if(p->gdet)
    mzd_write_01(p->file_gdet, p->mHe, 1, p->gdet, p->debug);
  if(p->gobs)
    mzd_write_01(p->file_gobs, p->mLe, 1, p->gobs, p->debug);
  
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
  mzd_t *mE0=NULL;

  switch(p->mode){
    long long int synd_fail; 
    int *status;                   
    mzd_t *srow;                  /** case 0, case 1 */
    qllr_t *ans;                  /** case 1 */
    size_t size;                  /** case 3 */
    char * comment;

  case 0: // mode=0 internal vecdec decoder
        // At least one round always
        synd_fail = 0;
        srow = mzd_init(1, p->nchk);

        if (p->submode & 1) { // submode 1
            // In submode1, we use Maximum Likelihood decoder
            if (!p->mG || !p->mK) {
                for (int i = 0; i < p->nchk; i++) {
                    p->vP[i] = P_from_llr(p->vLLR[i]);
                }
                p->rankH = rank_csr(p->mH);
                p->rankL = rank_csr(p->mL);
                p->rankG = p->nvar - p->rankH - p->rankL;
                if (p->debug & 1) {
                    printf("rankG=%d", p->rankG);
                }
                p->mG = do_G_matrix(p->mHt, p->mLt, p->vLLR, p->rankG, p->debug);
                int rankG = rank_csr(p->mG);
                if (p->debug & 1) {
                    printf("G matrix created with: rankG=%d \n", rankG);
                }
                p->mK = Lx_for_CSS_code(p->mG, p->mH);
                int rankK = rank_csr(p->mK);
                if (p->debug & 1) {
                    printf("K matrix created with: rankK=%d \n", rankK);
                }
            }
            for (long long int iround = 1; iround <= rounds; iround++) {
                if (p->debug & 1) {
                    printf("# starting round %lld of %lld", iround, rounds);
                    if (cnt[TOTAL])
                        printf(" pfail=%g fail=%lld out of total=%lld\n",
                               (double)(cnt[TOTAL] - cnt[SUCC_TOT]) / cnt[TOTAL],
                               cnt[TOTAL] - cnt[SUCC_TOT], cnt[TOTAL]);
                    else
                        printf("\n");
                }
            
                if (!(ierr_tot = do_err_vecs(p)))
                    break; // no more rounds
                if (p->uW >= 0) { // pre-decoding enabled
                    status = calloc(ierr_tot, sizeof(int)); // non-zero value = success pre_dec
                    if (!status) ERROR("memory allocation");
                    p->mHeT = mzd_transpose(p->mHeT, p->mHe);
                    mE0 = mzd_init(p->nvec, p->nvar);
                    long long int cnt_pre = 0;
                    for (long long int ierr = 0; ierr < ierr_tot; ierr++) { // cycle over errors
                        mzd_copy_row(srow, 0, p->mHeT, ierr); // syndrome row in question       
                        int res_pre = dec_ufl_one(srow, p);
                        if (res_pre) { // pre-decoder success
                            mzd_row_add_vec(mE0, ierr, p->ufl->error, 1);
                            status[ierr] = res_pre;
                            cnt_pre++;
                        }
                    }
                    if (cnt_pre < ierr_tot) { // some pre-decoder failures
                        long long int num = ierr_tot - cnt_pre;
                        mzd_t *mST = mzd_init(num, p->nchk);
                        for (long long int ierr = 0, row = 0; ierr < ierr_tot; ierr++) {
                            if (!status[ierr])
                                mzd_copy_row(mST, row++, p->mHeT, ierr);
                        }
                        mzd_t *mS = mzd_transpose(NULL, mST);
                        mzd_t *mE2 = do_decode_BAR(mS, p); // each row a decoded error vector
                        for (long long int ierr = 0, row = 0; ierr < ierr_tot; ierr++) {
                            if (!status[ierr]) {
                                mzd_copy_row(mE0, ierr, mE2, row++);
                                status[ierr] = 4;
                                cnt[CONV_RIS]++;
                            }
                        }
                        mzd_free(mE2);
                        mzd_free(mS);
                        mzd_free(mST);
                    }
                } else { // no pre-decoding      
                    status = NULL;
                    // Actually decode and generate error vectors
                    #ifndef NDEBUG  // need `mHe` later
                    mzd_t *mS = mzd_copy(NULL, p->mHe);
                    mE0 = do_decode_BAR(mS, p); // each row a decoded error vector
                    mzd_free(mS);
                    mS = NULL;
                    #else
                    mE0 = do_decode_BAR(p->mHe, p); // each row a decoded error vector
                    #endif /* NDEBUG */
                    cnt[CONV_RIS] += ierr_tot;
                }  
                mzd_t *mE0t = mzd_transpose(NULL, mE0);
                mzd_free(mE0);
                mE0 = NULL;
                
                #ifndef NDEBUG
                mzd_t *prodHe = csr_mzd_mul(NULL, p->mH, mE0t, 1);
                if (p->pdet)
                    mzd_write_01(p->file_pdet, prodHe, 1, p->pdet, p->debug);
                
                mzd_add(prodHe, prodHe, p->mHe);
                if ((p->steps) && (!mzd_is_zero(prodHe))) {
                    if ((p->debug & 512) || (p->nvec <= 64)) {
                        printf("syndromes difference:\n");
                        mzd_print(prodHe);
                    }
                    if (p->steps > 0) // otherwise we do not care
                        ERROR("some syndromes are not matched!\n");
                }
                mzd_free(prodHe);
                prodHe = NULL;
                #else /* NDEBUG defined */
                if (p->pdet) {
                    mzd_t *prodHe = csr_mzd_mul(NULL, p->mH, mE0t, 1);
                    mzd_write_01(p->file_pdet, prodHe, 1, p->pdet, p->debug);
                    mzd_free(prodHe);
                }
                #endif

                if (p->perr)
                    mzd_write_01(p->file_perr, mE0t, 1, p->perr, p->debug);
                mzd_t *prodLe = csr_mzd_mul(NULL, p->mL, mE0t, 1);
                if (p->pobs)
                    mzd_write_01(p->file_pobs, prodLe, 1, p->pobs, p->debug);

                mzd_add(prodLe, prodLe, p->mLe);
                int fails = 0;
                if (p->uW >= 0) { // pre-decoding enabled
                    assert(status);
                    rci_t j = 0;
                    for (rci_t ic = 0; ic < ierr_tot; ic++) {
                        rci_t ir = 0;
                        if (mzd_find_pivot(prodLe, ir, ic, &ir, &ic)) {
                            cnt[SUCC_TOT] += ic - j;
                            while (j < ic) { // success
                                switch (status[j++]) {
                                    case 1: cnt[SUCC_TRIVIAL]++; break;
                                    case 2: cnt[SUCC_LOWW]++; break;
                                    case 3: cnt[SUCC_CLUS]++; break;
                                    case 4: cnt[SUCC_RIS]++; break;
                                    default: ERROR("unexpected");
                                }
                            }
                            j++;
                            fails++;
                        } else // no more pivots
                            break;            
                    }
                    cnt[SUCC_TOT] += ierr_tot - j;      
                    while (j < ierr_tot) { // success
                        switch (status[j++]) {
                            case 1: cnt[SUCC_TRIVIAL]++; break;
                            case 2: cnt[SUCC_LOWW]++; break;
                            case 3: cnt[SUCC_CLUS]++; break;
                            case 4: cnt[SUCC_RIS]++; break;
                            default: ERROR("unexpected");
                        }
                    }
                } else { // no pre-decoding
                    for (rci_t ic = 0; ic < ierr_tot; ic++) {
                        rci_t ir = 0;
                        if (mzd_find_pivot(prodLe, ir, ic, &ir, &ic))
                            fails++;    
                        else // no more pivots
                            break;
                    }      
                    cnt[SUCC_RIS] += ierr_tot - fails;
                    cnt[SUCC_TOT] += ierr_tot - fails;
                }
                // Update the global counts
                synd_fail += fails;
                cnt[TOTAL] += ierr_tot;
                mzd_free(prodLe);
                prodLe = NULL;
                if (status) {
                    free(status);
                    status = NULL;
                }
                if ((p->nfail > 0) && (synd_fail >= p->nfail))
                    break;
            }
        } else {
            // Original code for non-submode 1
            for (long long int iround = 1; iround <= rounds; iround++) {
                if (p->debug & 1) {
                    printf("# starting round %lld of %lld", iround, rounds);
                    if (cnt[TOTAL])
                        printf(" pfail=%g fail=%lld out of total=%lld\n",
                               (double)(cnt[TOTAL] - cnt[SUCC_TOT]) / cnt[TOTAL],
                               cnt[TOTAL] - cnt[SUCC_TOT], cnt[TOTAL]);
                    else
                        printf("\n");
                }
            
                if (!(ierr_tot = do_err_vecs(p)))
                    break; // no more rounds
                if (p->uW >= 0) { // pre-decoding enabled
                    status = calloc(ierr_tot, sizeof(int)); // non-zero value = success pre_dec
                    if (!status) ERROR("memory allocation");
                    p->mHeT = mzd_transpose(p->mHeT, p->mHe);
                    mE0 = mzd_init(p->nvec, p->nvar);
                    long long int cnt_pre = 0;
                    for (long long int ierr = 0; ierr < ierr_tot; ierr++) { // cycle over errors
                        mzd_copy_row(srow, 0, p->mHeT, ierr); // syndrome row in question       
                        int res_pre = dec_ufl_one(srow, p);
                        if (res_pre) { // pre-decoder success
                            mzd_row_add_vec(mE0, ierr, p->ufl->error, 1);
                            status[ierr] = res_pre;
                            cnt_pre++;
                        }
                    }
                    if (cnt_pre < ierr_tot) { // some pre-decoder failures
                        long long int num = ierr_tot - cnt_pre;
                        mzd_t *mST = mzd_init(num, p->nchk);
                        for (long long int ierr = 0, row = 0; ierr < ierr_tot; ierr++) {
                            if (!status[ierr])
                                mzd_copy_row(mST, row++, p->mHeT, ierr);
                        }
                        mzd_t *mS = mzd_transpose(NULL, mST);
                        mzd_t *mE2 = do_decode(mS, p); // each row a decoded error vector
                        for (long long int ierr = 0, row = 0; ierr < ierr_tot; ierr++) {
                            if (!status[ierr]) {
                                mzd_copy_row(mE0, ierr, mE2, row++);
                                status[ierr] = 4;
                                cnt[CONV_RIS]++;
                            }
                        }
                        mzd_free(mE2);
                        mzd_free(mS);
                        mzd_free(mST);
                    }
                } else { // no pre-decoding      
                    status = NULL;
                    // Actually decode and generate error vectors
                    #ifndef NDEBUG  // need `mHe` later
                    mzd_t *mS = mzd_copy(NULL, p->mHe);
                    mE0 = do_decode(mS, p); // each row a decoded error vector
                    mzd_free(mS);
                    mS = NULL;
                    #else
                    mE0 = do_decode(p->mHe, p); // each row a decoded error vector
                    #endif /* NDEBUG */
                    cnt[CONV_RIS] += ierr_tot;
                }  
                mzd_t *mE0t = mzd_transpose(NULL, mE0);
                mzd_free(mE0);
                mE0 = NULL;
                
                #ifndef NDEBUG
                mzd_t *prodHe = csr_mzd_mul(NULL, p->mH, mE0t, 1);
                if (p->pdet)
                    mzd_write_01(p->file_pdet, prodHe, 1, p->pdet, p->debug);
                
                mzd_add(prodHe, prodHe, p->mHe);
                if ((p->steps) && (!mzd_is_zero(prodHe))) {
                    if ((p->debug & 512) || (p->nvec <= 64)) {
                        printf("syndromes difference:\n");
                        mzd_print(prodHe);
                    }
                    if (p->steps > 0) // otherwise we do not care
                        ERROR("some syndromes are not matched!\n");
                }
                mzd_free(prodHe);
                prodHe = NULL;
                #else /* NDEBUG defined */
                if (p->pdet) {
                    mzd_t *prodHe = csr_mzd_mul(NULL, p->mH, mE0t, 1);
                    mzd_write_01(p->file_pdet, prodHe, 1, p->pdet, p->debug);
                    mzd_free(prodHe);
                }
                #endif

                if (p->perr)
                    mzd_write_01(p->file_perr, mE0t, 1, p->perr, p->debug);
                mzd_t *prodLe = csr_mzd_mul(NULL, p->mL, mE0t, 1);
                if (p->pobs)
                    mzd_write_01(p->file_pobs, prodLe, 1, p->pobs, p->debug);

                mzd_add(prodLe, prodLe, p->mLe);
                int fails = 0;
                if (p->uW >= 0) { // pre-decoding enabled
                    assert(status);
                    rci_t j = 0;
                    for (rci_t ic = 0; ic < ierr_tot; ic++) {
                        rci_t ir = 0;
                        if (mzd_find_pivot(prodLe, ir, ic, &ir, &ic)) {
                            cnt[SUCC_TOT] += ic - j;
                            while (j < ic) { // success
                                switch (status[j++]) {
                                    case 1: cnt[SUCC_TRIVIAL]++; break;
                                    case 2: cnt[SUCC_LOWW]++; break;
                                    case 3: cnt[SUCC_CLUS]++; break;
                                    case 4: cnt[SUCC_RIS]++; break;
                                    default: ERROR("unexpected");
                                }
                            }
                            j++;
                            fails++;
                        } else // no more pivots
                            break;            
                    }
                    cnt[SUCC_TOT] += ierr_tot - j;      
                    while (j < ierr_tot) { // success
                        switch (status[j++]) {
                            case 1: cnt[SUCC_TRIVIAL]++; break;
                            case 2: cnt[SUCC_LOWW]++; break;
                            case 3: cnt[SUCC_CLUS]++; break;
                            case 4: cnt[SUCC_RIS]++; break;
                            default: ERROR("unexpected");
                        }
                    }
                } else { // no pre-decoding
                    for (rci_t ic = 0; ic < ierr_tot; ic++) {
                        rci_t ir = 0;
                        if (mzd_find_pivot(prodLe, ir, ic, &ir, &ic))
                            fails++;    
                        else // no more pivots
                            break;
                    }      
                    cnt[SUCC_RIS] += ierr_tot - fails;
                    cnt[SUCC_TOT] += ierr_tot - fails;
                }
                // Update the global counts
                synd_fail += fails;
                cnt[TOTAL] += ierr_tot;
                mzd_free(prodLe);
                prodLe = NULL;
                if (status) {
                    free(status);
                    status = NULL;
                }
                if ((p->nfail > 0) && (synd_fail >= p->nfail))
                    break;
            }
        }


        if (!((p->fdet) && (p->fobs == NULL) && (p->perr))) { // except in the case of partial decoding
            if (p->steps > 0) {  // otherwise results are invalid as we assume syndromes to match
                cnt_out(p->debug & 1, p);
            }   
        } else if (p->debug & 1) {
            printf("# all finished\n");
        }

        break;

  case 1: /** `mode=1` various BP flavors */    
    ans = calloc(p->nvar, sizeof(qllr_t));
    if(!ans) ERROR("memory allocation failed!"); 

    const int do_file_output = (p->perr) || (p->pdet) || (p->pobs);
    mzd_t *pErr=NULL, *pHerr=NULL, *pLerr=NULL;
    srow=mzd_init(1,p->nchk);
    if(do_file_output){
      pErr=mzd_init(1,p->nvar);
      if(p->pdet)
	pHerr=mzd_init(1,p->mHt->cols);
      if(p->pobs)
	pLerr=mzd_init(1,p->mLt->cols);
    }
    
    for(long long int iround=1; iround <= rounds; iround++){
      if(p->debug&1){
	printf("# starting round %lld of %lld pfail=%g fail=%lld total=%lld\n", iround, rounds,
	       cnt[TOTAL] ? (double) (cnt[TOTAL]-cnt[SUCC_TOT])/cnt[TOTAL] : 0.5, cnt[TOTAL]-cnt[SUCC_TOT], cnt[TOTAL]);
	fflush(stdout);
      }
    
      if( !(ierr_tot = do_err_vecs(p)))
	break; /** no more rounds */
      p->mHeT = mzd_transpose(p->mHeT,p->mHe);
      p->mLeT = mzd_transpose(p->mLeT,p->mLe);
      for(long long int ierr = 0; ierr < ierr_tot; ierr++){ /** cycle over errors */
	cnt[TOTAL]++;
#ifndef NDEBUG	  
	if((p->debug&8)&&(p->nvar <= 256)&&(p->debug&512)){
	  printf("non-trivial error %lld of %lld:\n",ierr+1,ierr_tot);
	  if(p->mE) /** print column as row */	      
	    for(int i=0; i<p->nvar; i++)
	      printf("%s%d%s",i==0?"[":" ",mzd_read_bit(p->mE,i,ierr),i+1<p->nvar?"":"]\n");
	  mzd_print_row(p->mHeT,ierr);
	  mzd_print_row(p->mLeT,ierr);
	  out_llr("i",p->nvar,p->vLLR);
	}
#endif 	  
	//	mzd_t * const srow = mzd_init_window(p->mHeT, ierr,0, ierr+1,p->nchk); /* syndrome row */
	//! TODO: why does the above fail? 
	mzd_copy_row(srow,0,p->mHeT,ierr); /** syndrome row in question */
	int res_pre=0;
	if(p->uW >=0)
	  res_pre = dec_ufl_one(srow,p);
	if(res_pre){ /** pre-decoder success */
	  if(p->perr) 
	    write_01_vec(p->file_perr, p->ufl->error, p->mH->cols, p->perr); /** error prediction */
	  if(p->pdet){	      
	    write_01_vec(p->file_pdet, p->ufl->syndr, p->mH->rows, p->pdet); /** syndrome prediction */
	  }
	  if(!p->err)
	    p->err=vec_init(p->nvar); /** temporary storage */
	  if(!p->obs)
	    p->obs=vec_init(p->mL->rows);
	  assert(p->nvar >= p->mL->rows);
	  vec_t *vobs = csr_vec_mul(p->err, p->obs, p->mLt, p->ufl->error, 1);
	  if(p->pobs){
	    write_01_vec(p->file_pobs, vobs, p->mL->rows, p->pobs); /** obs prediction */
	  }
	  //	  cnt_update(CONV_TRIVIAL,0); /** trivial convergence after `0` steps */
	  if(mzd_row_vec_match(p->mLeT,ierr,vobs)){
	      cnt[SUCC_TOT]++;
	      switch(res_pre){
	      case 1: cnt[SUCC_TRIVIAL]++; break;
	      case 2: cnt[SUCC_LOWW]++; break;
	      case 3: cnt[SUCC_CLUS]++; break;
	      default: ERROR("unexpected"); break;
	      }
	    continue ; /** next error / syndrome vector pair */     	    
	  }
	}
	else{ /** failed `pre`, do actual BP */
	  cnt[NUMB_BP]++;
	  int conv = do_dec_bp_one(ans,srow,p);
	  int conv_BP = conv>>1, conv_OSD = conv%2;
	  if(do_file_output){ /** output predicted values */
	    mzd_row_clear_offset(pErr,0,0);
	    for(int i=0; i < p->nvar; i++)
	      if(ans[i]<0)
		mzd_flip_bit(pErr,0,i);
	    if(p->perr)
	      mzd_write_01(p->file_perr, pErr, 0, p->perr, p->debug);
	    if(p->pdet){
	      mzd_row_csr_mul(pHerr,0, pErr,0, p->mHt, 1);
	      mzd_write_01(p->file_pdet, pHerr, 0, p->pdet, p->debug);
	    }
	    if(p->pobs){
	      mzd_row_csr_mul(pLerr,0, pErr,0, p->mLt, 1);
	      mzd_write_01(p->file_pobs, pLerr, 0, p->pobs, p->debug);
	    }
	  }
	  if (!((p->fdet)&&(p->fobs==NULL)&&(p->perr))){ /** except the case of partial decoding */  
	    if((conv_BP)||(conv_OSD)){ /** `convergence success`  */	      
	      mzd_t * const obsrow = mzd_init_window(p->mLeT, ierr,0, ierr+1,p->mLeT->ncols);
	      if(syndrome_check(ans, obsrow, p->mL, p)){
		cnt[SUCC_TOT]++;
		if(conv_BP)
		  cnt[SUCC_BP]++;
		else
		  cnt[SUCC_OSD]++;
	      }
	      mzd_free(obsrow);
	    }	  
	  }
	}
	//	mzd_free(srow);
	if((p->nfail) && cnt[TOTAL]-cnt[SUCC_TOT] >= p->nfail)
	  break;
      }
    }    
    if(p->debug&1)
      ufl_cnt_print(p);
    cnt_out(p->debug&1,p);
    if(srow)
      mzd_free(srow);
    if(pErr)
      mzd_free(pErr);
    if(pHerr)
      mzd_free(pHerr);
    if(pLerr)
      mzd_free(pLerr);
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
    //    do_hash_fail_prob(p); /** output results */
    if(p->steps){
      /** TODO: make more intelligent distance estimate or param */
      int dx= p->fdem == NULL ? 1 : 3 ;      
      if (p->debug&1)
	printf("# try to remove reducible codewords dx=%d ...\n",dx);
      do_hash_remove_reduc(dx,p);
      if (p->debug&1)
	printf("# done\n");
    }
    do_hash_fail_prob(p); /** output results */    
    if(p->outC){
      char * name;
      if(p->fdem)
	name=p->fdem;
      else if (p->finH)
	name=p->finH;
      else
	name="(unknown source)";
      size_t size = 1 + snprintf(NULL, 0, "# codewords computed from '%s', maxW=%d ", name, p->maxW);
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
    if(p->submode&64){ /** write DEM */
      if((p->classical)||(!(p->mL)))
	ERROR("mode=%d submode=64 (bit 6 set) must be a quantum code",p->mode);
      if(!p->mHt)
	p->mHt = csr_transpose(NULL,p->mH);
      if(!p->mLt)
	p->mLt = csr_transpose(NULL,p->mL);
      size_t size = snprintf(NULL, 0, "DEM model from H file %s", p->finH);
      if(!(comment = malloc(size + 1)))
	ERROR("memory allocation");
      sprintf(comment, "DEM model from H file %s", p->finH);

      write_dem_file(p->fout,"D.dem",p->mHt, p->mLt, p->vP, comment);
      free(comment);
      break;
    }
    else if(p->submode&32){ /** matrix transformation */
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
	  if(p->finC){
	    p->mG = do_G_from_C(p->mLt,p->codewords,p->rankG,p->minW, p->maxW ? p->maxW : p->minW + p->dW, p->debug);
	  }
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
	    p->mK = do_K_from_C(p->mLt, p->codewords, p->ncws, p->nvar,
				p->minW, p->maxW ? p->maxW : p->minW+p->dW, p->debug);
	    //	  csr_out(p->mK);
	  }
	  else{
	    if ((p->mG)&&(p->mH)){
	      p->mK=Lx_for_CSS_code(p->mG,p->mH);
	    }
	    else
	      ERROR("use mode=3.3 or specify the codewords file `finC=...` or generator matrix `finG=...`");
	  }
	  comment[0]='K';
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

    if((p->debug &16) &&(p->mG) && (p->vLLR)){ /** calculate partition functions */
      printf("debug&16: calculating partition functions (may be slow)\n");
      mzd_t * err = mzd_init(1,p->nvar);
      double prob0 = do_Z(p->vLLR, err, p->mG, p->debug);
      printf("Z(0)=%g\n",prob0);
      if(p->mK){
	int max = p->ncws > 2 ? 1<<2 : 1 << p->ncws ; 
	uint32_t prev=0;
	for(int bin=1; bin < max; bin++){
	  uint32_t gray =binary_to_gray(bin);
	  uint32_t idx = getMsb(gray^prev); /** row to add */
	  for(int i=p->mK->p[idx]; i < p->mK->p[idx+1]; i++){
	    int j = p->mK->i[i]; /** position */
	    mzd_flip_bit(err,0,j);      
	  }	  
	  printf("Z(%d)=",gray); fflush(stdout);
	  double prob = do_Z(p->vLLR, err, p->mG, p->debug);
	  printf("%g ratio=%g\n", prob, prob/prob0);
	  prev=gray;
	}
      }
      mzd_free(err);
    }
    break;
    
  default:
    ERROR("mode=%d not supported\n",p->mode);
    break;
  }

  var_kill(p);
  return 0;
}
