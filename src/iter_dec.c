/**
 *  @file iter_dec.c
 *
 * @brief iter_dec - iterative decoder for quantum and classical codes
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
#include <m4ri/m4ri.h>
#include "utils.h"
#include "util_m4ri.h"
#include "qllr.h"
#include "vecdec.h"


void cnt_out(int print_banner){
  if(print_banner)
    printf("# FAIL_FRAC TOTAL  C_TRIVIAL S_TRIVIAL  C_BP C_BP_AVG C_BP_TOT S_BP   S_OSD S_TOT\n");
  printf(" %10g %lld  %lld %lld      %lld %lld %lld %lld  %lld %lld\n",
	 (double ) (cnt[TOTAL]-cnt[SUCC_TOT])/cnt[TOTAL],cnt[TOTAL],
	 cnt[CONV_TRIVIAL], cnt[SUCC_TRIVIAL],
	 cnt[CONV_BP], cnt[CONV_BP_AVG], cnt[CONV_BP_TOT], cnt[SUCC_BP],
	 cnt[SUCC_OSD], cnt[SUCC_TOT]);
}

/** @brief print entire `one_vec_t` structure by pointer */
void print_one_vec(const one_vec_t * const pvec){
  printf(" w=%d E=%g cnt=%d [",pvec->weight, dbl_from_llr(pvec->energ),pvec->cnt);
  for(int i=0; i < pvec->weight; i++)
    printf("%d%s",pvec->arr[i], i+1 < pvec->weight ? " " :"]\n");
}

void cnt_update(extr_t which, int istep){
  const long long int max = (LLONG_MAX >> 3);
  cnt[which]++;
  if((which == CONV_BP) || (which == CONV_BP_AVG))
    cnt[CONV_BP_TOT]++;
  iter1[which]+=abs(istep);
  iter2[which]+=istep*istep;
  if((cnt[which] > max) ||(iter1[which] > max) || (iter2[which]>max)){
    cnt_out(1);
    ERROR("too many steps\n"); 
  }
}

/** @brief helper function for binary search */
int cmpfunc(const void * a, const void * b) {
   return ( *(int *) a - *(int *) b );
}

void out_llr(const char str[], const int num, const qllr_t llr[]){
  const int max= num > 20 ? 20 : num-1 ;
  printf("%s[",str);
  for(int iv = 0; iv <= max; iv++)
    printf("%+6.1g%s", dbl_from_llr(llr[iv]), iv< max ? " " : "]\n");
}

/** @brief Check the syndrome for the LLR vector given (`0`th row of `syndrome`).
 * @return 1 if codeword is valid, 0 otherwise */
int syndrome_check(const qllr_t LLR[], const mzd_t * const syndrome,
		   const csr_t * const H,
		   _maybe_unused const params_t * const p){
  const int nchk = H->rows;
  //  const csr_t * const H = p->mH;
  for(int ic=0; ic<nchk; ic++){ /** target `check` node */
    int synd = mzd_read_bit(syndrome,0,ic);
    for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
      const int iv = H->i[j]; /** origin `variable` node index */
      if(LLR[iv] < 0)
	synd ^= 1;
    }
    if(synd)
      return 0; /** invalid codeword */
  }
  return 1; /** valid codeword */
}

/** 
 * NOTE: `message indexing conventions` 
 * total number of messages the same as non-zero elements in `H` (or `Ht`) matrix.
 * for `V->C` messages: `j` is index of element `v` in row `c` of `H`
 * for `C->V` messages: `i` is index of element `c` in row `v` of `Ht`.
 */ 

/** @brief for given `ic` compute V to C messages.  WARNING: `xLLR` need to be updated just before run */
static inline void bp_do_VCc_opt(const int ic, qllr_t * const mVtoC, const qllr_t * const mCtoV,
			     const qllr_t xLLR[], 
			     const csr_t * const H, const csr_t * const Ht){  
  for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
    const int iv = H->i[j]; /** origin `variable` node index */
    const int j0 = Ht->p[iv];
    const int *beg = &(Ht->i[j0]);
    const int *pos = (int *) bsearch (& ic, beg, Ht->p[iv+1] - j0, sizeof(int), cmpfunc);
    int j1 = (pos-beg) + j0;
    assert(Ht->i[j1] == ic);
    mVtoC[j] = xLLR[iv] - mCtoV[j1];
  }
}

/** @brief for given `ic` compute V to C messages.  
    This version only uses bare LLR values */
static inline void bp_do_VCc(const int ic, qllr_t * const mVtoC, const qllr_t * const mCtoV,
			     const csr_t * const H, const csr_t * const Ht, const qllr_t LLR[]){  
  for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
    const int iv = H->i[j]; /** origin `variable` node index */
    qllr_t msg = LLR[iv]; /** `bare LLR` need an extra parameter */
    for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++){
      const int ic1 = Ht->i[j1];/** aux `check` node index */
      if(ic1 != ic)
	msg += mCtoV[j1];
    }
    mVtoC[j] = msg;
  }
}

/** @brief for given `iv` compute V to C messages.  
    This version only uses bare LLR values */
static inline void bp_do_VCv(const int iv, qllr_t * const mVtoC, const qllr_t * const mCtoV,
			     const csr_t * const H, const csr_t * const Ht, const qllr_t LLR){
  for(int jv = Ht->p[iv]; jv < Ht->p[iv+1]; jv++){
    const int ic = Ht->i[jv]; /** destination `check` node index */

    const int j0 = H->p[ic];
    const int *beg = &(H->i[j0]);
    const int *pos = (int *) bsearch (& iv, beg, H->p[ic+1] - j0, sizeof(int), cmpfunc);
    int j = (pos-beg) + j0;
    assert((iv == H->i[j]) && (j< H->p[ic+1]));    
    qllr_t msg = LLR; /** `bare LLR` at `iv` need an extra parameter */
    for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++){
      const int ic1 = Ht->i[j1];/** aux `check` node index */
      if(ic1 != ic)
	msg += mCtoV[j1];
    }
    mVtoC[j] = msg;
  }
}


/** @brief init all V->C messages to bare LLR */
static inline void bp_init_VC(qllr_t mesVtoC[], const csr_t * const H, const qllr_t LLR[]){
  const int nchk = H->rows;
  for(int ic=0; ic<nchk; ic++){ /** target `check` node index */
    for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
      const int iv = H->i[j]; /** `variable` node index */
      mesVtoC[j]= LLR[iv];   /** `j` is the edge (message) index */
    }
  }
}

/** @brief for given `iv` compute C to V messages */
static inline void bp_do_CVv(const int iv, const qllr_t * const mVtoC, qllr_t * const mCtoV, 
		     const csr_t * const H, const csr_t * const Ht, 
		     const mzd_t * const srow){  
  for(int j = Ht->p[iv]; j < Ht->p[iv+1] ; j++){
    /** TODO: optimize this loop as in `LDPC_Code::bp_decode()` of `itpp` package */
    const int ic = Ht->i[j];/** origin `check` node index */
    int sbit = mzd_read_bit(srow,0,ic); /** syndrome bit */
    qllr_t msg=0; 
    int new=1;
    for(int j1 = H->p[ic]; j1 < H->p[ic+1] ; j1++){
      const int iv1 = H->i[j1]; /** aux `variable` node index */
      if(iv1!=iv){
	if(new){
	  msg = mVtoC[j1];
	  new=0;
	}
	else
	  msg = boxplus(msg, mVtoC[j1]);
      }
    }
    mCtoV[j] = (sbit ? -msg : msg);
  }
}

/** @brief for given `ic` compute C to V messages */
static inline void bp_do_CVc(const int ic, const qllr_t * const mVtoC, qllr_t * const mCtoV, 
		     const csr_t * const H, const csr_t * const Ht, const int sbit){  
  //  int sbit = mzd_read_bit(srow,0,ic); /** syndrome bit */
  /** TODO: optimize this loop as in `LDPC_Code::bp_decode()` of `itpp` package (`maybe`) */
  for(int jv = H->p[ic]; jv < H->p[ic+1] ; jv++){
    int iv=H->i[jv]; /** initial `v` node */
      
    /** find j: `Ht->p[iv]` <= `j` < `Ht->p[iv+1]` s.t. `Ht->i[j]`==`ic` */
    const int j0 = Ht->p[iv];
    const int *beg = &(Ht->i[j0]);
    const int *pos = (int *) bsearch (& ic, beg, Ht->p[iv+1] - j0, sizeof(int), cmpfunc);
    int j = (pos-beg) + j0;
    assert(Ht->i[j] == ic);
    
    qllr_t msg=0; 
    int new=1;
    for(int j1 = H->p[ic]; j1 < H->p[ic+1] ; j1++){
      const int iv1 = H->i[j1]; /** `variable` node summation index */
      if(iv1!=iv){
	if(new){
	  msg = mVtoC[j1];
	  new=0;
	}
	else
	  msg = boxplus(msg, mVtoC[j1]);
      }
    }
    mCtoV[j] = (sbit ? -msg : msg);
  }
}


/** @brief for given `iv` compute expected LLR */
static inline qllr_t bp_do_LLR(const int iv, const qllr_t * const mCtoV, const csr_t * const Ht, const qllr_t LLR){  
  qllr_t val= LLR; /** WARNING: do we need this? */
  for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++)
    val += mCtoV[j1];
  
  return val; /** need this in any case */
}


/** @brief baseline parallel BP algorithm (sum-product or max-product depending on LLRs used).  
 * 
 * Decode `e` vector using syndrome bits in `srow = e*Ht` and apriori
 * bit LLR values in `LLR`.  Assume `srow` non-zero (no check for a
 * trivial solution).  Non-zero return value: (+ number of steps) for regular
 * LLR, and (- number of steps) for average LLR.
 * 
 * @param[out] outLLR the resulting LLR (whether converged or not)
 * @return 0 if failed to converge or (+/- the number of steps) if converged successfully. 
 *  
 */
int do_parallel_BP(qllr_t * outLLR, const mzd_t * const srow,
		   const csr_t * const H, const csr_t * const Ht,
		   const qllr_t LLR[], const params_t * const p){
  assert(outLLR != NULL);
  _maybe_unused const int nchk = p->nchk; /** number of check nodes */
  const int nvar = p->nvar;     /** number of variable nodes */
  const int nz = p->mH->p[p->nchk]; /** number of messages */
  qllr_t *mesVtoC = malloc(nz*sizeof(qllr_t));
  qllr_t *mesCtoV = malloc(nz*sizeof(qllr_t));
  if((!mesVtoC)||(!mesCtoV)) ERROR("memory allocation failed");
  qllr_t * xLLR=NULL, *aLLR=NULL;
  int succ_BP=0;
  /** convenience setting: by default, use both LLR and average LLR */
  const int submode = ((p->submode & 3) == 0) ? (p->submode | 3) : p->submode;

  if((submode & 3) == 1) /** use xLLR only */
    xLLR = outLLR; /** return soft vector of LLR here */
  else{                /** use aLLR or both LLR and aLLR */
    aLLR = outLLR; 
    xLLR = calloc(nvar,sizeof(qllr_t)); /** needed for calculations in any case */
    if(!xLLR) ERROR("memory allocation failed");
  }
  
  /** init V->C messages to bare LLR */
  bp_init_VC(mesVtoC,H,LLR);

  for (int istep=1; istep <= p->steps  ; istep++){ /** main decoding cycle */
    /** C -> V messages */
#if 1 
    for(int iv=0; iv<nvar; iv++) /** target `variable` node index */
      bp_do_CVv(iv, mesVtoC, mesCtoV, H, Ht, srow);
#else
    for(int ic=0; ic<nchk; ic++) /** source `check` node index */
      bp_do_CVc(ic, mesVtoC, mesCtoV, H, Ht, mzd_read_bit(srow,0,ic));
#endif     
            
    for(int iv = 0; iv< nvar; iv++) /** calculate expected LLR for variable nodes */
      xLLR[iv] = bp_do_LLR(iv,mesCtoV,Ht,LLR[iv]);
    
    if(submode&2) /** use average LLR */
      for(int iv = 0; iv< nvar; iv++)
	aLLR[iv] = llr_from_dbl(p->bpalpha * dbl_from_llr(aLLR[iv]) +(1.0 - p->bpalpha) * dbl_from_llr(xLLR[iv])); 

#ifndef NDEBUG    
    if(p->debug & 8){
      if(submode&1) /** use regular LLR */
	out_llr("x",p->nvar, xLLR);
      if(submode&2)
	out_llr("a",p->nvar, aLLR);
    }
#endif
    /** do convergence check */
    if((submode&2) && (syndrome_check(aLLR,srow,p->mH, p))){ 
      cnt_update(CONV_BP_AVG, istep);
      succ_BP=-istep;
      break;
    }
    else if((submode&1) && syndrome_check(xLLR,srow,p->mH, p)){
      cnt_update(CONV_BP, istep);
      succ_BP=istep;
      if(submode&2)/** need to copy */
	for(int iv=0; iv < p->nvar; iv++)
	  outLLR[iv] = xLLR[iv];
      break;
    }
    /* V -> C messages */ 
#if 0
    for(int ic=0; ic<nchk; ic++) /** target `check` node */
      bp_do_VCc_opt(ic, mesVtoC, mesCtoV, xLLR, H, Ht);
#else 
    for(int iv=0; iv<nvar; iv++) /** source  `variable` node */
      bp_do_VCv(iv, mesVtoC, mesCtoV, H, Ht, LLR[iv]);
#endif     
    
  }
  
  /** clean up ********************/
  if(submode&2)
    free(xLLR);
  free(mesVtoC);
  free(mesCtoV);

  /** default value is returned */

  return succ_BP;
}

/** @brief serial-C BP algorithm (sum-product or max-product depending on LLRs used).  
 * 
 * Same as `do_parallel_BP()` function above but using `check-based` serial order: 
 * for each `ic` do:
 *    update V->C messages to `ic`;
 *    update C->V messages from `ic`.
 * 
 * @param[out] outLLR the resulting LLR (whether converged or not)
 * @return 0 if failed to converge or (+/- the number of steps) if converged successfully. 
 *  
 */
int do_serialC_BP(qllr_t * outLLR, const mzd_t * const srow,
		   const csr_t * const H, const csr_t * const Ht,
		   const qllr_t LLR[], const params_t * const p){
  assert(outLLR != NULL);
  const int nchk = p->nchk; /** number of check nodes */
  const int nvar = p->nvar;     /** number of variable nodes */
  const int nz = p->mH->p[p->nchk]; /** number of messages */
  qllr_t *mesVtoC = malloc(nz*sizeof(qllr_t));
  qllr_t *mesCtoV = malloc(nz*sizeof(qllr_t));
  if((!mesVtoC)||(!mesCtoV)) ERROR("memory allocation failed");
  qllr_t * xLLR=NULL, *aLLR=NULL;
  int succ_BP=0;
  /** convenience setting: by default, use both LLR and average LLR */
  const int submode = ((p->submode & 3) == 0) ? (p->submode | 3) : p->submode;

  if((submode & 3) == 1) /** use xLLR only */
    xLLR = outLLR; /** return soft vector of LLR here */
  else{                /** use aLLR or both LLR and aLLR */
    aLLR = outLLR; 
    xLLR = calloc(nvar,sizeof(qllr_t)); /** needed for calculations in any case */
    if(!xLLR) ERROR("memory allocation failed");
  }

  mzp_t *pivots = mzp_init(nchk); 
  mzp_t *perm   = perm_p(NULL, pivots,0); /* actual permutation */

  for (int itry=0; itry<p->bpretry; itry++){

    //  if(p->debug &1)    printf("randomizing initial node order\n");
    pivots = mzp_rand(pivots); /* LAPAC-style random node permutation */
    perm   = perm_p(perm, pivots,0); /* actual permutation */
    
    /** init V->C messages to bare LLR */
    bp_init_VC(mesVtoC,H,LLR);

    /** one round of parallel to init C->V messages */
    for (int iv=0; iv < nvar; iv++)
      bp_do_CVv(iv, mesVtoC, mesCtoV, H, Ht, srow);

    for (int istep=1; istep <= p->steps  ; istep++){ /** main decoding cycle */
      if(p->submode & 16){
	pivots = mzp_rand(pivots); /* LAPAC-style random node permutation */
	perm   = perm_p(perm, pivots,0); /* actual permutation */    
      }
      
      for(int ii=0; ii<nchk; ii++){      
	int ic = perm->values[ii]; /** use the random permutation */
	int sbit = mzd_read_bit(srow,0,ic); /** syndrome bit at `ic` */
	//! compute  all V->C messages into `ic`;
	bp_do_VCc(ic, mesVtoC, mesCtoV, H, Ht, LLR);
	//! compute all C->V messages from `ic`;      
	bp_do_CVc(ic, mesVtoC, mesCtoV, H, Ht, sbit);
      }

      for(int iv = 0; iv< nvar; iv++) /** calculate expected LLR for variable nodes */
	xLLR[iv] = bp_do_LLR(iv,mesCtoV,Ht,LLR[iv]);
    
      if(submode&2) /** use average LLR */
	for(int iv = 0; iv< nvar; iv++)
	  aLLR[iv] = llr_from_dbl(p->bpalpha * dbl_from_llr(aLLR[iv]) +(1.0 - p->bpalpha) * dbl_from_llr(xLLR[iv])); 

#ifndef NDEBUG    
      if(p->debug & 8){
	if(submode&1) /** use regular LLR */
	  out_llr("x",p->nvar, xLLR);
	if(submode&2)
	  out_llr("a",p->nvar, aLLR);
      }
#endif
      /** do convergence check */
      if((submode&2) && syndrome_check(aLLR, srow,p->mH, p)){ 
	cnt_update(CONV_BP_AVG, istep);
	succ_BP=-istep;
	break;
      }
      else if((submode&1) && syndrome_check(xLLR,srow,p->mH, p)){
	cnt_update(CONV_BP, istep);
	succ_BP=istep;
	if(submode&2) /** need to copy */
	  for(int iv=0; iv < p->nvar; iv++)
	    outLLR[iv] = xLLR[iv];
	break;
      }
    }
    if (succ_BP)
      break; /* from the re-try loop */
  }
  /** clean up ********************/
  if(submode&2)
    free(xLLR);
  free(mesVtoC);
  free(mesCtoV);
  mzp_free(perm);
  mzp_free(pivots);

  /** default value is returned */

  return succ_BP;
}


/** @brief serial-V BP algorithm (sum-product or max-product depending on LLRs used).  
 * 
 * Same as `do_serialC_BP()` function above but using `variable-based` serial order: 
 * for each `iv` do:
 *    update C->V messages to `iv`.
 *    update V->C messages from `iv`;
 * 
 * @param[out] outLLR the resulting LLR (whether converged or not)
 * @return 0 if failed to converge or (+/- the number of steps) if converged successfully. 
 *  
 */
int do_serialV_BP(qllr_t * outLLR, const mzd_t * const srow,
		  const csr_t * const H, const csr_t * const Ht,
		  const qllr_t LLR[], const params_t * const p){
  assert(outLLR != NULL);
  //  const int nchk = p->nchk; /** number of check nodes */
  const int nvar = p->nvar;     /** number of variable nodes */
  const int nz = p->mH->p[p->nchk]; /** number of messages */
  qllr_t *mesVtoC = malloc(nz*sizeof(qllr_t));
  qllr_t *mesCtoV = malloc(nz*sizeof(qllr_t));
  if((!mesVtoC)||(!mesCtoV)) ERROR("memory allocation failed");
  qllr_t * xLLR=NULL, *aLLR=NULL;
  int succ_BP=0;
  /** convenience setting: by default, use both LLR and average LLR */
  const int submode = ((p->submode & 3) == 0) ? (p->submode | 3) : p->submode;
 
  if((submode & 3) == 1) /** use xLLR only */
    xLLR = outLLR; /** return soft vector of LLR here */
  else{                /** use aLLR or both LLR and aLLR */
    aLLR = outLLR; 
    xLLR = calloc(nvar,sizeof(qllr_t)); /** needed for calculations in any case */
    if(!xLLR) ERROR("memory allocation failed");
  }

  mzp_t *pivots = mzp_init(nvar); 
  mzp_t *perm   = perm_p(NULL, pivots,0); /* actual permutation */

  for (int itry=0; itry<p->bpretry; itry++){

    //  if(p->debug &1)    printf("randomizing initial node order\n");
    pivots = mzp_rand(pivots); /* LAPAC-style random node permutation */
    perm   = perm_p(perm, pivots,0); /* actual permutation */
    
    /** init V->C messages to bare LLR */
    bp_init_VC(mesVtoC,H,LLR);

    for (int istep=1; istep <= p->steps  ; istep++){ /** main decoding cycle */
      if(p->submode & 16){
	pivots = mzp_rand(pivots); /* LAPAC-style random node permutation */
	perm   = perm_p(perm, pivots,0); /* actual permutation */    
      }
      
      for(int ii=0; ii<nvar; ii++){      
	int iv = perm->values[ii]; /** use the random permutation */
	
	//! compute  all C->V messages into `iv`;
	bp_do_CVv(iv, mesVtoC, mesCtoV, H, Ht, srow);
	//! compute all V->C messages from `iv`;      
	bp_do_VCv(iv, mesVtoC, mesCtoV, H, Ht, LLR[iv]);
      }

      for(int iv = 0; iv< nvar; iv++) /** calculate expected LLR for variable nodes */
	xLLR[iv] = bp_do_LLR(iv,mesCtoV,Ht,LLR[iv]);
    
      if(submode&2) /** use average LLR */
	for(int iv = 0; iv< nvar; iv++)
	  aLLR[iv] = llr_from_dbl(p->bpalpha * dbl_from_llr(aLLR[iv]) +(1.0 - p->bpalpha) * dbl_from_llr(xLLR[iv])); 

#ifndef NDEBUG    
      if(p->debug & 8){
	if(submode&1) /** use regular LLR */
	  out_llr("x",p->nvar, xLLR);
	if(submode&2)
	  out_llr("a",p->nvar, aLLR);
      }
#endif
      /** do convergence check */
      if((submode&2) && syndrome_check(aLLR, srow,p->mH, p)){ 
	cnt_update(CONV_BP_AVG, istep);
	succ_BP=-istep;
	break;
      }
      else if((submode&1) && syndrome_check(xLLR,srow,p->mH, p)){
	cnt_update(CONV_BP, istep);
	succ_BP=istep;
	if(submode&2) /** need to copy */
	  for(int iv=0; iv < p->nvar; iv++)
	    outLLR[iv] = xLLR[iv];
	break;
      }
    }
    if (succ_BP)
      break; /* from the re-try loop TODO: fix retry counters!!!!! */
  }
  /** clean up ********************/
  if(submode&2)
    free(xLLR);
  free(mesVtoC);
  free(mesCtoV);
  mzp_free(perm);
  mzp_free(pivots);

  /** default value is returned */

  return succ_BP;
}



/** @brief recursive part of OSD */
int do_osd_recurs(const int minrow, const rci_t jstart, const int lev, qllr_t vE[], mzd_t *mE,
	      const csr_t * const sHt, const qllr_t LLR[],
	      const mzp_t * const pivs, const mzp_t * const skip_pivs, const params_t * const p){
  assert(lev<=p->lerr);
  assert(lev>0);

#ifndef NDEBUG  
  if(p->debug & 128){
    printf("starting OSD lev=%d of %d jstart=%d maxosd=%d\n",lev,p->lerr,jstart, p->maxosd);
    if((p->nvar <= 256)&&(p->debug &512))   mzd_print(mE);
    printf("E[%d]=%g  ",minrow,dbl_from_llr(vE[minrow]));
    if((p->nvar <= 256)&&(p->debug &512))   mzd_print_row(mE,minrow);
    for(int i=0; i<lev; i++){
      printf("E[%d]=%g  ",i,dbl_from_llr(vE[i]));
      if((p->nvar <= 256)&&(p->debug &512))   mzd_print_row(mE,i);
    }
    printf("\n");
  }
#endif   

  int last_lev = lev < p->lerr ? 0 : 1;
  int ich_here=0, ich_below=0;

  int knum = skip_pivs->length; /** number of non-pivot cols in `mH` to go over */
  if ((lev>1) && (knum > p->maxosd))
    knum = p->maxosd;

  for(rci_t j=jstart; j<knum; j++){ /** `j`th `non-pivot` entry */
    rci_t jj=skip_pivs->values[j]; /** actual `non-pivot` column we are looking at */
    mzd_copy_row(mE,lev,mE,lev-1); /** fresh copy of error vector to update */
    vE[lev]=vE[lev-1];
#ifdef NDEBUG
    if(mzd_read_bit(mE,lev,jj))
      ERROR("set bit non-pivot position lev=%d jstart=%d jj=%d\n",lev,jstart,jj);
#endif 	    
    vE[lev] += LLR[jj];
    mzd_flip_bit(mE,lev,jj);
    for(int ii=sHt->p[jj]; ii < sHt->p[jj+1] ; ii++){
      int ir=pivs->values[sHt->i[ii]]; /** pivot for `ii`-th non-zero elt in row `jj` */
      if(mzd_read_bit(mE,lev,ir))
        vE[lev] -= LLR[ir]; /** `1->0` flip */      
      else
	vE[lev] += LLR[ir]; /** `0->1` flip */
      mzd_flip_bit(mE,lev,ir);
    }
    if (vE[lev] < vE[minrow] - 1e-10){ /** update min-energy vector */
#ifndef NDEBUG
      if(p->debug & 128){/** inf set decoding */
	printf("lev=%d j=%d jj=%d E0=%g -> E=%g success\n", lev,j,jj,
	       dbl_from_llr(vE[minrow]),dbl_from_llr(vE[lev]));
	if((p->nvar <= 256)&&(p->debug &512))
	  mzd_print_row(mE,lev);
      }
#endif
      vE[minrow]=vE[lev];
      mzd_copy_row(mE,minrow, mE,lev);
      ich_here++;
    }

    if((!last_lev) && (vE[minrow] >= 2*p->LLRmin)){ /** go up one recursion level */
      if(j+1<knum){
        ich_below += do_osd_recurs(minrow,j+1,lev+1,vE,mE, sHt, LLR,pivs,skip_pivs, p);
      }
    }
  }
#ifndef NDEBUG
  if(p->debug & 128)
    if(ich_here + ich_below)
      printf("exiting lev=%d of recursion, here ch=%d below ch=%d\n",
             lev,ich_here, ich_below);
#endif 
  return 0;
}

/** @brief run OSD up to `p->lerr` using up to `p->maxosd` columns at high levels.
    @param[in] LLR input LLR values (e.g., as returned by PB).
    @param[out] LLR output LLR values (+/- 1) from OSD. (`???`)
    @param srow syndrome row to match or NULL `to compute codewords` (???)
    @param Check matrix to use
    @output always `1` (???) 
*/
int do_osd_start(qllr_t * LLR, const mzd_t * const srow,
		 const csr_t * const H, const params_t * const p){
  assert(p->lerr >= 0); /** use `-1` for no OSD */

  /** generate permutation, create binary matrix */
  const int nvar = H->cols;
  mzp_t * perm = mzp_init(nvar);    /* initialize the permutation */
  mzp_t * pivs = mzp_init(nvar); /** list of pivot columns */
  if((!pivs) || (!perm))
    ERROR("memory allocation failed!\n");
  mzd_t * mH = mzd_from_csr(NULL, H);
  mzd_t * mSrow = mzd_copy(NULL,srow);
  const int minrow = p->lerr +1; /** WARNING: the minimum energy vector at this row */
  mzd_t * mE = mzd_init(minrow +1, nvar); /* storage for error vectors */
  //  mzd_print(mSrow);
  perm = sort_by_llr(perm, LLR, p);   /** order of decreasing `p` */

  /** full row echelon form (gauss elimination) using the order of `p`,
   * on the block matrix `[H|S]` (matrix / vector).
   */
  int rank=0;
  for(int i=0; i< p->nvar; i++){
    int col=perm->values[i];
    int ret=matvec_gauss_one(mH, mSrow, col, rank);
    if(ret)
      pivs->values[rank++]=col;
  }
  /** create error vector and calculate energy (`OSD0`) */
  qllr_t * vE = calloc(minrow+1,sizeof(qllr_t));
  if(vE==NULL) ERROR("memory allocation failed");
  for(int i=0;i< rank; i++)
    if(mzd_read_bit(mSrow,0,i)){
      mzd_write_bit(mE,minrow,pivs->values[i],1);
      vE[minrow] += p->vLLR[pivs->values[i]];
    }

#ifndef NDEBUG
  if(p->debug & 128){/** inf set decoding */
    printf("lev=0 E0=%g \n", dbl_from_llr(vE[minrow]));
    if((p->nvar <= 256)&&(p->debug &512))
      mzd_print_row(mE,minrow);
  }
#endif
  
  
  if((p->lerr>0) && (vE[minrow] >= 2*p->LLRmin)){
    /** TODO: `(later, maybe)` do Gauss while updating `s` if non-NULL */
  
    /** prepare CSR version of modified `Ht` */
    mzd_t * mHt = mzd_transpose(NULL, mH);
    csr_t * sHt = csr_from_mzd(NULL, mHt);
    mzd_free(mHt);
    
    mzp_t * skip_pivs = do_skip_pivs(rank, pivs);    
    vE[0]=vE[minrow];
    mzd_copy_row(mE,0, mE,minrow); /** initial row for OSD */

    /** with current `Emin`,`vmin` and `E`, `v`, launch recursion */
    do_osd_recurs(minrow, 0, 1, vE, mE, sHt, p->vLLR, pivs, skip_pivs, p);
    mzp_free(skip_pivs);
    csr_free(sHt);
  }

  /** copy bits to LLR vector WARNING: `this is a hack!` */
  for (int i=0; i<nvar; i++)
    LLR[i] = (mzd_read_bit(mE,minrow,i)) ? -5000 : 5000;
  if((p->debug & 128) && (nvar <= 256))
    out_llr(":",nvar,LLR);
  /** clean up */
  free(vE);
  mzd_free(mE);
  mzd_free(mSrow);
  mzd_free(mH);
  mzp_free(pivs);
  mzp_free(perm);
  return 1;
}
