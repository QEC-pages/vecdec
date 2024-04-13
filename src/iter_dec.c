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
    printf("# FAIL_FRAC TOTAL CONV_TRIVIAL CONV_BP CONV_BP_AVG SUCC_TRIVIAL SUCC_BP SUCC_TOT\n");
  printf(" %g %lld %lld %lld %lld %lld %lld %lld\n",
	 (double ) (cnt[TOTAL]-cnt[SUCC_TOT])/cnt[TOTAL],
	 cnt[TOTAL], cnt[CONV_TRIVIAL],cnt[CONV_BP], cnt[CONV_BP_AVG],
	 cnt[SUCC_TRIVIAL], cnt[SUCC_BP], cnt[SUCC_TOT]);
}

void cnt_update(extr_t which, int iteration){
  const long long int max = (LLONG_MAX >> 3);
  cnt[which]++;
  iter1[which]+=iteration;
  iter2[which]+=iteration*iteration;
  if((cnt[which] > max) ||(iter1[which] > max) || (iter2[which]>max)){
    cnt_out(1);
    ERROR("too many iterations\n"); 
  }
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
		   [[maybe_unused]] const params_t * const p){
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

/** @brief baseline parallel BP algorithm (sum-product or max-product depending on LLRs used).  
 * 
 * Decode `e` vector using syndrome bits in `srow = e*Ht` and apriori
 * bit LLR values in `LLR`.  Assume `srow` non-zero (no check for a
 * trivial solution).
 * 
 * @param[out] outLLR the resulting LLR (whether converged or not)
 * @return 0 if failed to converge or (the number of steps) if converged successfully
 */
int do_parallel_BP(qllr_t * outLLR, const mzd_t * const srow,
		   const csr_t * const H, const csr_t * const Ht,
		   const qllr_t LLR[], const params_t * const p){
  assert(outLLR != NULL);
  const int nchk = p->nchk; /** number of check nodes */
  const int nvar = p->nvar;     /** number of variable nodes */
  const int nz = p->mH->p[p->nchk]; /** number of messages */
  qllr_t *mesVtoC = malloc(nz*sizeof(qllr_t));
  qllr_t *mesCtoV = malloc(nz*sizeof(qllr_t));
  //  double * xLLR = malloc(nvar*sizeof(double));
  qllr_t * const xLLR = outLLR; /** return soft vector of LLR here */
  qllr_t * aLLR = calloc(nvar,sizeof(qllr_t)); /** averaged LLR */
  if((!aLLR)||(!mesVtoC)||(!mesCtoV))
    ERROR("memory allocation failed");
  int succ_BP=0;
  
    /** init V->C messages to bare LLR */
  for(int ic=0; ic<nchk; ic++){ /** target `check` node index */
    for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
      const int iv = H->i[j]; /** `variable` node index */
      mesVtoC[j]= LLR[iv];   /** `j` is the edge (message) index */
    }
  }
  
  for (int istep=1; istep <= p->steps  ; istep++){ /** main decoding cycle */
    //    cnt[TOTAL]++;
    /** C -> V messages */
    for(int iv=0; iv<nvar; iv++){ /** target `variable` node index */
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
	      msg = mesVtoC[j1];
	      new=0;
	    }
	    else
	      msg = boxplus(msg, mesVtoC[j1]);
	  }
	}
	mesCtoV[j] = (sbit ? -msg : msg);
      }
    }

        
    for(int iv = 0; iv< nvar; iv++){
      qllr_t val= LLR[iv]; /** HERE! (WARNING: do we need this?) */
      for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++){
	// const int ic1 = Ht->i[j1];/** aux `check` node index */
	val += mesCtoV[j1];
      }
      xLLR[iv] = val;
      /** TODO: play with the decay value `0.5` */
      if(!(p->submode&1))
	aLLR[iv] = llr_from_dbl(0.5 * dbl_from_llr(aLLR[iv]) + dbl_from_llr(val)); 
    }
    if(p->debug & 8){
      out_llr("x",p->nvar, xLLR);
      if(!(p->submode&1))
	out_llr("a",p->nvar, aLLR);
      
      //      for(int iv = 0; iv < max; iv++)
      //	printf(" %5.1g%s", dbl_from_llr(aLLR[iv]), iv< max ? "" : "\n");
    }


    if(syndrome_check(xLLR,srow,p->mH, p)){
      //      outLLR = xLLR; 
      cnt_update(CONV_BP, istep);
      succ_BP=istep;
      break;
    }
    else if((!(p->submode&1)) && (syndrome_check(aLLR,srow,p->mH, p))){
      for(int iv=0; iv< p->nvar; iv++)
	outLLR[iv] = aLLR[iv];
      cnt_update(CONV_BP_AVG, istep);
      succ_BP=-istep;
      break;
    }

    
    /* V -> C messages */ 
    for(int ic=0; ic<nchk; ic++){ /** target `check` node */
      for(int j = H->p[ic]; j < H->p[ic+1] ; j++){
	const int iv = H->i[j]; /** origin `variable` node index */
#if 0 /** TODO: verify these give identical results! */
	qllr_t msg = LLR[iv];
	for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++){
	  const int ic1 = Ht->i[j1];/** aux `check` node index */
	  if(ic1 != ic)
	    msg += mesCtoV[j1];
	}
	mesVtoC[j] = msg;
#else /* simplified calculation */
	//	assert(ic == Ht->i[j]);
	for(int j1 = Ht->p[iv]; j1 < Ht->p[iv+1] ; j1++){
	  const int ic1 = Ht->i[j1];/** aux `check` node index */
	  if(ic1 == ic)
	    mesVtoC[j] = xLLR[iv] - mesCtoV[j1];
	  /** TODO: use binary search here instead */
	}
#endif /* 0 */	
      }
    }

  }
  if(!succ_BP)
    outLLR = xLLR; /** TODO: make an option to change this to aLLR */
  /** clean up ********************/
  //  free(xLLR);
  if(aLLR){
    free(aLLR);
    aLLR=NULL;
  }
  free(mesVtoC);
  free(mesCtoV);

  return succ_BP;
}


/** @brief recursive part of OSD */
int do_osd_recurs(int minrow, rci_t jstart, int lev, qllr_t vE[], mzd_t *mE,
	      const csr_t * const sHt, const qllr_t LLR[],
	      const mzp_t * const skip_pivs, const params_t * const p){
  assert(lev<=p->lerr);
  assert(lev>0);
  int last_lev = lev < p->lerr ? 0 : 1;
  int ich_here=0, ich_below=0;

  int knum = skip_pivs->length; /** number of non-pivot cols in `mH` to go over */
  if ((lev>1) && (knum > p->maxosd))
    knum = p->maxosd;

  for(rci_t j=jstart; j<knum; j++){ /** `j`th `non-pivot` entry */
    rci_t jj=skip_pivs->values[j]; /** actual `non-pivot` column we are looking at */
    mzd_copy_row(mE,lev,mE,lev-1); /** fresh copy of error vector to update */
    vE[lev]=vE[lev-1];
#ifndef NDEBUG      
    if(mzd_read_bit(mE,lev,jj)) /** sanity check */
      ERROR("set bit found at lev=%d jj=%d\n",lev,jj);
#endif
    vE[lev] += LLR[jj];

    for(int ii=sHt->p[jj]; ii < sHt->p[jj+1] ; ii++){
      int ir=sHt->i[ii]; /** pos of `ii`-th non-zero elt in row `jj` */
      if(mzd_read_bit(mE,lev,ir)){
	mzd_flip_bit(mE,lev,ir);
        vE[lev] -= LLR[ir]; /** `1->0` flip */
      }
      else
	vE[lev] += LLR[ir]; /** `0->1` flip */
    }
    if (vE[lev] < vE[minrow] - 1e-10){ /** update min-energy vector */
#ifndef NDEBUG
      if(p->debug & 128){/** inf set decoding */
	printf("lev=%d j=%d jj=%d E0=%g E=%g success\n", lev,j,jj,
	       dbl_from_llr(vE[minrow]),dbl_from_llr(vE[lev]));
      }
#endif
      vE[minrow]=vE[lev];
      mzd_copy_row(mE,minrow, mE,lev);
      ich_here++;
    }
    if(!last_lev){ /** go up one recursion level */
      if(j+1<knum){
        ich_below+=do_osd_recurs(minrow,j+1,lev+1,vE,mE, sHt, LLR,skip_pivs, p);
      }
    }
  }

  if(p->debug & 128)
    if(ich_here + ich_below)
      printf("exiting lev=%d of recursion, here ch=%d below ch=%d\n",
             lev,ich_here, ich_below);

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
  mzd_print(mSrow);
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
  for(int i=0;i< rank; i++)
    if(mzd_read_bit(mSrow,0,i)){
      mzd_write_bit(mE,minrow,pivs->values[i],1);
      vE[minrow] += LLR[pivs->values[i]];
    }
  
  if(p->lerr>0){
    /** TODO: `(later, maybe)` do Gauss while updating `s` if non-NULL */
  
    /** prepare CSR version of modified `Ht` */
    mzd_t * mHt = mzd_transpose(NULL, mH);
    mzd_free(mH);
    csr_t * sHt = csr_from_mzd(NULL, mHt);
    mzd_free(mHt);
    
    mzp_t * skip_pivs = do_skip_pivs(rank, pivs);    
    vE[0]=vE[minrow];
    mzd_copy_row(mE,0, mE,minrow); /** initial row for OSD */

    /** with current `Emin`,`vmin` and `E`, `v`, launch recursion */
    do_osd_recurs(minrow, 0, 1, vE, mE, sHt, LLR, skip_pivs, p);
  }

  /** copy bits to LLR vector WARNING: `this is a hack!` */
  for (int i=0; i<nvar; i++)
    LLR[i] = (mzd_read_bit(mE,minrow,i)) ? -1 : 1;
  
  /** clean up */
  mzd_free(mE);
  mzd_free(mSrow);
  mzd_free(mH);
  mzp_free(pivs);
  mzp_free(perm);
  return 1;
}
