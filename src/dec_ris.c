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
  mzd_set_ui(mE,0); /** zero matrix to store `errors by column` */

  /** first pass ******************************************* */
  int rank=0;
  if (p->steps > 0){
    perm = sort_by_llr(perm, p->vLLR, p);   /** order of decreasing `p` */
    /** full row echelon form (gauss elimination) using the order of `p`,
     * on the block matrix `[H|S]` (in fact, two matrices).
     */
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
    for(int i=0;i< rank; i++)
      mzd_copy_row(mE,pivs->values[i],mS,i);
  }
  mzd_t *mEt0 = mzd_transpose(NULL,mE);
  for(int i=0; i< mS->ncols; i++)
    vE[i]=mzd_row_energ(p->vLLR,mEt0,i);

  if(p->debug & 512){
    printf("mEt0 after round 0:\n");
    mzd_print(mEt0);
  }

  if((p->steps > 0) &&(p->lerr > 0)){  /** do information-set decoding `**********************` */
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
