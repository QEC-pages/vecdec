/**
 *  @file vecdec.h
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
#include "utils.h"
#include "util_m4ri.h"
#include "vecdec.h"

params_t prm={ .nrows=0, .n=0, .ncws=0, .steps=1, .lerr=0, .swait=0,
  .nvec=16, .ntot=1, .nfail=0, .seed=0, 
  .debug=1, .fin=NULL,
  .colw=10, .mode=0, .maxJ=20, 
  .pmin=-1, .pmax=-1, .pstep=1e-3,
  .vP=NULL, .vLLR=NULL, .mH=NULL, .mHt=NULL, .mL=NULL, .mLt=NULL};

/** calculate the energy of the row `i` */
double mzd_row_energ(double *coeff, const mzd_t *A, const int i){
  double ans=0;
  for(rci_t j = 0; j < A->ncols; ++j)
    if(mzd_read_bit(A, i, j))
      ans += coeff[j];
  return (ans);
  /** todo: rewrite in terms of `__builtin_clzll` */
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

/** @brief create a sample of errors to play with.
 *  @param mHe matrix with `maxJ` columns to return the syndrome `H*e`
 *  @param mLe matrix with `maxJ` columns for logical error `L*e`
 *  @param p params_t structure with error model information
 * @return max number of generated errors of any one kind 
 *  todo: perhaps allow for an independent error model (matrices `mHt`, `mLt` and vector `vP`) ? 
 */
int do_errors(mzd_t *mHe, mzd_t *mLe,
               const params_t * const p){
  
  assert((mHe!=NULL) && (mLe!=NULL)); /** sanity check */
  
  int max = p->maxJ;  /** current size of `vec` */
  int * vec = malloc(p->maxJ * sizeof(int));
  
  mzd_set_ui(mHe,0); /** zero matrix */
  mzd_set_ui(mLe,0); /** zero matrix */
  int nvec = mHe->ncols; 
  /** for each error type (column of `mH` = row of `mHt`) */
  for(int i=0; i < p->mHt->rows; i++){
    int ivec=0;
    /** prepare the list of syndrome columns to deal with */
    double onebyL = -1.0/log(1.0 - p->vP[i]);
    int j =(int )floor(onebyL * rnd_exponential());
#ifndef NDEBUG    
    if(p->debug & 256) /** details of error setup */
      printf("p[%d]=%g onebyL=%g j=%d nvec=%d\n",i,p->vP[i],onebyL,j,nvec);
#endif /* NDEBUG */
    if(j < nvec){/** otherwise we are to skip this error altogether */
      do{
        if(ivec >= max){
          max=2*max;
          vec=realloc(vec,max * sizeof(int));
        }
        vec[ivec++]=j;
        j += (int )ceil(onebyL * rnd_exponential());
#ifndef NDEBUG            
        if(p->debug & 256) /** details of error setup */
          printf("p[%d]=%g onebyL=%g j=%d nvec=%d\n",i,p->vP[i],onebyL,j,nvec);
#endif /* NDEBUG */        
      }
      while(j < nvec);
      
#ifndef NDEBUG    
    if(p->debug & 256) /** details of error setup */
      for (j=0; j<ivec; j++)
        printf("vec[%d]=%d%s",j,vec[j],j+1==ivec?" ###\n":" ");
#endif /* NDEBUG */
    
      /** flip the bits in `mHe` row by row to speed it up */
      for(int ir = p->mHt->p[i]; ir < p->mHt->p[i+1]; ir++){
        int irow = p->mHt->i[ir]; /** column index in `mH` = row in `mHt` */
        for(j=0; j < ivec; j++)
          mzd_flip_bit(mHe,irow,vec[j]);
      }
#ifndef NDEBUG      
      if(p->debug & 256){ /** details of error setup */
        printf("i=%d j=%d\n",i,j);
        printf("matrix mHe:\n");  mzd_print(mHe);
      }
#endif /* NDEBUG */
      
      /** flip the bits in `mLe` row by row */
      for(int ir = p->mLt->p[i]; ir < p->mLt->p[i+1]; ir++){
        int irow = p->mLt->i[ir]; /** column index in `mL` */
        for(j=0; j<ivec; j++)
          mzd_flip_bit(mLe,irow,vec[j]);
      }
    }
  }
  //  p->maxJ = max; /** to reduce the number of `realloc`s next time */
  free(vec);
  return max;
}

/** @brief helper function to sort `ippair_t` 
 *  use `qsort(array, len, sizeof(ippair_t), cmp_ippairs);`
 */
static inline int cmp_ippairs(const void *a, const void *b){
  const double pa=((ippair_t *) a) -> prob;
  const double pb=((ippair_t *) b) -> prob;
  if (pa<pb)
    return +1;
  else if (pa>pb)
    return -1;
  return 0;
}

/** @brief helper function to sort `int` 
 *  use `qsort(array, len, sizeof(ippair_t), cmp_int);`
 */
static inline int cmp_rci_t(const void *a, const void *b){
  const rci_t va= *((rci_t *) a);
  const rci_t vb= *((rci_t *) b);
  return va-vb;
#if 0  
  if (va<vb)
    return +1;
  else if (va>vb)
    return -1;
  return 0;
#endif
}



/** @brief return permutation = decreasing probabilities */
mzp_t * sort_by_prob(mzp_t *perm, params_t const * const p){
  /** prepare array of ippairs */
  ippair_t * pairs = malloc(p->n * sizeof(ippair_t));
  if (!pairs)
    ERROR("memory allocation failed\n");
  for(int i=0; i<p->n; i++){
    pairs[i].index = i;
    pairs[i].prob = p->vP[i];
  }
  qsort(pairs, p->n, sizeof(ippair_t), cmp_ippairs);
  for(int i=0; i<p->n; i++)
    perm->values[i] = pairs[i].index;

  //  for(int i=0; i<p->n; i++) printf("i=%d p[perm[i]]=%g\n",i,p->vP[perm->values[i]]);
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

/**
 * @brief do local search up to `p->lerr` inclusive recursively 
 * @param vE0 best energy values for each `e` a row in `mE0`
 * @param mE0 best error vectors so far for each syndrome (rows)
 * @param jstart start local search from this column in `mH`
 * @param lev recusion level, must not exceed `p->lerr`
 * @param mE input error vectors (rows)
 * @param mH check matrix in row echelon form (pivots in `pivs`)
 * @param skip_pivs `sorted` list of `n-rank` non-pivot positions in `mH`
 * @param pivs list of `rank` pivots returned by gauss (length = `n`)
 * @param p pointer to the remaining program parameters
 * @return number of updated error vectors 
 * @todo: see if binary ops + transposition is faster 
 */
int do_local_search(double *vE0, mzd_t * mE0, const rci_t jstart, const int lev,
                     const mzd_t * const mE, const mzd_t * const mH,
                     const mzp_t * const skip_pivs, const mzp_t * const pivs,
                     const params_t * const p){  
  assert(lev<=p->lerr);
  if(p->debug&128)
    printf("entering lev=%d of recursion jstart=%d\n",lev,jstart);
  int ich_here=0, ich_below=0;
  rci_t knum = skip_pivs->length; /** number of non-pivot cols in `mH` to go over */
  rci_t rank = p->n - knum; /** number of valid pivot cols */
  rci_t rnum; /** number of non-zero entries in rlis (for each `j`) */
  int * rlis = malloc(rank * sizeof(int)); 
  if(!rlis) ERROR("memory allocation failed!");
  mzd_t *mE1 = NULL; //mzd_copy(NULL,mE); /** error vectors to update */
  
  for(rci_t j=jstart; j<knum; j++){ /** `j`th `non-pivot` entry */
    mE1 = mzd_copy(mE1,mE); /** fresh copy of error vectors to update */

    rci_t jj=skip_pivs->values[j]; /** actual `non-pivot` column we are looking at */
    rnum=0; /** prepare list of positions to update for `jj`th col of `mH` */
    for(rci_t ir=0; ir<rank; ir++)
      if(mzd_read_bit(mH,ir,jj)) /** non-zero bit */
        rlis[rnum++] = pivs->values[ir]; /** column of `mH` to update */
    if (p->debug & 128){
      printf("jj=%d rlis: ",jj);
      for(int ir=0; ir< rnum; ir++)
        printf(" %d%s",rlis[ir],ir+1==rnum?"\n":"");
    }
    
    for(rci_t is=0; is < mE1->nrows; is++){ /** syndrome rows */
      if(mzd_read_bit(mE1,is,jj)) /** sanity check */
        ERROR("bit found at is=%d jj=%d\n",is,jj);
      
      if(vE0[is]>1e-15){ /** substitute for zero syndrome check */
        double energ = mzd_row_energ(p->vLLR,mE1,is);
        mzd_flip_bit(mE1,is,jj);
        energ += p->vLLR[jj];
        for(rci_t ir = 0 ;  ir < rnum; ++ir){
          const int ii = rlis[ir]; /** position to update */
          if(mzd_read_bit(mE1,is,ii))
            energ -= p->vLLR[ii]; /** `1->0` flip */
          else
            energ += p->vLLR[ii]; /** `0->1` flip */
          mzd_flip_bit(mE1,is,ii);
        }        
        if(energ < vE0[is]-1e-10){
          if(p->debug & 128){/** inf set decoding */
            printf("lev=%d j=%d jj=%d is=%d E0=%g E=%g success:\n",
                   lev,j,jj,is,vE0[is],energ);
            assert(fabs(energ - mzd_row_energ(p->vLLR,mE1,is))<1e-8);
            mzd_print_row(mE0,is);
            mzd_print_row(mE1,is);            
          }
          vE0[is]=energ;
          mzd_copy_row(mE0,is, mE1,is);
          ich_here++;
        }
        else{
          if(p->debug & 128){/** inf set decoding */
            printf("lev=%d j=%d jj=%d is=%d E0=%g E=%g no change:\n",
                   lev,j,jj,is,vE0[is],energ);
            assert(fabs(energ - mzd_row_energ(p->vLLR,mE1,is))<1e-8);
            mzd_print_row(mE0,is);
            mzd_print_row(mE1,is);            
          }
        }
      }
    }
    
    if(lev+1 < p->lerr){ /** go up one recursion level */
      if(j+1<knum){
        ich_below+=do_local_search(vE0,mE0,j+1,lev+1,mE1, mH, skip_pivs, pivs,p);
      }
    }
  }
  if(p->debug & 128)
    if(ich_here + ich_below)
      printf("exiting lev=%d of recursion, here ch=%d below ch=%d\n",
             lev,ich_here, ich_below);
  free(rlis);
  mzd_free(mE1);
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
  double *vE = calloc(mS->ncols,sizeof(double)); /**< best energies */
  if((!mE) || (!vE))
    ERROR("memory allocation failed!\n");

  mzp_t * perm=mzp_init(p->n); /** identity column permutation */
  mzp_t * pivs=mzp_init(p->n); /** list of pivot columns */
  if((!pivs) || (!perm))
    ERROR("memory allocation failed!\n");
  
  /** first pass ******************************************* */
  perm = sort_by_prob(perm, p);   /** order of decreasing `p` */
  /** full row echelon form (gauss elimination) using the order of `p`,
   * on the block matrix `[H|S]` (in fact, two matrices).
   */  
  int rank=0;
  for(int i=0; i< p->n; i++){
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

  if(p->lerr){  /** do information-set decoding `**********************` */
    mzp_t * skip_pivs = do_skip_pivs(rank, pivs);
    mzd_t * mEt1 = mzd_copy(NULL, mEt0);    
    do_local_search(vE, mEt0, 0, 0, mEt1, mH, skip_pivs, pivs, p);
    if(p->debug & 512){
      printf("mEt0 after local search:\n");
      mzd_print(mEt0);
    }
    mzd_free(mEt1);
    free(skip_pivs);
  }
  
  int iwait=0, ichanged=0;
  /** main loop over permutations * `***************************` */
  for (int ii=1; ii< p->steps; ii++){
    pivs=mzp_rand(pivs); /** random pivots LAPAC-style */
    mzp_set_ui(perm,1); perm=perm_p_trans(perm,pivs,0); /**< corresponding permutation */
    /** todo: make it probability-dependent */
    rank=0;
    for(int i=0; i< p->n; i++){
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
    for(int i=0; i< mEt->nrows; i++){
      double energ=mzd_row_energ(p->vLLR,mEt,i);
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
    if(p->lerr){  /** do information-set decoding `**********************` */
      mzp_t * skip_pivs = do_skip_pivs(rank, pivs);
      do_local_search(vE, mEt0, 0, 0, mEt, mH, skip_pivs, pivs, p);
      if(p->debug & 512){
        printf("mEt0 after local search:\n");
        mzd_print(mEt0);
      }
      free(skip_pivs);
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

/** @brief read detector error model (DEM) created by `stim`.
 * Immediately create CSR matrices `p->mL` and `p->mH` and vector `p->vP`.
 * @param fnam file name for reading DEM from 
 * @param p structure to store produced matrices 
 */ 
void read_dem_file(char *fnam, params_t * const p){
  ssize_t linelen, col=0;
  size_t lineno=0, bufsiz=0; /**< buffer size for `readline` */
  char *buf = NULL;          /** actual buffer for `readline` */
  p->nzH = p->nzL = 0;  /** count non-zero entries in `H` and `L` */
  int maxH=100, maxL=100, maxN=100;
  double *inP = malloc(maxN*sizeof(double));
  int_pair * inH = malloc(maxH*sizeof(int_pair));
  int_pair * inL = malloc(maxL*sizeof(int_pair));
  if ((!inP)||(!inH)||(!inL))
    ERROR("memory allocation failed\n");

  if(p->debug & 1)
    printf("# opening DEM file %s\n",fnam);  
  FILE *f = fopen(fnam, "r");
  if(f==NULL)
    ERROR("can't open the file %s for reading\n",fnam);

  int r=-1, k=-1, n=0;
  int iD=0, iL=0; /** numbers of `D` and `L` entries */
  do{ /** read lines one-by-one until end of file is found *************/
    lineno++; col=0; linelen = getline(&buf, &bufsiz, f);
    if(linelen<0)
      break;
    //    printf("# %s",buf);
    char *c=buf;
    double prob;
    int num=0,val;
    while(isspace(*c)){ c++; col++; } /** `skip` white space */
    if((*c != '\0')&& (*c != '#') &&(col < linelen)){
      if(sscanf(c,"error( %lg ) %n",&prob,&num)){
        if((prob<=0)||(prob>=1))
          ERROR("probability should be in (0,1) exclusive p=%g\n"
                "%s:%zu:%zu: '%s'\n", prob,fnam,lineno,col+1,buf);              
        c+=num; col+=num;
        if(n>=maxN){
          maxN=2*maxN;
          inP=realloc(inP,maxN*sizeof(*inP));
        }
        inP[n]=prob;
        do{/** deal with the rest of the line */
          num=0;
          if(sscanf(c," D%d %n",&val, &num)){/** `D` entry */
            c+=num; col+=num;
            assert(val>=0);
            if(val>r)
              r=val+1;  /** update the number of `D` pairs */
            if(iD>=maxH){
              maxH=2*maxH;
              inH=realloc(inH,maxH*sizeof(*inH));
            }
            inH[iD].a=val;   /** add a pair */
            inH[iD++].b=n;
          }
          else if(sscanf(c," L%d %n",&val, &num)){/** `L` entry */
            c+=num; col+=num;
            assert(val>=0);
            if(val>k)
              k=val+1;  /** update the number of `L` pairs */
            if(iL>=maxL){
              maxL=2*maxL;
              inL=realloc(inL,maxL*sizeof(*inL));
            }
            inL[iL].a=val;   /** add a pair */
            inL[iL++].b=n;
          }
          else
            ERROR("unrecognized entry %s"
                "%s:%zu:%zu: '%s'\n",c,fnam,lineno,col+1,buf);              
        }
        while((c[0]!='#')&&(c[0]!='\n')&&(c[0]!='\0')&&(col<linelen));
        n++;
      }
      else if (sscanf(c,"detector( %d %n",&val,&num)){
        /** do nothing */
        //        printf("# ignoring row[%zu]=%s\n",lineno,c);
      }
      else if (sscanf(c,"shift_detectors( %d %n",&val,&num)){
        /** do nothing */
        //        printf("# ignoring row[%zu]=%s\n",lineno,c);
      }
      else
        ERROR("unrecognized DEM entry %s"
              "%s:%zu:%zu: '%s'\n",c,fnam,lineno,col+1,buf);   

    }
    /** otherwise just go to next row */
  }
  while(!feof(f));
  p->nrows = r;
  p->ncws = k;
  p->n = n;
  if(p->debug &1)
    printf("# read DEM: r=%d k=%d n=%d\n",r,k,n);

  p->mH = csr_from_pairs(p->mH, iD, inH, r, n);  
  p->mL = csr_from_pairs(p->mL, iL, inL, k, n);

  p->mHt = csr_transpose(p->mHt, p->mH);    
  p->mLt = csr_transpose(p->mLt,p->mL);
  /** todo: fix reallocation logic to be able to reuse the pointers model */
  if(p->vP)
    free(p->vP);
  p->vP=inP;
  p->vLLR = malloc(n*sizeof(double));
  assert(p->vLLR !=0);
  for(int i=0;  i < n; i++)
    p->vLLR[i] = p->vP[i] > MINPROB ? log((1.0/p->vP[i] -1.0)) : log(1/MINPROB - 1);

  if(p->debug & 2)/** `print` out the entire error model ******************** */
    printf("# error model read: r=%d k=%d n=%d\n",p->nrows, p->ncws, p->n);
    
  if (buf) 
    free(buf);
  free(inH);
  free(inL);
}

/** @brief given the error model read, prepare for decoding
    Create vectors `p->vP`, `p->LLR` and sparse matrices `p->mH` and `p->mL`
    @param in the error model array in
    @param p contains error model parameters and place to store vecs and matrices
    @param prob if positive, alternative global probability to use 
    @output nothing (modified data in `p`)
*/
void mat_init(one_prob_t **in, params_t *p){
  //  int init_mat = (p->vP == NULL ? 1 : 0 );  
  p->vP = malloc(p->n * sizeof(double));
  p->vLLR = malloc(p->n * sizeof(double));
  p->mH = csr_init(NULL, p->nrows, p->n, p->nzH);
  p->mL = csr_init(NULL, p->ncws,  p->n, p->nzL); /** transposed */
  if((!p->vP) || (!p->vLLR) || (!p->mH) || (!p->mL))
    ERROR("memory allocation failed!\n");  
  int ipair1=0, ipair2=0;
  for (int i=0; i< p->n; i++){
    one_prob_t *row = in[i];

    double pp = row->p;
    p->vP[i] = pp;
    p->vLLR[i] = pp > MINPROB ? log((1.0/pp -1.0)) : log(1/MINPROB - 1);
    
    int j=0;
    for( ; j< row->n1; j++){
      p->mH->i[ipair1]   = i;           /** column */
      p->mH->p[ipair1++] = row->idx[j]; /** row */
    }
    for( ; j< row->n2; j++){
      p->mL->i[ipair2]   = i;            /** column */
      p->mL->p[ipair2++] = row->idx[j]; /** row */
    }
  };
  p->mH->nz  = p->nzH;
  csr_compress(p->mH);
  
  p->mHt = csr_transpose(p->mHt, p->mH);
  
  p->mL->nz = p->nzL;
  csr_compress(p->mL);
  
  p->mLt = csr_transpose(p->mLt,p->mL);

#ifndef NDEBUG  
  if((p->debug & 2)&&(p->debug &512)){ /** print resulting vectors and matrices */    
    for(int i=0; i< p->n; i++)
      printf("%g%s",p->vP[i],i+1<p->n?" ":"\n");
    for(int i=0; i< p->n; i++)
      printf("%g%s",p->vLLR[i],i+1<p->n?" ":"\n");
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

void prob_init(params_t *p, double prob){
  double pp = prob;
  double LLR = pp > MINPROB ? log((1.0/pp -1.0)) : log(1/MINPROB - 1);
  for(int i=0; i< p->n; i++){
    p->vP[i] = pp;
    p->vLLR[i] = LLR;
  }
}

int var_init(int argc, char **argv, params_t *p){

  int dbg=0;
  double val;
  
  if(argc<=1)
    ERROR("try \"%s -h\" for help",argv[0]);

  for(int i=1; i<argc; i++){  /** `debug` */
    if(sscanf(argv[i],"debug=%d",& dbg)==1){
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
    else if(sscanf(argv[i],"mode=%d",& dbg)==1){
      if(dbg==0)
	p->mode = 0;
      else{
	p->mode ^= dbg;
        if(p->debug&1)
          printf("# read %s, mode=%d octal=%o\n",argv[i],p->mode,p->mode);
      }
    }
    else if (sscanf(argv[i],"nvec=%d",&dbg)==1){ /** `nvec` */
      p -> nvec = dbg;
      if (p->debug&1)
	printf("# read %s, nvec=%d\n",argv[i],p-> nvec);
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
    else if (sscanf(argv[i],"ntot=%d",&dbg)==1){ /** `ntot` */
      p -> ntot = dbg;
      if (p->debug&1)
	printf("# read %s, ntot=%d\n",argv[i],p-> ntot);
    }    
    else if (sscanf(argv[i],"nfail=%d",&dbg)==1){ /** `nfail` */
      p -> nfail = dbg;
      if (p->debug&1)
	printf("# read %s, nfail=%d\n",argv[i],p-> nfail);
    }    
    else if (sscanf(argv[i],"seed=%d",&dbg)==1){ /** `seed` */
      p->seed=dbg;
      if (p->debug&1)
	printf("# read %s, seed=%d\n",argv[i],p->seed);
    }
    else if (sscanf(argv[i],"pmin=%lg",&val)==1){ /** `pmin` */
      p -> pmin = val;
      if (p->debug&1)
	printf("# read %s, pmin=%g\n",argv[i],p-> pmin);
    }
    else if (sscanf(argv[i],"pmax=%lg",&val)==1){ /** `pmax` */
      p -> pmax = val;
      if (p->debug&1)
	printf("# read %s, pmax=%g\n",argv[i],p-> pmax);
    }
    else if (sscanf(argv[i],"pstep=%lg",&val)==1){ /** `pstep` */
      p -> pstep = val;
      if (p->debug&1)
	printf("# read %s, pstep=%g\n",argv[i],p-> pstep);
    }
    else if (0==strncmp(argv[i],"f=",2)){/** back compatibility */
      if(strlen(argv[i])>2)
        p->fin = argv[i]+2;
      else
        p->fin = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, (fdem) f=%s\n",argv[i],p->fin);
    }
    else if (0==strncmp(argv[i],"fdem=",2)){
      if(strlen(argv[i])>2)
        p->fin = argv[i]+2;
      else
        p->fin = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fdem=%s\n",argv[i],p->fin);
    }
    else if (0==strncmp(argv[i],"fdet=",2)){
      if(strlen(argv[i])>2)
        p->fdet = argv[i]+2;
      else
        p->fdet = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fdet=%s\n",argv[i],p->fdet);
    }
    else if (0==strncmp(argv[i],"fobs=",2)){
      if(strlen(argv[i])>2)
        p->fobs = argv[i]+2;
      else
        p->fobs = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fobs=%s\n",argv[i],p->fobs);
    }
    else if((strcmp(argv[i],"--help")==0)
            ||(strcmp(argv[i],"-h")==0)
            ||(strcmp(argv[i],"help")==0)){
      printf( USAGE , argv[0],argv[0]);
      exit (-1);
    }
    else{ /* unrecognized option */
      printf("# unrecognized parameter \"%s\" at position %d\n",argv[i],i);
      ERROR("try \"%s -h\" for help",argv[0]);
    }
      
  }
  
  if (p->seed == 0){
    p->seed=time(NULL)+1000000ul*getpid(); /* ensure a different seed */
    if(p->debug)
      printf("# initializing seed=%d from time(NULL)\n",p->seed);
    /** use `tinymt64_generate_double(&pp.tinymt)` for double [0,1] */
  }  
  // srand(time(p->seed));
  tinymt64_init(&tinymt,p->seed);

  if((p->mode &2)&&(p->pmax-p->pmin<1e-7)){
      p->pmin=1e-2;
      p->pmax=0.101;
      p->pstep=1e-2;
  }
  
  return 0;
};

void var_kill(params_t *p){
  free(p->vP);
  free(p->vLLR);
  p->vP = p->vLLR = NULL;
  p->mH =  csr_free(p->mH);
  p->mHt = csr_free(p->mHt);
  p->mL =  csr_free(p->mL);
  p->mLt = csr_free(p->mLt);
}

int main(int argc, char **argv){
  params_t * const p=&prm;
  var_init(argc,argv, & prm); /* initialize variables */
  /** read in the error model file, initialize sparse matrices */
  //  one_prob_t **err_mod = NULL;
  read_dem_file(p->fin,p); /** creates all matrices */
  
#ifndef NDEBUG  
  if(p->debug & 64){ /** print matrices */
    mzd_t *mH0 = mzd_from_csr(NULL,p->mH);
    printf("matrix mH0:\n");  mzd_print(mH0);
    mzd_free(mH0); 
    
    mzd_t *mL0 = mzd_from_csr(NULL,p->mL);
    printf("matrix mL0:\n");  mzd_print(mL0);
    mzd_free(mL0); 
    
  }
  
#endif 
    
  double pmax=p->pmax +0.005*p->pstep;
  for(double prob = p->pmin; prob < pmax; prob += p->pstep){    
    if(p->mode&2)
      prob_init(p, prob);
    /** at least one round always */
    int rounds=(int )ceil((double) p->ntot / (double) p->nvec); 
    long int synd_tot=0, synd_fail=0;
    for(int iround=0; iround < rounds; iround++){
      if(p->debug &1)
        printf("# starting round %d of %d\n", iround, rounds);
      /** todo: add code for reading syndrome information from file */
      
      // decoder_init( & prm);
      mzd_t *mHe = mzd_init(p->nrows, p->nvec); /** each column a syndrome vector `H*e` */
      mzd_t *mLe = mzd_init(p->ncws,  p->nvec); /** each column `L*e` vector */
      p->maxJ = do_errors(mHe,mLe,p);

      if((p->debug & 512)&&(p->debug &4)){ /** print matrices */
        printf("matrix mLe:\n");  mzd_print(mLe);
        printf("matrix mHe:\n");  mzd_print(mHe);
      }
    
      // actually decode and generate error vectors 
      mzd_t *mE0=NULL;
#ifndef NDEBUG  /** need `mHe` later */
      mzd_t *mS=mzd_copy(NULL,mHe);
      mE0=do_decode(mS, p); /** each row a decoded error vector */
      mzd_free(mS); mS=NULL;
#else
      mE0=do_decode(mHe, p); /** each row a decoded error vector */
#endif /* NDEBUG */
      mzd_t *mE0t = mzd_transpose(NULL, mE0);
      mzd_free(mE0); mE0=NULL;
    
#ifndef NDEBUG 
      mzd_t *prodHe = csr_mzd_mul(NULL,p->mH,mE0t,1);
      mzd_add(prodHe, prodHe, mHe);
      if(!mzd_is_zero(prodHe)){
        if((p->debug&512)||(p->nvec <=64)){
          printf("syndromes difference:\n");
          mzd_print(prodHe);
        }
        ERROR("some syndromes are not matched!\n");
      }
      mzd_free(prodHe); prodHe = NULL;
      mzd_free(mHe);    mHe    = NULL;
#endif

      mzd_t *prodLe = csr_mzd_mul(NULL,p->mL,mE0t,1);

      if(p->debug & 512){ /** print matrices */
        printf("prodLe:\n");
        mzd_print(prodLe); 
        printf("mLe:\n");
        mzd_print(mLe); 
      }
  
      mzd_add(prodLe, prodLe, mLe);
      mzd_free(mLe); mLe=NULL;
    
      int fails=0;
      for(rci_t ic=0; ic< prodLe->ncols; ic++){
        rci_t ir=0;
        if(mzd_find_pivot(prodLe, ir, ic, &ir, &ic)){
          fails++;
          //      printf("ir=%d ic=%d fails=%d\n",ir,ic,fails);
        }
        else /** no more pivots */
          break;          
      }
      /** update the global counts */
      synd_tot  += prodLe->ncols;
      synd_fail += fails;
      mzd_free(prodLe); prodLe=NULL;
      if((p->nfail>0) && (synd_fail >= p->nfail))
        break;
    }
    
    if(p->mode&2)
      printf(" %g %ld %ld # %s\n",prob, synd_fail, synd_tot, p->fin);
    else 
      printf(" %ld %ld # %s\n",synd_fail, synd_tot, p->fin);
  }
  // clean up
  //  free(err_mod); err_mod=NULL;
  var_kill(p);
  if(p->mode&2)
    printf("\n\n");
  return 0;
}

