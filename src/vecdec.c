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

params_t prm={ .nrows=0, .n=0, .ncws=0, .steps=1, .debug=1, .fin=NULL,
  .seed=0, .colw=10, .mode=0, .maxJ=20, .nvec=16, 
  .vP=NULL, .vLLR=NULL, .mH=NULL, .mHt=NULL, .mLt=NULL};

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
/* note: force-inlining actually slows it down (???) */
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

void do_errors(mzd_t *mHe, mzd_t *mLe,
               const double vp[],
               const params_t * const p){
  int max = p->maxJ;  /** current size of `vec` */
  int * vec = malloc(p->maxJ * sizeof(int));
  //  for(int i=0; i < p->n; i++)
  //    printf("vP[%d]=%g%s",i,vp[i],i+1==p->n? "\n":"\t");

  /** for each error type (column of `mH` = row of `mHt`) */
  for(int i=0; i < p->n; i++){
    int ivec=0;
    /** prepare the list of syndrome columns to flip */
    double onebyL = -1.0/log(1.0-vp[i]);
    int j =(int )floor(onebyL * rnd_exponential());
    printf("p=%g onebyL=%g j=%d nvec=%d\n",vp[i],onebyL,j,p->nvec);
    if(j < p->nvec){/** otherwise we are to skip this error altogether */
      do{
        if(ivec >= max){
          max=2*max;
          vec=realloc(vec,max);
        }
        vec[ivec++]=j;
        j += (int )ceil(onebyL * rnd_exponential());         
      }
      while(j < p->nvec);
      
      for (j=0; j<ivec; j++)
        printf("vec[%d]=%d%s",j,vec[j],j+1==ivec?" ###\n":" ");
      
      /** flip the bits in `mHe` row by row to speed it up */
      for(int ir = p->mHt->p[i]; ir < p->mHt->p[i+1]; ir++){
        int irow = p->mHt->i[ir];
        for(j=0; j<ivec; j++)
          mzd_flip_bit(mHe,irow,vec[j]);
      }
      printf("i=%d j=%d\n",i,j);
      printf("matrix mHe:\n");  mzd_print(mHe);

      
      /** flip the bits in `mLe` row by row */
      for(int ir = p->mLt->p[i]; ir < p->mLt->p[i+1]; ir++){
        int irow = p->mLt->i[ir];
        for(j=0; j<ivec; j++)
          mzd_flip_bit(mLe,irow,vec[j]);
      }
    }
  }
}

/** todo: reuse `mE` matrix */
mzd_t *do_decode(mzd_t *mS, params_t const * const p){
  mzd_t * mH = mzd_from_csr(NULL, p->mH);
  //  printf("here:\n"); mzd_print(mH);

  mzd_t * mE = mzd_init(mH->ncols,mS->ncols);
  mzd_t * mEt= mzd_init(mS->ncols,mH->ncols);
  mzd_t * mEt0=mzd_init(mS->ncols,mH->ncols); /** best errors to output */
  double *vE = calloc(mS->ncols,sizeof(double));
  if((!mE) || (!mEt) || (!mEt0) || (!vE))
    ERROR("memory allocation failed!\n");
  //  printf("here:\n"); mzd_print(mH);

  mzp_t * perm=mzp_init(p->n); /** identity column permutation */
  mzp_t * pivs=mzp_init(p->n); /** list of pivot columns */
  if((!pivs) || (!perm))
    ERROR("memory allocation failed!\n");

  
  /** first pass (todo: order by decreasing `p`) */
  int rank=0;
  for(int i=0; i< p->n; i++){
    int col=perm->values[i];
    int ret=twomat_gauss_one(mH,mS, col, rank);
    if(ret)       
      pivs->values[rank++]=col;      
  }
  if(p->debug &4){ /** debug gauss */
    printf("rank=%d\n",rank);
    printf("perm: "); mzp_out(perm);
    printf("pivs: "); mzp_out(pivs);
    mzd_print(mH);
    mzd_print(mS);
  }
  // for each syndrome, calculate error vector and energy
  mzd_set_ui(mE,0); /** zero matrix */
  for(int i=0;i< rank; i++)     
    mzd_copy_row(mE,pivs->values[i],mS,i);
  mEt0 = mzd_transpose(mEt0,mE);
  for(int i=0; i< mS->ncols; i++)
    vE[i]=mzd_row_energ(p->vLLR,mEt0,i);

  /** main loop over permutations * ****************** */
  for (int ii=1; ii< p->steps; ii++){
    pivs=mzp_rand(pivs); /** random pivots LAPAC-style */
    mzp_set_ui(perm,1); perm=perm_p_trans(perm,pivs,0); /**< corresponding permutation */

    //    printf("permutation: "); mzp_out(perm);
    // gauss, record pivots
    rank=0;
    for(int i=0; i< p->n; i++){
      int col=perm->values[i];
      int ret=twomat_gauss_one(mH,mS, col, rank);
      if(ret)       
        pivs->values[rank++]=col;      
    }
    // for each syndrome, calculate errow vector and energy; update minima
    mzd_set_ui(mE,0); /** zero matrix */
    for(int i=0;i< rank; i++)     
      mzd_copy_row(mE,pivs->values[i],mS,i);
    mEt = mzd_transpose(mEt, mE);
    for(int i=0; i< p->n; i++){
      double energ=mzd_row_energ(p->vLLR,mEt0,i);
      if(energ< vE[i]){
        vE[i]=energ;
        mzd_copy_row(mEt0,i,mEt,i);
      }
    }
  }
  mzp_free(perm);
  mzp_free(pivs);
  mzd_free(mH);
  mzd_free(mE);
  mzd_free(mEt);
  free(vE);
  return mEt0;
}

one_prob_t ** read_file(char *fnam, params_t *p){
  ssize_t linelen, col=0;
  size_t lineno=0, bufsiz=0;
  char *buf = NULL;
  p->numH = p->numL = 0;

  if(p->debug&1)
    printf("opening file %s\n",fnam);  
  FILE *f = fopen(fnam, "r");
  if(f==NULL)
    ERROR("can't open the file %s\n",fnam);

  int cnt;
  do{ /** read lines one-by-one until `r k d` are found ****************************/
    int num;
    lineno++; col=0; linelen = getline(&buf, &bufsiz, f);
    if(linelen<0)
      ERROR("error reading line %zu of file %s\n",lineno,fnam);
    char *c=buf;
    while(isspace(*c)){ c++; col++; }
    if((*c == '\0')||(*c == '#')||(col >= linelen))
      cnt=0; /**  try next line */
    else{
      int r,k,n;
      num=0;
      cnt=sscanf(c," %d %d %d %n",&r,&k,&n,&num);      
      if(cnt!=3)
        ERROR("expected three integers ' r k n ' in a row\n"
              "%s:%zu:%zu: '%s'\n", fnam,lineno,col+1,buf);
      col+=num;
      p->nrows=r;
      p->ncws =k;
      p->n =n;
      if(p->debug &1)
        printf("# r=%d k=%d n=%d\n",r,k,n);
    }
  }  
  while(cnt==0);

  size_t max = (p->colw*sizeof(int)+ sizeof(  one_prob_t)); /** memory for a single error */
  one_prob_t **s = calloc(p->n,max);
  one_prob_t *pos = ( one_prob_t *) (s + p->n); /** leave space for row pointers */
  for(int i=0; i < p->n; i++){
    s[i] = pos;
    pos = (one_prob_t *) (( char *) pos + max);  
  }
  
  cnt = 0;
  int nread=0;
  do{ /** read lines one-by-one until `n` rows are read */
    int cnt=0, num, val;
    lineno++; col=0; linelen = getline(&buf, &bufsiz, f);
    if(linelen<0)
      ERROR("error reading line %zu of file %s\n",lineno,fnam);
    char *c=buf;
    while(isspace(*c)){ c++; col++; }
    if((*c == '\0')||(*c == '#')||(col >= linelen))
      cnt=0; /**  try next line */
    else{/** this row `must` contain a valid entry! */
      num=0;
      cnt = sscanf(c," %lg %n",& s[nread]->p, &num);
      if (cnt!=1)
        ERROR("expected a double followed by [[# ]...] ; [[# ]...], nread=%d\n"
              "%s:%zu:%zu: '%s'\n", nread,fnam,lineno,col+1,buf);
      if ((s[nread]->p <=0) || (s[nread]->p>=1))
        ERROR("expect probability=%g in (0,1) range, nread=%d\n"
              "%s:%zu:%zu: '%s'\n",s[nread]->p,nread,fnam,lineno,col+1,buf);
      c+=num; col+=num;
      int i=0;
      do{
        if(i >= p->colw)
          ERROR("too many entries in a row, increase colw=%d on command line\n"
                "%s:%zu:%zu: '%s'\n", p->colw,fnam,lineno,col+1,buf);
        cnt = sscanf(c," %d %n",&val, &num);
        if(cnt==1){
          //          printf("read %d num=%d\n",val,num);
          s[nread]->idx[i++] = val;
          if (s[nread]->n1)
            p->numL ++;
          else
            p->numH ++; 
          col+=num; c+=num;
          if(c[0]==';'){
            //            printf("semicolon i=%d num=%d\n",i,num);
            if(s[nread]->n1)
              ERROR("only one ';' in a row expected, nread=%d i=%d\n"
                    "%s:%zu:%zu: '%s'\n", nread,i,fnam,lineno,col+1,buf);
            else{
              s[nread]->n1=i;
              col++;
              c++;
            }
          }
        }
        else if((c[0]!='#')&&(c[0]!='\n')&&(c[0]!='\0'))
          ERROR("unexpected entry, nread=%d i=%d\n"
                "%s:%zu:%zu: '%s'\n", nread,i,fnam,lineno,col+1,buf);
        else{/** end of the line or comment starting */
          if(s[nread]->n1){
            s[nread++]->n2=i;
            break;
          }
          else
            ERROR("missing ';' in this row, nread=%d i=%d\n"
                  "%s:%zu:%zu: '%s'\n", nread,i,fnam,lineno,col+1,buf);
        }          
      }
      while(1);
    }    
  }
  while (nread < p->n);

  if(p->debug & 2){/** `print` out the entire error model ************************** */
    printf("# error model read: r=%d k=%d n=%d\n",p->nrows, p->ncws, p->n);
    for( int i=0; i < p->n; i++){
      one_prob_t *row = s[i];
      printf("# n1=%d n2=%d\n",row-> n1, row->n2);
      printf("%g ",row-> p);
      int j=0;
      for(  ; j< row->n1; j++)
        printf(" %d",row->idx[j]);
      printf(" ;");
      for(  ; j< row->n2; j++)
        printf(" %d",row->idx[j]);
      printf("\n");
    }          
  }
  
  if(p->debug&1)
    printf("closing file %s\n",fnam);
  fclose(f);
  if (buf) 
    free(buf);
  return s;
}

void mat_init(one_prob_t **in, params_t *p){
  double *vP = p->vP = malloc(p->n * sizeof(double));
  double *vLLR = p->vLLR = malloc(p->n * sizeof(double));
  csr_t *mH = p->mH   = csr_init(NULL, p->nrows, p->n,    p->numH);
  csr_t *mLt = p->mLt = csr_init(NULL, p->n,     p->ncws, p->numL); /** transposed */
  if((!vP) || (!vLLR) || (!mH) || (!mLt))
    ERROR("memory allocation failed!\n");
  int ipair1=0, ipair2=0;
  for (int i=0; i< p->n; i++){
    one_prob_t *row = in[i];
    vP[i] = row->p;
    vLLR[i] = row->p > MINPROB ? log((1.0/ row->p -1)) : log(1/MINPROB - 1);
    int j=0;
    for( ; j< row->n1; j++){
      mH->i[ipair1]  = i;    /** column */
      mH->p[ipair1++]= row->idx[j]; /** row */
    }
    for( ; j< row->n2; j++){
      mLt->p[ipair2]  = i;    /** column */
      mLt->i[ipair2++]= row->idx[j]; /** row */
    }
  };
  mH->nz  = p->numH;
  csr_compress(mH);

  p->mHt = csr_transpose(p->mHt, mH);
  
  mLt->nz = p->numL;
  csr_compress(mLt);
  
  if(p->debug & 2){ /** print resulting vectors and matrices */    
    for(int i=0; i< p->n; i++)
      printf("%g%s",vP[i],i+1<p->n?" ":"\n");
    for(int i=0; i< p->n; i++)
      printf("%g%s",vLLR[i],i+1<p->n?" ":"\n");
    csr_out(mH);
    csr_out(mLt);
  }
}

int local_init(int argc, char **argv, params_t *p){

  int dbg=0;
  
  for(int i=1; i<argc; i++){  /** `debug` */
    if(sscanf(argv[i],"debug=%d",& dbg)==1){
      if(dbg==0)
	p->debug = 0;
      else{
	p->debug ^= dbg;
	printf("# read %s, debug=%d octal=%o\n",argv[i],p->debug,p->debug);
      }
    }
    else if(sscanf(argv[i],"mode=%d",& dbg)==1){
      if(dbg==0)
	p->mode = 0;
      else{
	p->mode ^= dbg;
	printf("# read %s, mode=%d octal=%o\n",argv[i],p->mode,p->mode);
      }
    }
    else if (sscanf(argv[i],"nrows=%d",&dbg)==1){ /** `nrows` */
      p -> nrows = dbg;
      if (p->debug)
	printf("# read %s, nrows=%d\n",argv[i], p-> nrows);
    }
    else if (sscanf(argv[i],"ncws=%d",&dbg)==1){ /** `ncws` */
      p -> ncws = dbg;
      if (p->debug)
	printf("# read %s, ncws=%d\n",argv[i],p-> ncws);
    }
    else if (sscanf(argv[i],"steps=%d",&dbg)==1){ /** `steps` */
      p -> steps = dbg;
      if (p->debug)
	printf("# read %s, steps=%d\n",argv[i],p-> steps);
    }
    else if (sscanf(argv[i],"seed=%d",&dbg)==1){ /** `seed` */
      p->seed=dbg;
      if (p->debug)
	printf("# read %s, seed=%d\n",argv[i],p->seed);
    }
    else if (0==strncmp(argv[i],"f=",2)){
      p->fin = argv[i]+2;
      if (p->debug)
	printf("# read %s, finH=%s\n",argv[i],p->fin);
    }
    else if((strcmp(argv[i],"--help")==0)||(strcmp(argv[i],"-h")==0)){
      printf( USAGE , argv[0],argv[0]);
      exit (-1);
    }
    else{ /* unrecognized option */
      printf("# unrecognized parameter \"%s\" at position %d\n",argv[i],i);
      ERROR("try \"%s -h\" for options",argv[0]);
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
  
  return 0;
};

void local_kill(one_prob_t **copy, params_t *p){
  if(copy)
    free(copy);

  free(p->vP);
  free(p->vLLR);
  p->mH = csr_free(p->mH);
  p->mHt = csr_free(p->mHt);
  p->mLt = csr_free(p->mLt);
  p->vP = p->vLLR = NULL;
}

int main(int argc, char **argv){
  params_t * const p=&prm;
  local_init(argc,argv, & prm); /* initialize variables */

  /** read in the model file, initialize sparse matrices */
  one_prob_t **copy = read_file(p->fin, p );  /** @todo: sort events and optimize decoder model */
  mat_init(copy, p);
  
  // decoder_init( & prm);
  /** choose how many syndromes to use (depending on columns in `H`) */

  // copy entries to main matrix
  mzd_t * mH0=mzd_from_csr(NULL,p->mH);
  printf("matrix mH0:\n");  mzd_print(mH0);
  printf("xxx %d %d \n",p->nrows, p->nvec);

  // get or prepare the syndromes
  mzd_t *mHe = mzd_init(p->nrows, p->nvec);
  mzd_set_ui(mHe,0); /** zero matrix */

  mzd_t *mLe = mzd_init(p->ncws,  p->nvec);
  mzd_set_ui(mLe,0); /** zero matrix */
  do_errors(mHe,mLe,p->vP,p);
  printf("matrix mLe:\n");  mzd_print(mLe);
  
  // actually decode and generate error vectors (sparse)
    mzd_t *mE0=NULL; 
  mE0=do_decode(mHe, p); 
  mzd_print(mE0);

  // give the answer `L*e` for each set of syndromes in turn

  // clean up
  mzd_free(mH0);
  local_kill(copy, p);
  return 0;
}

