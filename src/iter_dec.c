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
#include <time.h>
#include <unistd.h>
#include <m4ri/m4ri.h>
#include "utils.h"
#include "util_m4ri.h"
#include "iter_dec.h"

params_t prm={ .nrows=0, .n=0, .ncws=0, .steps=50,
  .lerr=0, .useP=0.0, //.swait=0, //  .nvec=16,
  .ntot=1, .nfail=0, .seed=0, 
  .debug=1, .fdem=NULL, .fdet=NULL, .fobs=NULL, .finH=NULL, .finP=NULL,
  .finG=NULL, .finL=NULL, .internal=0, 
  .mode=0, 
  .LLRmin=0, .LLRmax=0, 
  .vP=NULL, .vLLR=NULL, .mH=NULL, .mHt=NULL,
  .mL=NULL, .mLt=NULL,
  .file_det=NULL, .file_obs=NULL, .line_det=0, .line_obs=0 };

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
    else if (sscanf(argv[i],"steps=%d",&dbg)==1){ /** `steps` */
      p -> steps = dbg;
      if (p->debug&1)
	printf("# read %s, steps=%d\n",argv[i],p-> steps);
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
    else if (sscanf(argv[i],"useP=%lg",&val)==1){ /** `useP` */
      p -> useP = val;
      if (p->debug&1)
	printf("# read %s, useP=%g\n",argv[i],p-> useP);
    }
    else if (0==strncmp(argv[i],"f=",2)){/** back compatibility */
      if(strlen(argv[i])>2)
        p->fdem = argv[i]+2;
      else
        p->fdem = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, (fdem) f=%s\n",argv[i],p->fdem);
    }
    else if (0==strncmp(argv[i],"fdem=",5)){
      if(strlen(argv[i])>5)
        p->fdem = argv[i]+5;
      else
        p->fdem = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fdem=%s\n",argv[i],p->fdem);
    }
    else if (0==strncmp(argv[i],"finH=",5)){
      if(strlen(argv[i])>5)
        p->finH = argv[i]+5;
      else
        p->finH = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finH=%s\n",argv[i],p->finH);
    }
    else if (0==strncmp(argv[i],"finP=",5)){
      if(strlen(argv[i])>5)
        p->finP = argv[i]+5;
      else
        p->finP = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finP=%s\n",argv[i],p->finP);
    }
    else if (0==strncmp(argv[i],"finL=",5)){
      if(strlen(argv[i])>5)
        p->finL = argv[i]+5;
      else
        p->finL = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finL=%s\n",argv[i],p->finL);
    }
    else if (0==strncmp(argv[i],"finG=",5)){
      if(strlen(argv[i])>5)
        p->finG = argv[i]+5;
      else
        p->finG = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, finG=%s\n",argv[i],p->finG);
    }
    else if (0==strncmp(argv[i],"fdet=",5)){
      if(strlen(argv[i])>5)
        p->fdet = argv[i]+5;
      else
        p->fdet = argv[++i]; /**< allow space before file name */
      if (p->debug&1)
	printf("# read %s, fdet=%s\n",argv[i],p->fdet);
    }
    else if (0==strncmp(argv[i],"fobs=",5)){
      if(strlen(argv[i])>5)
        p->fobs = argv[i]+5;
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
    p->seed=time(NULL)+1000000ul*getpid(); /* ensure a different seed */
    if((p->debug)&&(p->mode!=3))
      printf("# initializing seed=%d from time(NULL)+1000000ul*getpid()\n",p->seed);
    /** use `tinymt64_generate_double(&pp.tinymt)` for double [0,1] */
  }
  // srand(time(p->seed));
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
    p->n = p->mH->cols;
    p->nrows = p->mH->rows;
    p->ncws = p->mL->rows;
  }
  
  if(p->finH){
    if(! p->finL){
      if((p->fobs) || (p->fdet))
	ERROR("Without L matrix, cannot specify fdet=%s or fobs=%s\n",
	      p->fdet ? p->fdet : "",  p->fobs ? p->fobs : "");
      else
	p->internal=1;
    }
    if((! p->finL) && (! p->finG)){ /** only `H` matrix specified */
      p->classical=1;
      p->internal=1;
    }
    p->mH=csr_mm_read(p->finH, p->mH, 0);
    p->n = p->mH->cols;
    p->nrows = p->mH->rows;
  }
  /** at this point we should have `H` matrix for sure */
  assert(p->mH);
  
  if(p->finL){
    p->mL=csr_mm_read(p->finL, p->mL, 0);
    p->ncws = p->mL->rows;
    if(p->mL->cols != p->n)
      ERROR("column number mismatch L[%d,%d] and H[%d,%d]\n",
	    p->mL->rows, p->mL->cols, p->mH->rows, p->mH->cols);
  }

  if(p->finG){
    p->mG=csr_mm_read(p->finG, p->mG, 0);
    if(p->mG->cols != p->n)
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

    /** verify row orthogonality of `G` and `L` */
    if((p->mL)&&(product_weight_csr_mzd(p->mL,MGt,0)))
      ERROR("rows of L=Lx and G=Hz should be orthogonal \n");
    mzd_free(MGt);        
  }

  if(p->finP){
    int nrows, ncols, siz=0;
    p->vP=dbl_mm_read(p->finP, &nrows, &ncols, &siz, NULL);
    if(nrows*ncols != p->n)
      ERROR("invalid dimension=%d of P vector for H[%d,%d]\n",
	    nrows*ncols,p->mH->rows, p->mH->cols);
  }

  if(p->useP > 0){/** override */
    if(!p->vP)
      p->vP=malloc(p->n * sizeof(double));
    if(!p->vP)
      ERROR("memory allocation failed, n=%d",p->n);
    double *ptr=p->vP;
    for(int i=0;i<p->n;i++, ptr++)
      *ptr = p->useP;    
  }
  if(!p->vP)
    ERROR("probabilities missing, specify 'fdem', 'finP', or 'useP'");
  
  // rci_t linedet=0, lineobs=0;
  if(((p->fdet)&&(!p->fobs))||((!p->fdet)&&(p->fobs)))
    ERROR("need both fdet=%s and fobs=%s or none\n",
	  p->fdet ? p->fdet : "",  p->fobs ? p->fobs : "");
  
  if(p->fdet){/** expect both `fdet` and `fobs` to be defined */
    p->file_det=fopen(p->fdet, "r");
    if(p->file_det==NULL)
      ERROR("can't open the (det) file %s for reading\n",p->fdet);
    p->file_obs=fopen(p->fobs, "r");
    if(p->file_obs==NULL)
      ERROR("can't open the (obs) file %s for reading\n",p->fobs);
  }
  else{
    p->internal=1;
    /** TODO: generate a sufficient sample of errors (here?) */
    
  }

  if((p->debug & 64)&&(p->n < 128)){ /** print matrices */   
    assert(p->n != 0);
    if(p->mH){
      mzd_t *mH0 = mzd_from_csr(NULL,p->mH);
      printf("matrix mH0:\n");  mzd_print(mH0);
      mzd_free(mH0);
    }
    if(p->mG){
      mzd_t *mG0 = mzd_from_csr(NULL,p->mG);
      printf("matrix mG0:\n");  mzd_print(mG0);
      mzd_free(mG0);
    }
    if(p->mL){
      mzd_t *mL0 = mzd_from_csr(NULL,p->mL);
      printf("matrix mL0:\n");  mzd_print(mL0);
      mzd_free(mL0);
    }   
    if(p->vP){      
      for(int i=0; i < p->n; i++)
	printf(" P[%d]=%g \n",i,p->vP[i]);
    }
  }
  return 0;
}

void var_kill(params_t *p){
  if(p->file_det)
    fclose(p->file_det);
  if(p->file_obs)
    fclose(p->file_obs);
  
  if(p->vP)
    free(p->vP);
  if(p->vLLR)
    free(p->vLLR);
  p->vP = p->vLLR = NULL;
  p->mH =  csr_free(p->mH); /** OK if `NULL` */
  p->mHt = csr_free(p->mHt);
  p->mL =  csr_free(p->mL);
  p->mLt = csr_free(p->mLt);
  p->mG = csr_free(p->mG);
}
  
int main(int argc, char **argv){
  params_t * const p=&prm;
  var_init(argc,argv, p); /* initialize variables, read matrices, and open file (if any) */
  
  
  var_kill(p);
}
