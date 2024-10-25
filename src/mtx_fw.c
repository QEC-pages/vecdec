/**
 *  @file mtx_sub.c
 *
 * @brief mtx_qc - generate single-weight errors
 * 
 * 01 strings with fixed weight are sent to `stdout`.
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2024 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 *
 */

#include <assert.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <m4ri/m4ri.h>
#include "mmio.h"
#include "utils.h"
#include "util_m4ri.h"

#include "vec.h"

#define MAX_W 400

typedef struct PAR_T {
  int num; /**< how many strings (default: 0, all) */
  int nz[MAX_W]; /**< explicitly specified non-zero positions, starting with 0 */
  int n; /**< string length (default: 0, must set) */
  int w; /**< row weight (default: 0, must set) */
  //  int random; /**< generate random strings (default: 0, not random) */
  long long int seed; /**< seed for rng (default: 0, use `time(NULL)` */
  int debug; /** default: 1 */
} par_t;

par_t params={
  .num=0,
  .n=0,   //! must specify 
  .w=0,   //! must specify
  .nz={ -1 },
  //  .random=0,
  .seed=0,
  .debug=1
};


#define USAGE                                                          \
  "%s: mtx_fw - generate fixed-weight 01 strings\n" \
  "  usage: %s param=value [[param=value] ... ]\n"                      \
  "\t Parameters are read in the order specified.\n"                    \
  "\t Supported parameters:\n"                                          \
  "\t --help\t\t: give this help (also '-h' or just 'help')\n"          \
  "\t n=int\t: length of each output string (must be specified)\n"      \
  "\t w=int\t: weight of each output string\n"				\
  "\t nz=int[,int ...]\t: explicitly specified nz positions\n"		\
  "\t num=int\t: how many strings (choose randomly if not all)\n"       \
  "\t seed=int\t: RNG seed (default: 0, use 'time(NULL)')\n"            \
  "\t debug=int\t: debug bitmap (default: 1)\n"				\
  "\t either 'nz' or 'w' has to be specified\n"


int do_scan_nz(int i, const char argvi[], int val, int pos, par_t * const p){

  int w=0, dbg=val;
  const char *c = argvi;
  do{
    if((dbg<0)||(dbg>=p->n))
      ERROR("arg[%d]='%s' : invalid position(%d)=%d , n=%d\n",
	    i,argvi,w,dbg,p->n);
    if((w>0)&&(dbg<=p->nz[w-1]))
      ERROR("arg[%d]='%s' pos(%d)=%d pos(%d)=%d :"
	    " must be an increasing sequence! ",
	    i,argvi,w-1,p->nz[w-1],w,dbg);
    if(w>=MAX_W)
      ERROR("arg[%d]='%s' : too many non-zero entries w>=%d MAX_W=%d\n", i,argvi,w,MAX_W);
    p->nz[w++]=dbg;
    c+=pos;
  }
  while (sscanf(c,",%d%n",& dbg, &pos)==1);
  if ((p->w!=0) && (p->w != w))
    ERROR("arg[%d]='%s' : explicitly specified w=%d must match weight=%d of nz vector",i,argvi,p->w,w);
  if(p->debug&4){
    printf("# read argv[%d]='%s' , setting nz=[ ",i,argvi);
    for(int ii=0; ii< w; ii++)
      printf("%d%s",p->nz[ii],  ii+1==w ? " ] ":" ");
    printf(" wght=%d\n",w);
  }
  return w;
}


int var_init(int argc, char **argv, par_t *p){
  
  if (argc==1)
    ERROR("try '%s --help'\n",argv[0]);

  for(int i=1; i<argc; i++){
    int dbg=0, cnt=0;
    if(sscanf(argv[i],"debug=%d",& dbg)==1){/** `debug` */
      cnt++;
      if(cnt==1)
        p->debug = dbg; /** just assign if encountered once */
      else
        p->debug ^= dbg; /** otherwise `XOR` */
      if(p->debug &4)
        fprintf(stderr,"# read %s, debug=%d octal=%o\n",argv[i],p->debug,p->debug);      
    }
  };

  for(int i=1; i<argc; i++){
    int dbg=0;
    if (sscanf(argv[i],"w=%d",& dbg)==1){ /** `w` */
      p -> w = dbg;
      if (p->debug&4)
        fprintf(stderr,"# read %s, setting w=%d\n",argv[i],p->w);
    }
  };

  for(int i=1; i<argc; i++){   /** scan for everything else */
    int dbg, p1;
    long long int lldbg;
    if (sscanf(argv[i],"n=%d",& dbg)==1){ /** `n` */
      p -> n = dbg;
      if (p->debug&4)
        fprintf(stderr,"# read %s, setting n=%d\n",argv[i],p->n);
    }
    else if (sscanf(argv[i],"nz=%d%n",& dbg,&p1)==1){ /** `nz` */
      //      p -> nz[0] = dbg;
      p->w = do_scan_nz(i,argv[i],dbg,p1,p);
      if (p->debug&4)
        fprintf(stderr,"# read %s, setting nz=\n",argv[i]);
    }
    else if (sscanf(argv[i],"num=%d",& dbg)==1){ /** `num` (random if set) */
      p -> num = dbg;
      if (p->debug&4)
        fprintf(stderr,"# read %s, setting num=%d\n",argv[i],p->num);
    }
    else if (sscanf(argv[i],"seed=%lld",& lldbg)==1){ /** `seed` */
      p -> seed = lldbg;
      if (p->debug&4)
        fprintf(stderr,"# read %s, setting seed=%lld\n",argv[i],p->seed);
    }
    else if((strcmp(argv[i],"--help")==0)
            ||(strcmp(argv[i],"-h")==0)
            ||(strcmp(argv[i],"help")==0)){
      fprintf(stderr, USAGE , argv[0],argv[0]);
      exit (-1);
    }
    else if ((strncmp(argv[i],"debug=",6)!=0) &&
	     (strncmp(argv[i],"w=",2)!=0)){
      fprintf(stderr,"parameter argv[%d]='%s' is not recognized.\n",i,argv[i]);
      ERROR("please run '%s -h ' or '%s --help' for help",argv[0],argv[0]);
    }
  }
    
  if ((p->seed <= 0)&&(p->num > 0)){
      long long int seed_old= - p->seed; 
      p->seed=(seed_old) + time(NULL)+1000000ul*getpid(); /* ensure a different seed even if started at the same time */
      if(p->debug)
        fprintf(stderr,"# initializing seed=%lld from time(NULL)+1000000ul*getpid()+%lld\n",p->seed, seed_old);
      /** use `tinymt64_generate_double(&pp.tinymt)` for double `[0,1)` */
    }
  
    tinymt64_init(&tinymt,p->seed);

    if((p->n<=0) || (p->w<=0) || (p->w > p->n)){
      ERROR("n=%d and w=%d must be positive and w<=w\n",p->n,p->w);
    }

    if(p->num<0)
      ERROR("num=%d must be non-negative\n",p->num);

    return 0;
}

int main(int argc, char **argv){

  par_t *p = &params;

  var_init(argc,argv,p);

  const int w = p->w;
  vec_t *v = vec_init(p->n);

  if (p->nz[0]>=0){
    for(int i=0; i<p->w; i++)
      vec_push(p->nz[i],v);
    write_01_vec(stdout,v,v->max,"stdout");
    if(p->debug&1)
      fprintf(stderr,"wrote one vector specified by nz coefficients to stdout\n");
  }
  else if(p->num){ /** generate `num` random vectors of weight `w` */
    mzp_t *pivots = mzp_init(p->n); 
    mzp_t *perm   = perm_p(NULL, pivots,0); /* actual permutation */
    
    for (int i=0; i<p->num; i++){
      pivots = mzp_rand(pivots); /* LAPAC-style random node permutation */
      perm   = perm_p(perm, pivots,0); /* actual random permutation */

      /* construct random vec and print */
      v->wei=0;
      for(int i=0; i<p->w; i++)
	vec_push(perm->values[i],v);
      qsort(v->vec, w, sizeof(rci_t), cmp_rci_t);
      write_01_vec(stdout,v,v->max,"stdout");
    }
    //    free(v);
    mzp_free(perm);
    mzp_free(pivots);
  }
  else{
    long long int cnt=0;
    for(int i=0; i<p->w; i++)
      vec_push(i,v);
    do{
      write_01_vec(stdout,v,v->max,"stdout");
      cnt++;
    }
    while(vec_next(v));
    if(p->debug&1)
      fprintf(stderr,"wrote %lld vectors to stdout\n",cnt);
  }  
  free(v);
}


