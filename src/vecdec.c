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
#include <m4ri/m4ri.h>
#include "utils.h"
#include "util_m4ri.h"
#include "vecdec.h"

params_t prm={ .nrows=0, .n=0, .ncws=0, .steps=1, .debug=1, .fin=NULL, .seed=0, .colw=10};

/** 
 * @brief one step of complete gauss on column `idx` 
 * @param st pointer to the structure in question 
 * @param idx index of the device to deal with 
 * @param begrow row to start with 
 * @return number of pivot points found, i.e., 1 or 2 
 */
static inline int mzd_gauss_one(mzd_t *M, const int idx, const int begrow){
/* note: force-inlining actually slows it down (???) */
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
  return pivots;
  // if one, need to update the current pivot list 
}

one_prob_t ** read_file(char *fnam, params_t *p){
  ssize_t linelen, col=0;
  size_t lineno=0, bufsiz=0;
  char *buf = NULL;

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

  if(p->debug & 2){/** print out the entire error model */
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

int local_init(int argc, char **argv, params_t *p){

  int dbg=0;
  
  for(int i=1; i<argc; i++){  /** `debug` */
    if(sscanf(argv[i],"debug=%d",& dbg)==1){
      if(dbg==0)
	p->debug = 0;
      else{
	p->debug |= dbg;
	printf("# read %s, debug=%d octal=%o\n",argv[i],prm.debug,prm.debug);
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
      prm.seed=dbg;
      if (prm.debug)
	printf("# read %s, seed=%d\n",argv[i],prm.seed);
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
    if(p->debug)
      printf("# initializing rng from time(NULL)\n");
    srand(time(NULL));
  }
  else {
    srand(p->seed);
    if(p->debug)
      printf("# setting srand(%d)\n",prm.seed);
  }
  return 0;
};

void local_kill(one_prob_t **copy){
  if(copy)
    free(copy);
}

int main(int argc, char **argv){
  local_init(argc,argv, & prm); /* initialize variables */

  /** read in the model file, initialize sparse matrices, sort events and optimize decoder model */
  one_prob_t **copy = read_file(prm.fin, &prm );
  // decoder_init( & prm);

  /** choose how many syndromes to use (depending on columns in `H`) */

  // copy entries to main matrix

  /** main loop */
  for (int ii=0; ii< prm.steps; ii++){
    // generate random permutation

    // gauss, record pivots 

    // for each syndrome, calculate e (sparse) and energy; update minima

  }

  // give the answer `L*e` for each set of syndromes in turn

  // clean up
  local_kill(copy);
  return 0;
}

