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
#include <time.h>
#include <m4ri/m4ri.h>
#include "utils.h"
#include "util_m4ri.h"
#include "vecdec.h"

params_t prm={ .nrows=0, .n=0, .ncws=0, .steps=1, .debug=1, .fin=NULL, .seed=0};

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
    else if (0==strncmp(argv[i],"fin=",5)){
      p->fin = argv[i]+5;
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

int main(int argc, char **argv){
  local_init(argc,argv, & prm); /* initialize variables */

  /** read in the model file, initialize sparse matrices, sort events and optimize decoder model */
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

  return 0;
}

