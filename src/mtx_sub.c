/**
 *  @file mtx_sub.c
 *
 * @brief mtx_cut - generate a submatrix of a binary mtx or alist matrix
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2023 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 *
 */

#include "utils.h"
#include "util_m4ri.h"

typedef struct PAR_T {
  char *fin, *out;
  int minR, maxR, minC, maxC;
  csr_t *mat;
  int debug;
} par_t;

#define USAGE							       \
  "%s: mtx_cut - generate a submatrix of a binary mtx or alist matrix\n" \
  "  usage: %s param=value [[param=value] ... ]\n"			\
  "\t Command line arguments are processed in the order given.\n"	\
  "\t Supported parameters:\n"						\
  "\t --help\t\t: give this help (also '-h' or just 'help')\n"		\
  "\t fin=[string]\t: name of the input file with the matrix (mm or alist)\n" \
  "\t out=[string]\t: name of the output file (or 'stdout', default value)\n" \
  "\t row=[int]\t: for this row, print the last non-zero column.\n"	\
  "\t rows=[int,int]\t: for this row block, print the last non-zero column.\n" \
  "\t minR, maxR, minC, maxC=[int]\t: min/max rows/columns (inclusive)\n" \
  "\t\t for the submatrix to generate.  Must be ordered minR <= maxR, minC <= maxC\n" \
  "\t\t and fit the dimensions of the original matrix.\n"		\
  "\t\t Missing parameters will be replaced by matrix dimensions or 0\n"	\
  "\t debug=[integer]\t: bitmap for aux information to output (default: 1)\n" \
  "\t Multiple 'debug' parameters are XOR combined except for 0.\n"	\
  "\t Space is allowed between 'fin=' / 'out=' and the file name string.\n"	\
  "\t Input file name 'fin' must be prior to parameters using the matrix.\n"	\
  "\t Output file name 'out' must be the last argument.\n"


//  "\t box=[int,int,int,int]\t: specify minR, maxR, minC, maxC (exclusive max)\n" 

par_t var={
  .fin=NULL,
  .out="std" "out",
  .mat=NULL,
  .minR=-1,
  .maxR=-1,
  .minC=-1,
  .maxC=-1,
  .debug=1
};

int var_init(int argc, char **argv, par_t *p){
  int dbg=0; 
  int cnt=0;

  for(int i=1; i<argc; i++){
    if(sscanf(argv[i],"debug=%d",& dbg)==1){/** `debug` */
      cnt++;
      if(cnt==1)
	p->debug = dbg; /** just assign if encountered once */
      else
	p->debug ^= dbg; /** otherwise `XOR` */
      if(p->debug &4)
	printf("# read %s, debug=%d octal=%o\n",argv[i],p->debug,p->debug);      
    }
  };

  if (argc==1)
    ERROR("try '%s --help'\n",argv[0]);
  
  for(int i=1; i<argc; i++){
    int pos=0;
    if(sscanf(argv[i],"debug=%d",& dbg)==1){
      /** `debug` - do nothing - already processed */
    }
    else if (0==strncmp(argv[i],"fin=",4)){ /** `fin` */
      if(strlen(argv[i])>4)
        p->fin = argv[i]+4;
      else
        p->fin = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, fin=%s\n",argv[i],p->fin);
      p->mat=csr_mm_read(p->fin, p->mat, 0, p->debug);
      if(p->debug&1)
	printf("# matrix dimensions %d by %d\n",p->mat->rows,p->mat->cols);
      //      if(p->debug&128)	csr_out(p->mat);
      if(p->debug&128){
	if(p->debug&1)
	  printf("# original matrix:\n");
	if (p->mat->cols < 100){
	  mzd_t *mmat = mzd_from_csr(NULL,p->mat);
	  mzd_print(mmat);
	  mzd_free(mmat);
	}
	else
	  csr_out(p->mat);
      }
    }
    else if (0==strncmp(argv[i],"out=",4)){ /** `out` */
      if(strlen(argv[i])>4)
        p->out = argv[i]+4;
      else
        p->out = argv[++i]; /**< allow space before file name */
      if (p->debug&4)
	printf("# read %s, out=%s\n",argv[i],p->out);
      if(i+1 != argc)
	ERROR("out=%s must be the last argument\n",p->out);
      if(!p->mat)
	ERROR("must provide the matrix, use fin=file_name\n");
      if (p->minR == -1)
	p->minR = 0;
      if (p->maxR == -1)
	p->maxR = p->mat->rows - 1;
      if (p->minC == -1)
	p->minC = 0;
      if (p->maxC == -1)
	p->maxC = p->mat->cols - 1;
      if(p->debug & 1)
	printf("creating submatrix rows %d to %d, cols %d to %d (inclusive)\n",  p->minR, p->maxR, p->minC, p->maxC);
      csr_t * sub = csr_submatrix(p->mat, p->minR, p->maxR+1, p->minC, p->maxC+1);
      //      if(p->debug&128)	csr_out(sub);
      if(p->debug&128){
	printf("# created submatrix:\n");
	if(sub->cols <100){
	  mzd_t *msub = mzd_from_csr(NULL,sub);
	  mzd_print(msub);
	  mzd_free(msub);
	}
	else
	  csr_out(sub);
      }
      char *comment;
      size_t size = snprintf(NULL, 0, "submatrix of %s minR=%d maxR=%d minC=%d maxC=%d\n", p->fin, p->minR, p->maxR, p->minC, p->maxC);
      if(!(comment = malloc(size + 1)))
	ERROR("memory allocation");
      snprintf(comment,size, "submatrix of %s minR=%d maxR=%d minC=%d maxC=%d",
	       p->fin, p->minR, p->maxR, p->minC, p->maxC);
      csr_mm_write(p->out,"",sub, comment);
      free(comment);
      csr_free(sub);      
    }    
    else if  (sscanf(argv[i],"minR=%d",&dbg)==1){ /** `minR` param */
      p -> minR = dbg;
      if (p->debug&4)
	printf("# read %s, minimum row %d\n",argv[i],p-> minR);
      if(!p->mat)
	ERROR("must provide the matrix first, use fin=file_name\n");
      if((dbg<0) || (dbg>= p->mat->rows))
	ERROR("minR=%d must be in the range from 0 to rows-1=%d\n", dbg,p->mat->rows-1);
      if((p->maxR !=-1) && (dbg > p->maxR))
	ERROR("minR must not exceed maxR=%d\n", p->maxR);
    }
    else if (sscanf(argv[i],"maxR=%d",&dbg)==1){ /** `maxR` param */
      p -> maxR = dbg;
      if (p->debug&4)
	printf("# read %s, maximum row %d\n",argv[i],p-> maxR);
      if(!p->mat)
	ERROR("must provide the matrix first, use fin=file_name\n");
      if((dbg<0) || (dbg>= p->mat->rows))
	ERROR("maxR=%d must be in the range from 0 to rows-1=%d\n", dbg, p->mat->rows-1);
      if((p->minR !=-1) && (dbg < p->minR))
	ERROR("maxR must not be smaller than minR=%d\n", p->minR);
    }
    else if  (sscanf(argv[i],"minC=%d",&dbg)==1){ /** `minC` param */
      p -> minC = dbg;
      if (p->debug&4)
	printf("# read %s, minimum col %d\n",argv[i],p-> minC);
      if(!p->mat)
	ERROR("must provide the matrix first, use fin=file_name\n");
      if((dbg<0) || (dbg>= p->mat->cols))
	ERROR("minC=%d must be in the range from 0 to cols-1=%d\n",dbg, p->mat->cols-1);
      if((p->maxC !=-1) && (dbg>p->maxC))
	ERROR("minC must not exceed maxC=%d\n", p->maxC);
    }
    else if  (sscanf(argv[i],"maxC=%d",&dbg)==1){ /** `maxC` param */
      p -> maxC = dbg;
      if (p->debug&4)
	printf("# read %s, maximum col %d\n",argv[i],p-> maxC);
      if(!p->mat)
	ERROR("must provide the matrix first, use fin=file_name\n");
      if((dbg<0) || (dbg>= p->mat->cols))
	ERROR("maxC=%d must be in the range from 0 to cols-1=%d\n", dbg,p->mat->cols-1);
      if((p->minC !=-1) && (dbg < p->minC))
	ERROR("maxC must not be smaller than minC=%d\n", p->minC);
    }
    else if  (sscanf(argv[i],"row=%d",&dbg)==1){ /** `row` to check */
      if (p->debug&4)
	printf("# read %s, single row to check %d\n",argv[i],dbg);
      if(!p->mat)
	ERROR("must provide the matrix first, use fin=file_name\n");
      if(dbg >= p->mat->rows)
	ERROR("should be a valid matrix row < %d\n",p->mat->rows);
      if (p->mat->p[dbg+1] == p->mat->p[dbg] )
	ERROR("row %d is zero\n",dbg);
      long long pair = csr_min_max_blk(p->mat, dbg, dbg);
      int min = pair % p->mat->cols;
      int max = pair / p->mat->cols;
      if(p->debug & 1)
	printf("# row %d, first and last non-zero columns:\n",dbg);
      printf("%d %d\n", min, max);      

    }
    else if (sscanf(argv[i],"rows=%d,%n",&dbg,&pos)==1){ /** `row1, row2` to check */
      int r1=dbg, r2;
      if(sscanf(argv[i]+pos,"%d",&dbg)!=1)
	ERROR("read argv[%d]=%s, must specify two integers w/o spaces, 'rows=[int],[int]'\n", i,argv[i]);
      r2=dbg;
      if (p->debug&4)
	printf("# read %s, row range to check %d, %d\n",argv[i],r1,r2);
      if(!p->mat)
	ERROR("must provide the matrix first, use fin=file_name\n");
      long long pair = csr_min_max_blk(p->mat, r1, r2);
      int min = pair % p->mat->cols;
      int max = pair / p->mat->cols;
      if(p->debug & 1)
	printf("# row range [%d .. %d], first and last non-zero columns:\n",r1,r2);
      printf("%d %d\n", min, max);      
    }
    else if((strcmp(argv[i],"--help")==0)
            ||(strcmp(argv[i],"-h")==0)
            ||(strcmp(argv[i],"help")==0)){      
      printf( USAGE , argv[0],argv[0]);
      exit (-1);
    }    
    else{
      printf("unrecognized parameter argv[%d]=%s\n",i,argv[i]);
      ERROR("try '%s --help'\n",argv[0]);	   
    }
  }
  return 0;
}
  

int main(int argc, char **argv){
  var_init(argc,argv, &var);

  if(var.mat) csr_free(var.mat);
  return 0;
}
