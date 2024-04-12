
/**
 * @file utils.c
 *
 * @brief Collection of utility functions
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2022 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 */
#include <strings.h>
#include <ctype.h>
#include <m4ri/m4ri.h>
#include "util_m4ri.h"
#include "utils.h"
#include "mmio.h"
#include "qllr.h"

tinymt64_t tinymt;

/** hashing storage helper functions ***** use `uthash.h` ***************************/


/** @brief print entire `one_vec_t` structure by pointer */
void print_one_vec(const one_vec_t * const pvec){
  printf(" w=%d E=%g cnt=%d [",pvec->weight, dbl_from_llr(pvec->energ),pvec->cnt);
  for(int i=0; i < pvec->weight; i++)
    printf("%d%s",pvec->arr[i], i+1 < pvec->weight ? " " :"]\n");
}

/** some extra io functions *******************************************************/

/** @brief open a new `NZLIST` file for writing. 
 * 
 *  `NZLIST file format:` 
 * header '%% NZLIST\n' (just one space)
 * zero or more comment lines starting with `%`
 * zero or more `entry` lines: 
 *   integer `w` (row weight) followed by exactly `w` strictly
 *   increasing positive integers representing non-zero positions
 *   of a binary vector starting with `1`.
 * `comment lines are only allowed at the top of the file.`
 */
FILE * nzlist_w_new(const char fnam[], const char comment[]){
  FILE *f=fopen(fnam,"w");
  if(!f)
    ERROR("can't open file %s for writing",fnam);
  fprintf(f,"%%%% NZLIST\n");
  if(comment)
    fprintf(f,"%% %s\n",comment);
  return f;
}

/** @brief write a `one_vec_t` entry to an open `NZLIST` file */
int nzlist_w_append(FILE *f, const one_vec_t * const vec){
  assert(vec && vec-> weight >0 );
  assert(f!=NULL);
  const int w=vec->weight;
  if(fprintf(f,"%d ",w)<=0)
    ERROR("can't write to `NZLIST` file");
  for(int i=0; i < w; i++)
    if(fprintf(f," %d%s", 1 + vec->arr[i], i+1 < w ? "" :"\n")<=0)
      ERROR("can't write to `NZLIST` file");
  return 0;
}

/** @brief prepare to read from an `NZLIST` file */
FILE * nzlist_r_open(const char fnam[], long int *lineno){
  int cnt;
  FILE *f=fopen(fnam,"r");
  if(!f)
    ERROR("can't open file %s for writing",fnam);
  if((EOF == fscanf(f,"%%%% NZLIST %n",&cnt)) || (cnt<9))
    ERROR("invalid signature line, expected '%%%% NZLIST'");
  *lineno=2;  
  char c=fgetc(f); /** are there comments to skip? */  //  putchar(c);
  while(c=='%'){ /** `comment line` starting with '%' */
    do{
      c=fgetc(f);      //      putchar(c);
      if(feof(f))
	return NULL;
    }
    while(c!='\n'); /** skip to the end of the comment line */
    (*lineno)++;
    c=fgetc(f);
  }
  ungetc(c,f); 
  /** we are at the start of a non-trivial entry or EOF */
  return f;
}

/** @brief read one item from an `NZLIST` file.
 * The structure in `vec` will be reallocated if necessary.
 * @param f an open file to read from.
 * @param[in] vec structure to hold the vector or NULL to allocate
 * @param fnam file name for debugging purposes
 * @param[in,out] lineno pointer to current line number in the file
 * @return the pointer to the structure containing the data or NULL. */
one_vec_t * nzlist_r_one(FILE *f, one_vec_t * vec, const char fnam[], long int *lineno){
  assert(f!=NULL);
  if ( ferror (f)|| feof(f) )
    return NULL; /** not an actuall error */
  //  printf("%s:%ld: start nzlist_r_one() here\n", fnam, *lineno);
  int w;
  if(!fscanf(f," %d",&w)){
    printf("%s:%ld: invalid NZLIST entry\n", fnam, *lineno);
    ERROR("expected an integer");
  }
  //  printf("read w=%d from line %ld\n",w,*lineno);
  if ((vec!=NULL) && (vec->weight<w)){
    free(vec);
    vec=NULL;
  }
  if(vec==NULL){
    vec = calloc(sizeof(one_vec_t)+w*sizeof(int), sizeof(char));
    if(!vec)
      ERROR("memory allocation");
  }
  vec->weight = w;
  vec->cnt = 1;
  for(int i=0; i<w; i++){
    if(!fscanf(f," %d ",vec->arr + i)){
      printf("%s:%ld: invalid entry of weight w=%d\n",fnam, *lineno, w);
      ERROR("expected an integer i=%d of %d",i,w);
    }    
    vec->arr[i]--; /** store as zero-based index */
  }

  for(int i=1; i<w; i++){ /** verify entry just read */
    if((vec->arr[i-1] < 0) || (vec->arr[i-1] >= vec->arr[i])){
      printf("%s:%ld: invalid entry of weight w=%d\n",fnam, *lineno, w);
      ERROR("expected strictly increasing positive entries");
    }   
  }
  (*lineno)++; /** in agreement with `NSLIST` file format */
  return vec;
}

/** @brief read an `MMX` array of doubles
 * With a column of doubles, `nrows` and `ncols` point to zeros on return.
 * @param[in] siz number of entries in the pre-allocated array 
 * @param[out] siz actual number of entries read
 * @param *arr must be allocated to size `siz` or NULL if `siz==0`.
 * @return the pointer to the array with the data
 */
double * dbl_mm_read(const char * const fin, int *nrows, int *ncols, int *siz, double *  arr){
  MM_typecode matcode;
  FILE *f;
  int num_items;
  *nrows = 0;
  *ncols = 0;

  if ((f = fopen(fin, "r")) == NULL) 
    ERROR("can't open file %s",fin);

  if (mm_read_banner(f, &matcode) != 0)
    ERROR("Could not process Matrix Market banner.");

  if (!(mm_is_matrix(matcode) && mm_is_dense(matcode) && 
	mm_is_real(matcode) && mm_is_general(matcode) )){
    printf("Sorry, this application does not support ");
    printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    ERROR("input file %s",fin);
    exit(1);
  }
  
  if(mm_read_mtx_array_size(f, nrows, ncols))
    ERROR("Can't read array size");
  if(arr != NULL){/** allocated array */
    if(*siz <= 0)
      ERROR("invalid input value siz=%d",*siz);
    if((*siz) < (*nrows) * (*ncols)){
      *siz = (*nrows) * (*ncols);
      arr = realloc(arr,sizeof(double)*(*siz));
      if(!arr)
	ERROR("memory allocation failed");
    }
  }
  else{
    if(*siz != 0)
      ERROR("invalid input value siz=%d, arr==NULL",*siz);
    *siz = (*nrows) * (*ncols);
    arr = malloc(sizeof(double)*(*siz));
    if(!arr)
      ERROR("memory allocation failed");
  }
  num_items = *siz;
  double *ptr = arr; 
  for(int i=0; i < num_items; i++, ptr++){
    if(1 != fscanf(f," %lg ", ptr))
      ERROR("failed to read entry %d of %d",i,num_items);
  }
  return arr;
}

void dbl_mm_write( char * const fout, const char fext[],
		   const int rows, const int cols, const double buf[],
		  const char comment[]){
  int result=0; /**< non-zero if write error */
  size_t len=strlen(fout)+strlen(fext)+1;
  char *str;
  
  FILE *f;
  if(strcmp(fout,"stdout")!=0){
    str=calloc(len, sizeof(char));
    sprintf(str,"%s%s",fout,fext);
    f=fopen(str,"w");
  }
  else{/** fout == "stdout" */
    str=fout;
    f=stdout;
  }
  
  if(!f)
    ERROR("can't open file '%s' for writing",str);
  if(fprintf(f,"%%%%MatrixMarket matrix array real general\n")<0)
    result++;
  if(comment!=NULL){
    if(fprintf(f,"%% %s\n",comment)<0)
      result++;
  }
  if(fprintf(f,"%d %d\n",rows,cols)<3)
    result++;
  const size_t num=rows*cols;
  for(size_t i=0;i<num;i++)
    if(fprintf(f,"%+18.12g\n", buf[i])<0)
      result++;
  if(result)
    ERROR("error writing to file '%s'",str);
  
  if(strcmp(fout,"stdout")!=0){
    fclose(f);
    free(str);
  }
}

/** @brief read detector error model (DEM) created by `stim`.
 * Immediately create CSR matrices `mH` and `mL` and vector `vP`; 
 * return the corresponding pointers via `ptrs` (in this order).
 * TODO: make it test for file type 
 * @param fnam file name for reading DEM from
 * @param p structure to store produced matrices
 */
void read_dem_file(char *fnam, void * ptrs[3], int debug){
  ssize_t linelen, col=0;
  size_t lineno=0, bufsiz=0; /**< buffer size for `readline` */
  char *buf = NULL;          /** actual buffer for `readline` */
  //  int nzH=0, nzL=0;  /** count non-zero entries in `H` and `L` */
  int maxH=100, maxL=100, maxN=100; 
  double *inP = malloc(maxN*sizeof(double));
  int_pair * inH = malloc(maxH*sizeof(int_pair));
  int_pair * inL = malloc(maxL*sizeof(int_pair));
  if ((!inP)||(!inH)||(!inL))
    ERROR("memory allocation failed\n");

  if(debug & 1)
    printf("# opening DEM file %s\n",fnam);
  FILE *f = fopen(fnam, "r");
  if(f==NULL)
    ERROR("can't open the (DEM) file %s for reading\n",fnam);

  int r=-1, k=-1, n=0;
  int iD=0, iL=0; /** numbers of `D` and `L` entries */
  do{ /** read lines one-by-one until end of file is found *************/
    lineno++; col=0; linelen = getline(&buf, &bufsiz, f);
    if(linelen<0)
      break;
    if(debug & 32) printf("# %s",buf);
    char *c=buf;
    double prob;
    int num=0, val;
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
            if(val>=r)
              r=val+1;  /** update the number of `D` pairs */
            if(iD>=maxH){
              maxH=2*maxH;
              inH=realloc(inH,maxH*sizeof(*inH));
            }
            inH[iD].a   = val;   /** add a pair */
            inH[iD++].b = n;
            if(debug & 32) printf("n=%d iD=%d val=%d r=%d\n",n,iD,val, r);
          }
          else if(sscanf(c," L%d %n",&val, &num)){/** `L` entry */
            c+=num; col+=num;
            assert(val>=0);
            if(val>=k)
              k=val+1;  /** update the number of `L` pairs */
            if(iL>=maxL){
              maxL=2*maxL;
              inL=realloc(inL,maxL*sizeof(*inL));
            }
            inL[iL].a   = val;   /** add a pair */
            inL[iL++].b = n;
            if(debug & 32) printf("n=%d iL=%d val=%d k=%d\n",n,iD,val,k);
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

  if(debug &1)
    printf("# read DEM: r=%d k=%d n=%d\n",r,k,n);

  ptrs[0] = csr_from_pairs(ptrs[0], iD, inH, r, n);
  ptrs[1] = csr_from_pairs(ptrs[1], iL, inL, k, n);
  ptrs[2] = inP;
  if (buf)
    free(buf);
  free(inH);
  free(inL);
  fclose(f);

}


#ifdef __MINGW32__ /** windows compiler */
int getline(char **line, size_t *n, FILE *fp){
  size_t size, oldsize;
  size_t len=0;
  if ((line==NULL) || (n==NULL) || (fp== NULL))
    ERROR("Null input parameters!");
  if (ferror (fp))
    return -1;
  if (feof(fp))
    return -1;
  for(oldsize=0, size=BUFSIZ ; ; size += BUFSIZ){
    /* BUFSIZ is "the optimal read size for this platform" */
    char *block;
      if (*n<=size){
	block = (char *) realloc(*line,size+1);
	/* realloc(NULL,n) equiv malloc(n) */
	if (block !=NULL){  /* reallocation worked */
	  *n=size+1;
	  *line=block;
	}
	else	   /* the original pointer may have been retained */
	  ERROR("memory allocation in getline");	
      }
      block = *line + oldsize; 
      if(NULL==fgets(block,BUFSIZ+1,fp)){
	if (oldsize>0)
	  break;  /* there is a string present */
	else
	  return -1;  /* nothing can be read - read error or EOF*/
      }
      len = oldsize+strlen(block);      
      if ((*line)[len-1]=='\n'){
	(*line)[--len]='\0';	
	break;  /* we are done */
      }
      if (feof(fp))
	break;
      oldsize=size; /* fgets puts a terminal '\0' on the end of the
		       string, so we make sure to overwrite this */
  } 
  return len;   
}

#endif /* __MINGW32__ */


