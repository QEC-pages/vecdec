
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
#include "utils.h"

tinymt64_t tinymt;

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
  if(fprintf(f,"%%MatrixMarket matrix array real general\n")<0)
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
