/**
 *  @file mtx_sub.c
 *
 * @brief mtx_01_to_nz - filter 01 strings to nz format (stdin to stdout)
 *
 * @author Leonid Pryadko (University of California, Riverside)
 *
 * Copyright (C) 2024 Leonid Pryadko
 * University of California, Riverside
 * All rights reserved.
 *
 *
 */
#include <stdio.h>

int main(void){
  int row=0, w=0, len=0;
  int ok=1;
  char c=getchar();
  while (c!=EOF){
    if(ok==0){
      putchar(c);
      if(c=='\n')
	ok=1;
    }
    else
      switch(c){
      case '0':
	len++;
	break;
      case '1':
	printf("%s%d",w>0?",":"", len);
	w++;
	len++;
	break;
      case '\n':
	//	putchar('\n');
	printf(" # len=%d w=%d\n",len,w);
	//      fprintf(stderr,"row=%d len=%d w=%d\n",row+1,len,w);
	w=0;
	len=0;
	row++;
	break;
      default:
	ok=0; //! just echo non-01 strings 
	putchar(c);
	break;
      }
    c=getchar();
  }
  

  return 0;
}
