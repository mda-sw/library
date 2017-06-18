#include <stdio.h>
#include <stdlib.h>


extern void HierClust(int *order, int *ia, int *ib, double *levels, 
			   const double *data, const double *wts, 
			   const int *nrows_ptr, const int *ncols_ptr);

int main(int argc, char *argv[])
{
  FILE *fp;
  char buffer[1024];
  char label[64];
  double x1, x2;
  int nrows, ncols = 2, index;
  double *data;
  double *wts;
  double *height;
  int i, *order, *ia, *ib;
  
  
  if (argc < 2) {
     printf("parameter: input_data\n");
     return 0;
  }
  
  fp = fopen(argv[1], "r");
  fgets(buffer, sizeof(buffer), fp);
  nrows = 0;
  while (fgets(buffer, sizeof(buffer), fp)!= NULL) {
    nrows++;
  }
  fclose(fp);
  
  
  data = (double*)malloc(sizeof(double)* 2 * nrows);
  wts = (double*)malloc(sizeof(double)* nrows);
  
  order = (int*)malloc(sizeof(int)* nrows);
  ia = (int*)malloc(sizeof(int)* (nrows-1));
  ib = (int*)malloc(sizeof(int)* (nrows-1));
  height = (double*)malloc(sizeof(double)* (nrows-1));
  
  for (i=0; i<nrows; i++) {
    wts[i] = 1.0;  
    order[i] = 0.0;       
  }

  for (i=0; i<(nrows-1); i++) {
        
    ia[i] = 0.0;  
    ib[i] = 0.0;
    height[i] = 0.0;
  }

  
  fp = fopen(argv[1], "r");
  fgets(buffer, sizeof(buffer), fp);
  index = 0;
  while (fgets(buffer, sizeof(buffer), fp)!= NULL) {
    sscanf(buffer, "%s %lf %lf", label, &x1, &x2);
 //   printf("%s\ %lf %lf\n", label, x1, x2);
    data[index++] = x1;
    data[index++] = x2;
  }
  fclose(fp);  
  printf("\n");
  
  HierClust(order, ia, ib, height, data, wts, &nrows, &ncols);
  
  printf("Order:\n");
  printf("======\n");
  for (i=0; i<nrows; i++) {
    printf("%d\n", order[i]);      
  } 
  printf("\n"); 
  
  printf("merge:\n");
  printf("======\n");
  for (i=0; i<(nrows-1); i++) {
    printf("%d %d\n", ia[i], ib[i]);      
  }  
  printf("\n");
   
  printf("height:\n");
  printf("=======\n");
  for (i=0; i<(nrows-1); i++) {
    printf("%lf or %e\n", height[i],  height[i]);      
  }  
  printf("\n");  
  
  free(data);
  free(wts);
  free(order);
  free(ia);
  free(ib);
  free(height);
  
  return 0;
}
