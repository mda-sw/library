/* 

Determine the gamma measure of ultrametricty of time series, 
for varying values of sliding window length.
Currently window lengths are: 5, 10, 15, ... 100, 105, 110.
We start at the beginning of the time series, and go as far as 
we have a full window length.
Squared Euclidean distances are ordinal coded.  

equil-time-series.c

gcc -lm equil-time-series.c -o equil-time-series
equil-time-series 1326 2 ftse1326.dat
equil-time-series 1326 2 r

Notes: 
Best - highest - gamma values have been found for 2 category 
coding.
Secondly, the number of input data values given is to be less than 
or equal to the number of values to be read.
Thirdly, above provides an example of a data file, and use of 
random data.

FM, 2004/5/16

*/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>    /* RAND_MAX is defined here; usually 32767 */
#define NRANSI
#include "nrutil.h"
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
#define NR_END 1
#define FREE_ARG char*
#define EPS 1.0E-10

/* nearest integer function */
# define NINT(x)  ((x >= 0.) ? (long) (x + 0.5) : (long) (x - 0.5))
# define MAX(a,b) (a >= b ? a : b)
# define MIN(a,b) ((a) <= (b) ? (a) : (b))

main(argc, argv)
int argc;
char *argv[];

{
FILE *stream;
unsigned long  n, m, ncat, j, k, ind, ind2, klo; 
unsigned long morig, ind0, zerocnt, nsam, totseg, nseg,loop, len; 
unsigned long i, i2, nnmin1div2; 
float tempmax, tempmin, ang1, ang2, umcount;
float mini,maxi,minj,maxj,mink,maxk;
float *data, **matrix(), **dist, *vector(), *catcount;
float **datax, xminval, xmaxval, inddist, ind2dist;
float in_value, med, d1, d2, d3, eucl();
float tripletVals[3], tripletTemp[3], tripletInd[3], tempval, *colsums;
double nnondegen, coeff, nsamples; 
void free_matrix(), free_vector(), piksr2(), sort2();
char option, *strncpy();
int check = 0, ival, jval, seed; 
/* Use check = 1 for more ouput diagnostics.  Else = 0. */

   if (argc < 4)  {
      printf("Syntax help: equil-time-series #vals ncat filename seed \n");
      exit(1);
      }
   n = atoi(argv[1]);              /* # values */
   ncat = atoi(argv[2]);           /* # no. of categories for each val */
   strncpy(&option,argv[3],1);     /* Simulated data option */
   seed = 1131;                    /* Used for random data generation */ 
   if (argc == 5) seed = atoi(argv[4]);  /* Seed value for random data */

   if (check == 1) {
      printf("No. vals., ncat = %d %d  \n", n, ncat);
      printf("Input file name or r = runif-simulate: %s.\n",argv[1]);
   }

   data = vector(n);  /* Storage allocation for input data */
   catcount = vector(ncat);    /* Cardinalities of categories */
   tempmax = -1.0; tempmin = +1.e5; 
   switch(option) {
       case 'r': 
             if (check == 1) printf("Will simulate...\n");
             srand(seed);  
             for (i = 1; i <= n; i++)  {
		 data[i] = rand();
                 if (data[i] > tempmax) tempmax = data[i]; 
                 if (data[i] < tempmin) tempmin = data[i]; 
                 data[i] /= (float)RAND_MAX;
             }
             if (check == 1) {
	     printf("Check on random no. generation: [0, 32767] assumed.\n");
	     printf("Max, min uniform deviates found = %12.4f %12.4f\n", 
                tempmax, tempmin);
             } 
             break; 
       default:
             if ((stream = fopen(argv[3],"r")) == NULL) {
                fprintf(stderr, "Program %s : cannot open file %s\n",
                       argv[0], argv[3]);
                exit(1);
             }
             /* ----------- Now read in data ------------- */
             for (i = 1; i <= n; i++)  {
                 fscanf(stream, "%f", &in_value);
                 data[i] = in_value;
             }
             break;
   }

   /* Check on (part of) input data. */
   if (check == 1) { 
       printf("Check on input data, first 30 values.\n"); 
           for (i = 1; i <= 30; i++) printf("%7.1f", data[i]);  
           printf("\n"); 
   }



   /* Take sets of len values in succession; and build distance matrix. */

   for (loop = 1; loop <= 22; loop++) {
       len = loop*5; 

       totseg = (int)((float)n/(float)len); 
       dist  = matrix(len, len); 
       datax = matrix(len, len); 

       coeff = 0.0;
       nnondegen = 0.0;
       nsamples = 0.0;

       for (nseg = 1; nseg <= totseg; nseg++) {

           /* Take values 1 to len for each window in turn */
           for (ival = len*(nseg-1)+1; ival <= len*(nseg-1)+len; ival++) {
               for (jval = len*(nseg-1)+1; jval <= len*(nseg-1)+len; jval++) {
                   dist[ival-(nseg-1)*len][jval-(nseg-1)*len] = 
                   pow(data[ival] - data[jval],2);
	       }
           }

           /* ----------- Recode -------------------------- */
           if (ncat > 1) {
              xmaxval = -1.e12;
              xminval = 1.e12;
              for (i=1; i<=len; i++) {
	          for (j=1; j<= len; j++) {
                      if (dist[i][j] < xminval) xminval = dist[i][j];
                      if (dist[i][j] > xmaxval) xmaxval = dist[i][j];
                  }
              }

              for (i=1; i<=len; i++) {
	          for (j = 1; j <= len; j++) {
                      tempval = 
                       (float)ncat*(dist[i][j]-xminval)/(xmaxval-xminval);
                      ind = (unsigned long)tempval;
                      /* Rare case of tempval=ncat */
                      if (ind == ncat) ind = ncat-1; 
	              ind += 1;                 /* So now ind = 1,2,...ncat */
                      datax[i][j] = ind;
                      if (i == j) datax[i][j] = 0.0;   /* Diagonal */
                  }
              }

              /* Return to use of 'dist', whether categorized like 
               we have just done here, or uncat'd, i.e. used as input. */
              for (i=1; i<=len; i++) {
	          for (j = 1; j <= len; j++) dist[i][j] = datax[i][j];
              }
	   }

           umcount = 0.0; 
           tempval = 0.0; 
           mini    = 1.e12; 
           maxi    = -1.e12;
           minj    = 1.e12; 
           maxj    = -1.e12;
           mink    = 1.e12; 
           maxk    = -1.e12;

           for (i=1; i<= len; i++) { 
               for (j=1; j<=len; j++) { 
                   for (k=1; k<=len; k++) { 

                       /* We need to consider (i,j), (j,k) and (i,k) */
                       d1 = dist[i][j];
                       d2 = dist[j][k];
                       d3 = dist[i][k];

                       nsamples += 1.0; 
                       if (d1 > EPS && d2 > EPS && d3 > EPS) {
                          nnondegen += 1.0; 

	                  /* Determine cosines of 3 angles in triangle */
                          tripletTemp[1] =(d1*d1 + d2*d2 - d3*d3)/(2.0*d1*d2);
                          tripletTemp[2] =(d2*d2 + d3*d3 - d1*d1)/(2.0*d2*d3);
                          tripletTemp[3] =(d1*d1 + d3*d3 - d2*d2)/(2.0*d1*d3);
                          tripletInd[1] = 1;
                          tripletInd[2] = 2;
                          tripletInd[3] = 3;
           
			  /* Sort the cosine values */
                          piksr2(3, tripletTemp, tripletInd);

                          /* Note: cos 60 deg = 0.5
                             60 deg = 1.0472 rad 
                             1 deg = 0.01745328 rad
                             2 deg = 0.03490656 */

			  /* if largest cosine angle is leq 60 deg. and
			     largest cosine angle is gt 0 then do */
                          if (tripletTemp[3] >= 0.5 && tripletTemp[3] < 1.0) {
		             /* So we have found the smallest angle, 
                                i.e. largest cosine 
		                Next, we will take diff. of angles (rad.) 
                                rather than cosines, of 2 other angles */
                             ang1 = acos(tripletTemp[1]);
                             ang2 = acos(tripletTemp[2]);
                             if (fabs(ang1 - ang2) < 0.03490656 ) 
                                coeff += 1.0;
                          }

                       }  /* End of: if d1, d2, d3 > EPS */

	           }  /* End of: k-loop */
	       }   /* End of: j-loop */
           }    /* End of: i-loop */

       } /* End of: nseg-loop */

       if (check == 1)
          printf("Len= %d gamma 1UM, 0NUM: %9.6f %14.0f %14.0f \n",
            len, coeff/nnondegen, nnondegen, nsamples); 
       printf("%9.6f \n", coeff/nnondegen); 


       free_matrix(datax, len, len); 
       free_matrix(dist, len, len);
   } /* End of len loop */


   /* ----------------- Free up some storage ------------------------*/
   free_vector(data, n);

}  /* End of: equil-time-series */


/** Euclidean distance: creation  *****************************/

float eucl(data, n, m, ival, jval)
float **data;
unsigned long ival, jval, n, m;
{
unsigned long j;
float distval;

     distval = 0.0; 
     for (j = 1; j <= m; j++) 
       distval +=(data[ival][j]-data[jval][j])*(data[ival][j]-data[jval][j]); 
     distval = sqrt(distval);

return distval;

}   /* End of: eucl */



/**  Error handler  **************************************************/
void erhand(err_msg)
char err_msg[];
/* Error handler */
{
    fprintf(stderr,"Run-time error:\n");
    fprintf(stderr,"%s\n", err_msg);
    fprintf(stderr,"Exiting to system.\n");
    exit(1);
}  /* End of: void erhand */

/**  Allocation of vector storage  ***********************************/
float *vector(n)
unsigned long n;
/* Allocates a float vector with range [1..n]. */
{
    float *v;

    v = (float *) malloc ((unsigned) n*sizeof(float));
    if (!v) erhand("Allocation failure in vector().");
    return v-1;
}  /* End of: float *vector */

/**  Allocation of float matrix storage  *****************************/
float **matrix(n,m)
unsigned long n, m;
/* Allocate a float matrix with range [1..n][1..m]. */
{
    unsigned long i;
    float **mat;

    /* Allocate pointers to rows. */
    mat = (float **) malloc((unsigned) (n)*sizeof(float*));
    if (!mat) erhand("Allocation failure 1 in matrix().");
    mat -= 1;

    /* Allocate rows and set pointers to them. */
    for (i = 1; i <= n; i++)
        {
        mat[i] = (float *) malloc((unsigned) (m)*sizeof(float));
        if (!mat[i]) erhand("Allocation failure 2 in matrix().");
        mat[i] -= 1;
        }

     /* Return pointer to array of pointers to rows. */
     return mat;
}  /* End of: float **matrix */

/**  Deallocate vector storage  *********************************/
void free_vector(v,n)
float *v;
unsigned long n;
/* Free a float vector allocated by vector(). */
{
   free((char*) (v+1));
}  /* End of: free_vector */

/**  Deallocate float matrix storage  ***************************/
void free_matrix(mat,n,m)
float **mat;
unsigned long n, m;
/* Free a float matrix allocated by matrix(). */
{
   unsigned long i;

   for (i = n; i >= 1; i--)
       {
       free ((char*) (mat[i]+1));
       }
   free ((char*) (mat+1));
} /* End of: free_matrix */

/**  Sort vector and also return sorted indexes ****************/
void piksr2(n, arr, brr)
unsigned long n;
float arr[], brr[];
/* Sorts an array arr[1..n] into ascending numerical order, 
   by straight insertion, while making the corresponding 
   rearrangement of the array brr[1..n] 
*/
{
   unsigned long i, j; 
   float a, b; 

   for (j=2; j<=n; j++) {   /* Pick out each element in turn */
       a=arr[j]; 
       b=brr[j]; 
       i=j-1;
       while (i > 0 && arr[i] > a) {    /* Look for the place to insert it */
             arr[i+1]=arr[i];
             brr[i+1]=brr[i];
             i--;
       }
       arr[i+1]=a;
       brr[i+1]=b;
   }
}  /* End of: piksr2 */

/**** Heapsort and also return sorted indexes *****************/
void sort2OLD(n, ra, rb)
     unsigned long n;
     float ra[], rb[]; 
/* Sorts an array ra[1..n] into ascending numerical order, 
   by heapsort, while making the corresponding 
   rearrangement of the array rb[1..n].  
   See Press et al., C, p. 247 */
{
   unsigned long i, j, l, ir; 
   float rra, rrb;

   l = (n >> 1)+1;
   ir = n;
   for (;;) {
      if (l > 1) {
	rra = ra[--l];
        rrb = rb[l];
      }
      else {
	rra = ra[ir];
	rrb = rb[ir];
	ra[ir] = ra[1];
	rb[ir] = rb[1];
	if (--ir == 1) {
	  ra[1] = rra; 
	  rb[1] = rrb; 
	  return;
	}
      }
      i = l;
      /* Should below be: j = l+1 << 1; ??? */ 
      j = l << 1;
      while (j <= ir) {
	if (j < ir && ra[j] < ra[j+1]) ++j;

	if (rra < ra[j]) {
	  ra[i] = ra[j];
	  rb[i] = rb[j];
	  j += (i=j); 
	}
	else  j = ir + 1;
      }
      ra[i] = rra;
      rb[i] = rrb;
   }
}  /* End of: sort2OLD */


void sort2(unsigned long n, float arr[], float brr[])
{
	unsigned long i,ir=n,j,k,l=1;
	int *istack,jstack=0;
	float a,b,temp;

	istack=ivector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
				arr[i+1]=a;
				brr[i+1]=b;
			}
			if (!jstack) {
				free_ivector(istack,1,NSTACK);
				return;
			}
			ir=istack[jstack];
			l=istack[jstack-1];
			jstack -= 2;
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			SWAP(brr[k],brr[l+1])
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
				SWAP(brr[l],brr[ir])
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
				SWAP(brr[l+1],brr[ir])
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1])
				SWAP(brr[l],brr[l+1])
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			b=brr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
				SWAP(brr[i],brr[j])
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			brr[l+1]=brr[j];
			brr[j]=b;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort2.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI

/*************** ivector ***********************************************/
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

/************** free_ivector ******************************************/
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

/************** nrerror ***********************************************/
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

