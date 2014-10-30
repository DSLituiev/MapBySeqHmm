/*******************************************************************
A MEX-adapted version of a C program "minfo.c" by Eric R. Weeks

Adaptation: Dmytro S. Lituiev, March 2013, University of Zurich, Switzerland

Handles a 2D arrays, with the signals running along the first dimension.

**************************************************
minfo.c -- Eric R. Weeks -- started 2/28/97

does the mutual information algorithm discussed by Fraser & Swinney
(Phys Rev A 33 (1986) p1134-1140)

v01:  2-28-97: taken from shell.c (7/1/96)
			quicksort routine taken from sane.c (original version)
v02:  2-28-97: revised sorting for s[] (different than q[])
			sped up math
v03:  2-28-97: add in tau loop
	 3-01-97: fix for variable number of input; add -b option
v04:  3-01-97: take out chi2 tests for substructure
v05:  3-01-97: realize that with chi2 tests taken out, number()
			function doesn't need to be called very often.  remove
			a[] and b[][] arrays!  Much faster now.

This program is public domain, although please leave my name and
email address attached.

email: weeks@physics.emory.edu
web: http://www.physics.emory.edu/~weeks/

explanation of how to use this program:
    http://www.physics.emory.edu/~weeks/software/minfo.html

 *******************************************************************/
#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MAXNUM 100000
#define KMAX 25
#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_3D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 3 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)

// mexCalcMutInfDelay
void mexFunction( mwIndex  nlhs, mxArray *plhs[],  mwIndex  nrhs, const mxArray *prhs[])
 {
 /* Macros for the ouput and input arguments */
    #define A_OUT plhs[0]
	#define E_IN prhs[0]
	#define T_IN prhs[1]
    mwIndex  Mx, Npx, n, k, m;

    // double *E, *T, *A;
    
	if(nrhs < 2 || nrhs > 3) /* Check the number of arguments */
	    mexErrMsgTxt("Wrong number of input arguments.");
	else if(nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");
		
     if (nlhs != 1) 
        mexErrMsgTxt( "One output argument required.");
    
    if (! (mxGetNumberOfDimensions(E_IN) == 2)){
        printf("number of dimensions: %u \n", mxGetNumberOfDimensions(E_IN)) ;
        mexErrMsgTxt("E is has other than 2 dims"); 
    }

	if (!IS_REAL_2D_FULL_DOUBLE(E_IN)) /* Check E */
		mexErrMsgTxt("E must be a real 2D full double array.");
    
    if (!IS_REAL_3D_FULL_DOUBLE(T_IN)) /* Check T */
		mexErrMsgTxt("T must be a real 3D full double array.");
	
    
	// get the input array dimensions and the pointer
	Mx = (mwIndex) mxGetM(E_IN); /* Get the dimensions */
	Npx = (mwIndex) mxGetN(E_IN);
    
   double *E = mxGetPr(E_IN);  /* Get the pointer to the data of E */   
    
   double *T = mxGetPr(T_IN);  /* Get the pointer to the data of T */

   mwSize A_SIZE[]={Npx, Npx, Mx};
   size_t nDimsA = 3;
   // mwSize *A_SIZE = mxGetDimensions(T_IN);
   // size_t nDimsA = mxGetNumberOfDimensions(T_IN);
   
   
  //  printf("dim[2]: %u \n", A_SIZE[2]);
 //  printf("trying to assign A \n ");

   A_OUT  = mxCreateNumericArray((size_t) 3, A_SIZE, mxDOUBLE_CLASS, mxREAL); /* Create the output matrix */
    
   double *A = mxGetPr(A_OUT); /* Get the pointer to the data of LALPHA */
   
 //  printf("A is assigned");
 	
 //  
	/*----- done reading in data -----*/

	

// mexPrintf("Nx:\t%u\tM:\t%u\n\n", Npx, Mx);

for(m=0;m<Mx;m++){
    for(n=0;n<Npx;n++){
      //   printf("m: %u, n: %u, i: %u \n", m, n, m + Mx* n);
        for(k=0; k<Npx; k++){
            // printf(" k: %u, n: %u, m: %u, i: %u  \n", k, n, m, n + Npx*m);
            A[n + Npx*k +  Npx*Npx*m] = T[n + Npx*k +  Npx*Npx*m]*E[m + Mx*n];             
          //  printf("i: %u, A[i]= %f \n", n + Npx*k +  Npx*Npx*m,  A[n + Npx*k +  Npx*Npx*m]);
        }
    }
}

return;
//	exit(0);
}	/* END OF MAIN */


