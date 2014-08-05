
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
#define MAX_M 100000
#define KMAX 25
#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_3D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 3 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == return array1)


// mexCalcMutInfDelay
void mexFunction( mwIndex  nlhs, mxArray *plhs[],  mwIndex  nrhs, const mxArray *prhs[])
 {
 /* Macros for the ouput and input arguments */
    #define LALPHA_OUT plhs[0]
    #define LBETA_OUT plhs[1]
    #define  A_OUT plhs[2]
	#define E_IN prhs[0]
	#define T_IN prhs[1]
    mwIndex  Mx, Npx, n, k, m;

    double INF = mxGetInf();
    
    // double *E, *T, *A;
    
	if(nrhs < 2 || nrhs > 3) /* Check the number of arguments */
	    mexErrMsgTxt("Wrong number of input arguments.");
	
    if(nlhs >3) 
		mexErrMsgTxt("Too many output arguments.");
		
     if (nlhs != 3) 
        mexErrMsgTxt( "Two output argument required.");
    
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


   // const size_t nDimsA = 3;
   // mwSize *A_SIZE = mxGetDimensions(T_IN);
   A_OUT = mxCreateNumericArray( (size_t) 3,  mxGetDimensions(T_IN), mxDOUBLE_CLASS, mxREAL); /* Create an intermediate matrix */
   double *A = mxGetPr(A_OUT); /* Get the pointer to the data of LALPHA */
  // double A[Mx*Npx*Npx];
           
    mwSize Out_SIZE[]={Mx, Npx};
    
	LALPHA_OUT = mxCreateDoubleMatrix( Mx, Npx, mxREAL); /* Create the output matrix */
	double *LALPHA = mxGetPr(LALPHA_OUT); /* Get the pointer to the data of LALPHA */

 	LBETA_OUT = mxCreateDoubleMatrix( Mx, Npx, mxREAL);/* Create the output matrix */
	double *LBETA = mxGetPr(LBETA_OUT); /* Get the pointer to the data of LALPHA */
    
	
	/*----- done reading in data -----*/
	
// mexPrintf("Nx:\t%u\tM:\t%u\n\n", Npx, Mx);

//--------------- A
for(m=0;m<(Mx-1);m++){
    for(n=0;n<Npx;n++){
      //   printf("m: %u, n: %u, i: %u \n", m, n, m + Mx* n);
        for(k=0; k<Npx; k++){
            A[n + Npx*k +  Npx*Npx*m] = T[n + Npx*k +  Npx*Npx*m]*E[m + Mx*n];             
           //   printf("i: %u, A[i]= %f \n", n + Npx*k +  Npx*Npx*m,  A[n + Npx*k +  Npx*Npx*m]);
        }
    }
}

// printf("E*T matrix has been calculated \n\n");
//--------------- Alpha
double aCum[Mx*Npx];
    
for(n=0;n<Npx;n++) {
    aCum[Mx-1 + Mx*n] = E[Mx-1 + Mx*n] ; 
    LALPHA[Mx-1 + Mx*n] = log10( aCum[Mx-1 + Mx*n] ); 
    
      //      printf(" i: %u, E: %f  \n",  Mx-1 + Mx*n, aCum[Mx-1 + Mx*n]  );
    }

    double scalePrev = 0.0;    
    double scaleNew = 0.0;
    double expScaleDiff = 1.0;
    double scaleMax = -DBL_MAX;

            
for(m=(Mx-2);m>=0;m--){
    expScaleDiff = pow(10, scaleNew - scalePrev);
    for(n=0;n<Npx;n++){
        aCum[m + Mx*n] = 0.0;
        for(k=0; k<Npx; k++){
         //   printf(" i: %u, aCum: %f  \n",  m + 1 + Mx*n, aCum[m + 1 + Mx*k]);
            aCum[m + Mx*n] += A[n + Npx*k +  Npx*Npx*m] * aCum[(m + 1) + Mx*k] * expScaleDiff;    
            LALPHA[m + Mx*n]= log10(aCum[m + Mx*n]) - scaleNew;
        }
        if (scaleMax < LALPHA[m + Mx*n] + scaleNew)
            scaleMax = LALPHA[m + Mx*n] + scaleNew;
    }
             scalePrev = scaleNew;
             scaleNew = - scaleMax;
         //  printf("m: %u, scaleMax: %f\n", m, scaleMax);
}

//--------------- Beta
    scalePrev = 0.0;    
    scaleNew = 0.0;
    scaleMax = -DBL_MAX;

for(n=0;n<Npx;n++) {
    aCum[ 0 + Mx*n] = 1.0; 
    LBETA[ 0 + Mx*n] = 0.0; 
    // printf(" i: %u, E: %f  \n",  Mx-1 + Mx*n, aCum[Mx-1 + Mx*n]  );
    }
    
for(m=1; m<Mx; m++){
   //  printf(" m: %u  =============== \n", m);
    printf("m: %u, scaleMax: %f\n", m, scalePrev);
    expScaleDiff = pow(10, scaleNew - scalePrev);
    for(n=0; n<Npx; n++){
        aCum[m + Mx*n] = 0.0;
        for(k=0; k<Npx; k++){           
         //   printf(" i: %u, aCum[m-1,k]: %f  \n",  m - 1 + Mx*k, aCum[m - 1 + Mx*k]);            
         //   printf(" i: %u, A[n + Npx*k +  Npx*Npx*m]: %f  \n",   n + Npx*k +  Npx*Npx*(m-1),  A[n + Npx*k +  Npx*Npx*(m-1)]);
            aCum[m + Mx*n] += A[n + Npx*k +  Npx*Npx*(m-1)] * aCum[(m - 1) + Mx*k] * expScaleDiff;  
            LBETA[m + Mx*n]= log10(aCum[m + Mx*n]) - scaleNew;
        }
      //  printf(" i: %u, aCum[m,n]: %f  \n",  m + Mx*n, aCum[m + Mx*n]);
        if (scaleMax < LBETA[m + Mx*n] + scaleNew)
            scaleMax = LBETA[m + Mx*n] + scaleNew;
    }
            scalePrev = scaleNew;
            scaleNew = - scaleMax;
}
    

   //LALPHA[m + Mx*n] = -INF;
return;
//	exit(0);
}	/* END OF MAIN */



