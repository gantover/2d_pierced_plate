#include <slepceps.h>
#include "interface_slepc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "time.h"

/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file was inspired from https://slepc.upv.es/documentation/current/src/eps/tutorials/ex1.c.html
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

/// @brief Solving with slepc the minimal eigen value problem
/// @param evals a pointer to a double to store the asked eigen value
/// @param evecs a pointer to an allocated space of n doubles to store the eigen vector
/// @return integer for error handling
int slepc(problem *s, double *evals, double *evecs)
{
    double ti, tf;
    PetscInt n = s->n;
    PetscInt nev = 1;

    Mat A;

    EPS eps;
    PetscInt its, nconv, i;
    PetscReal error;
    
    PetscScalar kr;
    Vec xr;

    PetscCall(SlepcInitialize(NULL,NULL,(char*)0,NULL));
    // NULL -> giving no external source of parameters, everything will be defined by calling functions

    /* Matrix initialisation */
    // since there is no matrix representation that takes ia,ja,a in petsc,
    // we have to copy every element to a new matrix in their format
    tic(mytimer_wall, ti);
    // nnz is an array containing the number of non-zero elements for each line of the matrix
    int *nnz = (int*)malloc(n*sizeof(int));
    for (int i = 0; i < n; i++) {
        nnz[i] = s->ia[i+1]-s->ia[i];
        // i did not find a way to make it faster with a buffer for ia[i+1]
    }
    PetscCall(MatCreateSeqAIJ(PETSC_COMM_WORLD, n, n, 0, nnz, &A));
    /* this is a sequential sparse matrix, sequential because we are not using a cluster of computers.
       preallocation of non-zeros with nnz allows us to be up to 3x faster to copy the matrix than simply 
       using MatCreate() and MatSetSizes() and let petsc do the work blindly */

    PetscCall(MatSetFromOptions(A));
    PetscCall(MatSetUp(A));
    for (int i = 0; i < n; i++) {
        for (int j = s->ia[i]; j < s->ia[i+1]; j++) {
            PetscCall(MatSetValue(A,i,s->ja[j],s->a[j],INSERT_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));
    free(nnz);
    tac(mytimer_wall, tf, "slepc matrix copy");

    /* vector allocation (real and imaginary part) */
    PetscCall(MatCreateVecs(A,NULL,&xr));

    /* Solver parameters */
    PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));
    PetscCall(EPSSetOperators(eps,A,NULL)); // A.x = lambda.B.x, B=NULL
    PetscCall(EPSSetType(eps, EPSBLOPEX));
    // this is the solution method, here we use blopex : specialized for minimal eigenvalue retrieval 
    PetscCall(EPSSetProblemType(eps,EPS_HEP)); // Hermitian eigenvalue problem

    PetscCall(EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL));
    // in the case of blopex, we coult omit this line since we are already searching for the minimal eigenvalue

    PetscCall(EPSSetFromOptions(eps));
    PetscCall(EPSSetDimensions(eps,nev,PETSC_DEFAULT,PETSC_DEFAULT));

    tictac(PetscCall(EPSSolve(eps)), "slepc to solve", mytimer_wall, ti, tf);
    
    #if SLEPC_CONFIG_PRINT
    vspace;
    broadcast("SLEPC Solver Parameters");
    EPSView(eps, NULL);
    vspace;
    broadcast("SLEPC Results");
    #endif

    PetscCall(EPSGetIterationNumber(eps,&its));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %" PetscInt_FMT "\n",its));
    PetscCall(EPSGetConverged(eps,&nconv));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %" PetscInt_FMT "\n",nconv));


    if (nconv>0) {
        /* Display eigenvalues and relative errors */
        for (i=0;i<nconv;i++) {
            /* Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) */ 
            PetscCall(EPSGetEigenpair(eps,i,&kr,NULL,xr,NULL));
            // NULL replaces both ki and xi, imaginary part of the eigen value and the eigen vector
            // They should always be = 0 because we have a symmetric real matrix

            /* Compute the relative error associated to each eigenpair */
            PetscCall(EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error));

            PetscCall(PetscPrintf(PETSC_COMM_WORLD," eigenvalue calculated : %12e, error : %12g\n",(double)kr,(double)error));
        }
    } else {
        printf("0 converged eigenvalues, closing program\n");
        return EXIT_FAILURE;
    }

    /* transfer eigenvalue data outside slepc */
    evals[0] = kr;

    /* transfer vector data outside slepc */
    PetscScalar *temp;
    PetscCall(VecGetArray(xr, &temp));
    // we get a pointer to a contiguous form of the vector array
    for (int i = 0; i < n; i++) {
        // by the method used and the check for nconv != 0, we know that nconv=1 and that we should not go beyond n
        evecs[i] = temp[i];
    }
    PetscCall(VecRestoreArray(xr, &temp));

    /* free heap memory */
    PetscCall(EPSDestroy(&eps));
    PetscCall(MatDestroy(&A));
    PetscCall(VecDestroy(&xr));
    PetscCall(SlepcFinalize());
    return EXIT_SUCCESS;
}