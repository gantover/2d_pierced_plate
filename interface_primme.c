#include <stdlib.h>
#include "interface_primme.h"
#include "config.h"
#include "time.h"

static double *a;
static int n, *ia, *ja;

/// @brief to initialize static varibles for primme
/// @param primme_n number of unknowns in the system
/// @param primme_ia array 'ia' of matrix A
/// @param primme_ja array 'ja' of matrix A
/// @param primme_a array 'a' of matrix A
/// @return integer for error handling
int init_primme(int primme_n, int *primme_ia, int *primme_ja, double *primme_a) 
{
    n = primme_n;
    a = primme_a;
    ja = primme_ja;
    ia = primme_ia;
    return EXIT_SUCCESS;
}

/// @brief Calculate the matrix-vector product vy = A*vx.
/// The A matrix has to be stored beforehand in static variables (n,ia,ja,a) corresponding to CSR format
/// @param vx input vector(s)
/// @param vy output vector(s) from A*vx
/// @param blockSize number of vectors
/// @param primme input parameters to optimise calculation when using with primme.
/// @note primme can be set to NULL
void matvec_primme(void *vx, void *vy, int *blockSize, primme_params *primme)
{
    int i, j, b;
    double *x = (double*)vx, *y=(double*)vy;

    for(b = 0; b < (*blockSize)*n; b+=n)
        for(i = 0; i < n; i++){
            y[b+i] = 0;
            for (j = ia[i]; j < ia[i + 1]; j++)
                y[b+i] += a[j] * x[b+ja[j]];
        }
} 


/// @brief Calculate the lowest and highest eigen value of the matrix A of dimenssions primme_n x primme_n.
/// Stored in the CSR format with the help of primme_ia, primme_ja, primme_a vectors
/// @param min_evals minimal eigen value
/// @param min_evecs eigen vector from minimal eigen value
/// @param max_evals maximal eigen value
/// @param max_evecs eigen vector from maximal eigen value
/// @return integer for error handling
int primme(double *min_evals, double *min_evecs, double *max_evals, double *max_evecs)
{
    /*  note : compared to the original version of this program,
        primme.numEvals = nev has disappeared since primme_largest and primme_lowest already asks for one eigenvalue 
        this can be verified with PRIMME_CONFIG_PRINT 1 in config.h */

    double ti, tf;
    int err;

    // residual norm buffer
    double *resn = (double*)malloc(sizeof(double));
    if (resn == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour un vecteur auxilier dans la fonction primme\n\n");
        return 1;
    }

    /* Min eigenvalue */
    primme_params primme;
    primme_initialize (&primme);
    primme.matrixMatvec = matvec_primme; 
    primme.target = primme_smallest;
    primme.n = n;
    primme.printLevel = 0; // we want to handle the results output ourselves
    if((err = primme_set_method (DEFAULT_MIN_TIME, &primme))) {
        printf("\nPRIMME: erreur N %d dans le choix de la methode \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }
  
    /* display PRIMME parameters, they will essentialy be the same for the maximum eigenvalue problem */
    #if PRIMME_CONFIG_PRINT
    primme_display_params (primme);
    broadcast("primme results")
    #endif /* PRIMME_PRINT */

    tic(mytimer_wall, ti);
    if((err = dprimme (min_evals, min_evecs, resn, &primme))) {
        printf("\nPRIMME: erreur N %d dans le calcul des valeurs propres \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return EXIT_FAILURE;
    }
    tac(mytimer_wall, tf, "primme to solve for mininal eigenvalue");

    printf("Minimal eigen value: %e, error : %e\n", min_evals[0], resn[0]);

    primme_Free (&primme);


    /* Max eigenvalue */
    /* /!\ we have to reconfigure primme since dprimme de-initialized some parameters */
    primme_initialize (&primme);
    primme.matrixMatvec = matvec_primme; 
    primme.n = n;
    primme.printLevel = 0; 
    primme.target = primme_largest; 
    if((err = primme_set_method (DEFAULT_MIN_TIME, &primme))) {
        printf("\nPRIMME: erreur N %d dans le choix de la methode \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }

    tic(mytimer_wall, ti);
    if((err = dprimme (max_evals, max_evecs, resn, &primme))) {
        printf("\nPRIMME: erreur N %d dans le calcul des valeurs propres \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }
    tac(mytimer_wall, tf, "primme to solve for maximal eigenvalue");

    printf("Maximal eigen value : %e, error : %e\n", max_evals[0], resn[0]);

    /* free memory */
    primme_Free (&primme); free(resn);

    return EXIT_SUCCESS;
}
