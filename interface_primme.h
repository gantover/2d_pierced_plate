#ifndef INTERFACE_PRIMME_H
#define INTERFACE_PRIMME_H

#include "./primme/PRIMMESRC/COMMONSRC/primme.h"

int init_primme(int primme_n, int *primme_ia, int *primme_ja, double *primme_a);

int primme(double *min_evals, double *min_evecs, double *max_evals, double *max_evecs);

void matvec_primme(void *vx, void *vy, int *blockSize, primme_params *primme);

#endif /* INTERFACE_PRIMME_H */