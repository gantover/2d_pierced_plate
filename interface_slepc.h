#ifndef INTERFACE_SLEPC_H
#define INTERFACE_SLEPC_H

#include "prob.h"

int slepc(problem *s, double *evals, double *evecs);

#endif // !INTERFACE_SLEPC_H