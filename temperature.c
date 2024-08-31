#include "temperature.h"
#include "gnuplot.h"
#include <stdlib.h>
#include "config.h"
double d = DIFFUSIVITY;

/// @brief iterates through u(k+1) = (I-DA)u(k) = u(k) -Dv(k)
///       to get the temperature evolution using the progressive euler method.
///       At the end u(k+1) is returned through uk.
/// @param uk temperature at any point of the grid
/// @param vk result of A*uk
/// @param n size of both uk and vk vector
/// @param dt time step of the progressive euler method
/// @param t time elapsed since starting the method 
void temperature_iterate(double *uk, double *vk, int n, double dt, double *t) {
    for (int i = 0; i < n; i++) {
        uk[i] -= dt*d*vk[i];
    }
    (*t)+=dt;
}
