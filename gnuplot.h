#ifndef H_GNUPLOT

#define H_GNUPLOT

#include "prob.h"
#include <stdio.h>
#include <stdlib.h>

typedef struct sGnuplot gnuplot;
struct sGnuplot {
    FILE *context;
    char *plot_cmd;
    problem *p;
    int (*write)(gnuplot*, problem*, double*);
    int (*open)(gnuplot*, problem*, char*);
    int (*close)(gnuplot*);
};

int init_gnuplot(gnuplot* self, char* config, char* plot_cmd, problem* p);

#endif // !H_GNUPLOT
