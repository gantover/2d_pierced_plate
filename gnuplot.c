#include "gnuplot.h"
#include "prob.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NEW_LINE fprintf(s->context, "\n")
#define ZERO fprintf(s->context, "%g ", 0.0)
#define WVAL(x,y,z) fprintf(f, "%g %g %g\n", (double)x, (double)y, z);

/// @brief this function write the data of v line by line in a gnuplot buffer, 
///        this data will have to be displayed right after with the open() function
/// @param s the gnuplot object handling its process
/// @param p the problem object containing its shape,...
/// @param v a vecctor holding the value for every point of the grid starting from the bottom left corner
/// @return integer for error handling
/// @attention This function also writes a boundary fixed at the value of 0
int write_for_gnuplot(gnuplot *s, problem *p, double *v) 
{
    int nx = p->nx;
    int ny = p->ny;

    Rectangle is = p->i_s;

    FILE *f = s->context;

    int ind = 0;

    // first horizontal boundary 
    for (int ix = 0; ix < nx + 2; ix ++) {
        WVAL(ix, 0, 0.0);
    }
    NEW_LINE;

    for (int iy = 1; iy < is.y[0]+1; iy ++) {
        // vertical boundary
        WVAL(0, iy, 0.0);
        // line data
        for (int ix = 1; ix < nx+1; ix ++) {
            WVAL(ix, iy, v[ind]);
            ind ++;
        }
        // vertical boundary
        WVAL(nx + 1, iy, 0.0);
        NEW_LINE;
    }

    for (int iy = is.y[0]+1; iy <= is.y[1]+1; iy ++) {
        // vertical boundary
        WVAL(0, iy, 0.0); 
        // line data
        for (int ix = 1; ix < is.x[0]+1; ix++) {
            WVAL(ix, iy, v[ind]);
            ind ++;
        }
        // the hole
        for (int ix = is.x[0]+1; ix <= is.x[1]+1; ix++) {
            WVAL(ix, iy, 0.0);
        }
        // line data
        for (int ix = is.x[1]+2; ix < nx+1; ix++) {
            WVAL(ix, iy, v[ind]);
            ind ++;
        }
        // vertical boundary
        WVAL(nx+1, iy, 0.0);
        NEW_LINE;
    }
    for (int iy = is.y[1]+2; iy < ny+1; iy++) {
        // vertical boundary
        WVAL(0, iy, 0.0);
        // line data
        for (int ix = 1; ix < nx+1; ix ++) {
            WVAL(ix, iy, v[ind]);
            ind ++;
        }
        // vertical boundary
        WVAL(nx+1, iy, 0.0); // ix = 0
        NEW_LINE;
    }

    // second horizontal boundary 
    for (int ix = 0; ix < nx + 2; ix ++) {
        WVAL(ix, ny+1, 0.0);
    }
    
    fprintf(s->context, "e\n"); 

    return EXIT_SUCCESS;
}

/// @brief This is a slower but shorter (in terms of lines of code) version of write_for_gnuplot()
int write_for_gnuplot_slow(gnuplot *s, problem *p, double *v) 
{
    int nx = p->nx; int ny = p->ny;

    FILE *f = s->context;

    int ind = 0;
    int *inds = p->inds;

    for (int ix = 0; ix < nx + 2; ix++) WVAL(ix, 0, 0.0);

    NEW_LINE;

    for (int iy = 1; iy < ny+1; iy++) {
        WVAL(0, iy, 0.0);
        for (int ix = 1; ix < nx+1; ix++) {
            if (inds[ind] != -1) {
                WVAL(ix, iy, v[inds[ind]]);
            } else {
                WVAL(ix, iy, 0.0);
            }
            ind ++;
        }
        WVAL(nx + 1, iy, 0.0);
        NEW_LINE;
    }

    for (int ix = 0; ix < nx + 2; ix++) WVAL(ix, ny+1, 0.0);
    fprintf(s->context, "e\n"); 

    return EXIT_SUCCESS;
}

/// @brief After writing to gnuplot the gnuplot context, 
///        we ask the gnuplot process to show that data on the screen
/// @param s the gnuplot object handling its process
/// @param p the problem object containing its shape,...
/// @param title the title to display on top of the plot
/// @return integer for error handling
int open_gnuplot_data(gnuplot *s, problem *p, char* title) 
/*
    */
{
    fprintf(s->context, "set title '%s'\n", title);
    fprintf(s->context, "%s\n", s->plot_cmd);
    return EXIT_SUCCESS;
}

/// @brief handles the exit off gnuplot by carefully closing the pipe
/// @param s the gnuplot object handling its process
/// @return integer for error handling
int close_gnuplot_context(gnuplot *s) 
{
    fprintf(s->context, "unset output; exit gnuplot\n");
    free(s->plot_cmd);
    fclose(s->context);
    return EXIT_SUCCESS;
}

/// @brief this function initializes the gnuplot object
/// @param self the yet uninitialized gnuplot object
/// @param config configures the gnuplot context
/// @param plot_cmd is the command that will be executed everytime we want to open the data
/// @param p the problem object containing its shape,...
/// @return integer for error handling 
int init_gnuplot(gnuplot* self, char* config, char* plot_cmd, problem* p) 
{
    FILE *context = popen("gnuplot -persist", "w");
    if (context == NULL) {
        return EXIT_FAILURE;
    }

    double ratio = (double)p->ny / p->nx;

    fprintf(context, "set size ratio %f\n", ratio);
    fprintf(context, "set view equal xy\n");
    fprintf(context, "set nokey\n");

    fprintf(context, "%s\n", config);

    self->context = context;

    self->plot_cmd = (char*)malloc((strlen(plot_cmd)+1) * sizeof(char));
    strcpy(self->plot_cmd, plot_cmd);

    self->open = open_gnuplot_data;
    self->close = close_gnuplot_context;
    self->write = write_for_gnuplot;
    return EXIT_SUCCESS;
}
