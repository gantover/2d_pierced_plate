#ifndef H_PROB
#define H_PROB

#include <stdlib.h>
#include <stdio.h>
#define line_sep printf("\n--------------------------------\n")

typedef struct {
    int x,y;
} pos2d;

typedef struct sRectangle Rectangle;
struct sRectangle {
    int x[2];
    int y[2];
};

void init_rectangle(Rectangle* self, int x1, int x2, int y1, int y2);

Rectangle get_sub_shape_indices(Rectangle *sub_shape, int m);

typedef struct sProblem problem;

struct sProblem {
    pos2d m_s; // main shape
    Rectangle s_s; // sub shape
    Rectangle i_s; // sub shape (index data)
    int *ia, *ja;
    double *a;
    int *inds; // indices for each point the the grid, -1 to indicate the hole
    int m, n, nnz, nx, ny;
    int nx_is, ny_is;
    int (*generate_mat)(problem*);
    void (*close)(problem*);
    int (*extract_mat)(problem*);
};

int init_problem(problem* self, int m, pos2d shape, Rectangle sub_shape);

double calc_res(problem *s, double *u, double w2);
double compare_vecs(double* u, double* v, int n);

#endif // !H_PROB
