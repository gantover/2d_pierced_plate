#include "prob.h"
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "interface_primme.h"
#define square(x) (x)*(x)

/// @brief initializes a rectangle object
/// @param x1 left 
/// @param x2 right 
/// @param y1 bottom 
/// @param y2 upper
void init_rectangle(Rectangle* self, int x1, int x2, int y1, int y2) {
    self->x[0] = x1;
    self->x[1] = x2;
    self->y[0] = y1;
    self->y[1] = y2;
}

/// @brief gives the coordinates in the grid system (defined by m)
/// @param sub_shape contains the coordinates in the metric system
/// @param m the number of grid points for the unit lenght
/// @return the Rectangle with adapted coordinates
Rectangle get_sub_shape_indices(Rectangle *sub_shape, int m) {
    Rectangle final;
    final.x[0] = sub_shape->x[0] * (m-1) - 1;
    final.x[1] = sub_shape->x[1] * (m-1) - 1;
    final.y[0] = sub_shape->y[0] * (m-1) - 1;
    final.y[1] = sub_shape->y[1] * (m-1) - 1;
    return final;
}

/// @brief takes a position on the grid and returns if it is within (not strictly)
/// the boundaries of the hole
/// @param is the hole in the grid coordinates system
/// @param ix the point the in grid coordinates system
/// @param iy the point the in grid coordinates system
/// @return true if in the zone
int in_zone(Rectangle *is, int ix, int iy) {
    if ((ix >= is->x[0]) && (ix <= is->x[1]) && (iy >= is->y[0]) && (iy <= is->y[1]))
        return true; 
    return false;
}

/// @brief generates the problem matrix with the help of an array to hold offset indices
/// to take into account the wall
/// @return integer for error handling
int generate_mat(problem *s) {
    int ind = 0;
    int u = 0;
    int *inds = s->inds;
    double invh2 = (s->m-1)*(s->m-1); // for unit lenght
    
    for (int iy = 0; iy < s->ny; iy++) {
        for (int ix = 0; ix < s->nx; ix++) {
            if (in_zone(&s->i_s,ix,iy)) {
                inds[u] = -1;
                // this will later indicate that we are in the wall
            } else {
                inds[u] = ind;
                ind ++;
            }
            u++;
        }
    }

    /*
        the inds array makes it easy to keep the function of the problem without an hole. 
        for instance, to have the south neighbor, instead of taking ind-nx, i convert it with
        inds[ind-nx] which gives me the right place in the matrix taking the hole into account
    */

    int nnz = 0;
    ind = 0;
    for (int iy = 0; iy < s->ny; iy++) {
        for (int ix = 0; ix < s->nx; ix++) {

            ind = ix + s->nx * iy;
            if (inds[ind] == -1) continue;
            s->ia[inds[ind]] = nnz;

            /* filling up the line : south neighbor */
            if (iy > 0 && (inds[ind - s->nx] != -1))  {
                s->a[nnz] = -invh2; /* for D=1 */
                s->ja[nnz] = inds[ind - s->nx];
                nnz++; 
            }

            /* filling up the line : west neighbor */
            if (ix > 0 && (inds[ind - 1] != -1))  {
                s->a[nnz] = -invh2; /* for D=1 */
                s->ja[nnz] = inds[ind - 1];
                nnz++;
            }
            /* filling up the line : diagonal element */
            s->a[nnz] = 4.0*invh2; /* for D=1 */
            s->ja[nnz] = inds[ind];
            nnz++;

            /* filling up the line : east neighbor */
            if (ix < s->nx - 1 && (inds[ind + 1] != -1)) {
                s->a[nnz] = -invh2; /* for D=1 */
                s->ja[nnz] = inds[ind + 1];
                nnz++;
            }

            /* filling up the line : north neighbor */
            if (iy < s->ny - 1 && (inds[ind + s->nx] != -1)) {
                s->a[nnz] = -invh2; /* for D=1 */
                s->ja[nnz] = inds[ind + s->nx];
                nnz++;
            }
        }
    }
    (s->ia)[inds[ind]+1] = nnz; 
    // we give a final ia element to know the number of elemnts of the last column

    return EXIT_SUCCESS;
}

/// @brief Calculates ||u-v||/||u|| using the euclidian norm
/// @param u should be the vector found by primme
/// @param v the vector to compare with
/// @param n the size of both vectors
/// @return result of calculation
/// @attention the absolute value of the vector elements are taken
double compare_vecs(double *u, double *v, int n) 
{
    double umv_buf = 0;
    double u_buf = 0;
    for (int i = 0; i < n; i ++) {
        umv_buf += square(abs(u[i])-abs(v[i]));
        // the absolute value has to be taken
        // sometimes those vectors are the same but opposed
        u_buf += square(u[i]);
    }
    return sqrt(umv_buf/u_buf);
}

/// @brief Calculing the norm of the residual from A u = w2 u with the formula ||A u - w2 u||/||u|| using euclidian norm
/// @param u the calculated eigen vec 
/// @param w2 the calculated eigen value
/// @return the norm of the residual
double calc_res(problem *s, double *u, double w2) {
    double result = 0;
    double u_norm2 = 0;
    int ia_buf = 0; // ia[0] == 0 is always true

    double line_result;
    for (int i = 0; i < s->n; i++) {
        line_result = -w2*u[i];
        u_norm2 += square(u[i]);
        for (int j = ia_buf; j < (ia_buf= s->ia[i+1]); j++) {
            line_result += s->a[j] * u[s->ja[j]];
        }
        result += square(line_result);
    }
    result /= u_norm2;
    return sqrt(result);
}

/// @brief Creates 3 files to extract the CSR matrix
/// @return integer for error handling
int extract_mat(problem *s) {
    line_sep;
    FILE *pia, *pja, *pa;
    pia = fopen("./compare_mat/ia_gen.txt","w");
    pja = fopen("./compare_mat/ja_gen.txt","w");
    pa = fopen("./compare_mat/a_gen.txt","w");
    if(pia == NULL || pja == NULL || pa == NULL) return EXIT_FAILURE;
    int nnz = s->ia[s->n]-s->ia[0];
    for (int i = 0; i < s->n + 1; i++) {
        fprintf(pia,"%i\n", s->ia[i]);
    }
    for (int i = 0; i < nnz; i++) {
        fprintf(pja, "%i\n", s->ja[i]);
    }
    for (int i = 0; i < nnz; i++) {
        fprintf(pa, "%.15lf\n", s->a[i]);
    }
    fclose(pia); fclose(pja); fclose(pa);
    return EXIT_SUCCESS;
}

/// @brief frees all heap allocated memory for problem object
void remove_problem(problem *s) {
    free(s->ia);
    free(s->ja);
    free(s->a);
    free(s->inds);
}

/// @brief Initializes the problem object
/// @param self The yet unitialized object
/// @param m The number of grid point for the unit lenght
/// @param shape The shape of the membrane
/// @param sub_shape The shape of the hole
/// @return integer for error handling
/// @attention This program will not work for holes that overlap the boundary of the membrane,
/// they should be stricly contained inside
int init_problem(problem *self, int m, pos2d shape, Rectangle sub_shape) {
    /* struct for storage of problem data, makes
     passing problem data in function arguments easier */

    self->m = m;

    // main shape
    self->m_s = shape;
    self->nx = shape.x * (m-1) - 1;
    self->ny = shape.y * (m-1) - 1;

    // sub shape (hole)
    self->s_s = sub_shape;
    self->i_s = get_sub_shape_indices(&self->s_s, m);
    int nx_is = (self->s_s.x[1]-self->s_s.x[0])*(m-1)+1; 
    int ny_is = (self->s_s.y[1]-self->s_s.y[0])*(m-1)+1; 
    self->nx_is = nx_is;
    self->ny_is = ny_is;

    /* allocations */
    self->inds = (int*)malloc(sizeof(int) * self->nx*self->ny);
    self->n = self->nx * self->ny - nx_is * ny_is;
    self->ia = (int*)malloc(((self->n)+1) * sizeof(int));
    // number of non-zero elements
    int nnz = 5 * self->nx * self->ny- 2*self->nx - 2*self->ny - 2*nx_is -2*ny_is;
    self->ja = (int*)malloc(nnz * sizeof(int));
    self->a = (double*)malloc(nnz * sizeof(double));
    if (self->inds == NULL || self->ia == NULL || self->ja == NULL || self->a == NULL ) {
        printf("\n ERROR : not enough memory to generate the matrix\n\n");
        return EXIT_FAILURE;
    }
    
    // function pointers
    self->generate_mat = generate_mat;
    self->close = remove_problem;
    self->extract_mat = extract_mat;

    return EXIT_SUCCESS;
}
