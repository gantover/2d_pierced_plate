#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <signal.h>
#include "prob.h"
#include "time.h"
#include "interface_primme.h"
#include "interface_slepc.h"
#include "gnuplot.h" 
#include "temperature.h"
#include "config.h"

static volatile bool running = true;

/// @brief this function will handle a SIGTERM sent from wrapper.sh and change the value of running to false.
/// When running goes to false, the temperature evolution loop stops and makes a goto: to the safe closing of the program
/// @param dummy the signal function needs a function that accepts an int but we don't need it here
void gnuplot_loop_handler(int dummy) {
  running = false;
}

int main(int argc, char *argv[])
{
  vspace;

  int m = M_UNIT_STEPS;
  double *min_evals, *max_evals, *min_evecs, *max_evecs;
  double ti,tf;

  pos2d shape = {4,5}; // size of membrane
  Rectangle sub_shape; init_rectangle(&sub_shape, 1, 2, 1, 3); // size of hole
  /* /!\ this program only support holes that don't overlap the membrane edges */
  broadcast("Problem Initialisation")

  problem p; if (init_problem(&p, m, shape, sub_shape)) return EXIT_FAILURE;

  tictac(p.generate_mat(&p),"generating problem matrix", mytimer_wall, ti, tf);

  #if EXTRACT_MAT
  // extracting the problem matrix to 3 files in ./compare_mat
  if(p.extract_mat(&p)) {printf("failed to open files to extract matrix\n"); return EXIT_FAILURE;}
  #endif

  printf("m = %5d   n = %8d  nnz = %9d\n", m, p.n, p.ia[p.n] );
  vspace;
  
  /* allocate memory for vectors & eigenvalues */
  min_evals = (double*)malloc(sizeof(double));
  max_evals = (double*)malloc(sizeof(double));
  min_evecs = (double*)malloc(p.n * sizeof(double));
  max_evecs = (double*)malloc(p.n * sizeof(double));

  if (min_evals == NULL || max_evals == NULL || min_evecs == NULL || max_evecs == NULL) {
      printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs propres\n\n");
      return EXIT_FAILURE;
  }

  #if SOLVING_WITH_SLEPC
  double *slepc_evals, *slepc_evecs;
  slepc_evals = (double*)malloc(sizeof(double));
  slepc_evecs = (double*)malloc(p.n * sizeof(double));
  if (slepc_evals == NULL || slepc_evecs == NULL) {
      printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs propres\n\n");
      return EXIT_FAILURE;
  }
  #endif

  /* primme solver */
  broadcast("solving with primme");
  init_primme(p.n, p.ia, p.ja, p.a);
  if(primme(min_evals, min_evecs, max_evals, max_evecs))
     return EXIT_FAILURE;
  vspace;

  /* alternative solver : slepc with blopex */
  #if SOLVING_WITH_SLEPC
  broadcast("solving with slepc for minimal eigenvalue");
  if (slepc(&p, slepc_evals, slepc_evecs)) printf("slepc failed\n");
  vspace;
  broadcast("comparing eigenvectors and eigenvalues from primme and slepc")
  double compare_vectors = compare_vecs(min_evecs, slepc_evecs, p.n);
  printf("||primme-slepc||/||primme|| = %e\n", compare_vectors);
  double compare_values = compare_vecs(min_evals, slepc_evals, 1);
  printf("|primme-slepc|/|primme| = %e\n", compare_values);
  vspace;
  #endif

  /* residual calculation */
  #if CALC_RESIDUAL
  broadcast("calculating the residual for primme min eigenvalue");
  tic(mytimer_wall, ti);
  double res = calc_res(&p, min_evecs, min_evals[0]);
  tac(mytimer_wall, tf, "it");
  printf("calculated residual was : %e\n", res);
  vspace;
  #endif

  /* minimal eigenvector display */
  #if DISPLAY_EIGENVEC
  broadcast("opening gnuplot context for minimal eigenvector")
  printf(ANSI_COLOR_RED "/!\\ you should close the gnuplot window in order to continue...\n" ANSI_COLOR_RESET);
  char mp_title[] = "eigen vector for min eigenvalue";
  char mp_plotcmd[] = "splot '-' using 1:2:3";
  char mp_config[] = "set pm3d; set hidden3d\npause mouse keypress;\n";
  gnuplot mp; init_gnuplot(&mp, mp_config, mp_plotcmd,&p);
  mp.open(&mp, &p, mp_title);
  tictac(mp.write(&mp, &p, min_evecs),"writing to gnuplot", mytimer_wall, ti,tf);
  mp.close(&mp);
  vspace;
  #endif /*DISPLAY_EIGENVEC*/
  
  #if SHOW_TEMPERATURE_EVOL
  broadcast("showing heat evolution")
  char hp_plotcmd[] = "plot '-' using 1:2:3 with image";
  char hp_config[] = "set palette rgb 33,13,10\nset cbrange [0:10]";
  gnuplot hp; init_gnuplot(&hp, hp_config, hp_plotcmd, &p);

  double *uk = (double*)malloc(sizeof(double) * p.n);
  for (int i = 0; i < p.n; i++) {
    uk[i] = INITIAL_TEMP;
  }
  char title[64]; // will be used to display the time on top of the graph
  double *vk = (double*)malloc(sizeof(double) * p.n); // array to store the product A*uk
  double t = 0;

  signal(SIGTERM, gnuplot_loop_handler);
  /* sigterm can be send by the wrapper program or by any task manger */

  int blockSize = 1; // needed for matvec_primme to work

  double dt_max = 2.0 / (max_evals[0]*DIFFUSIVITY);
  double dt = dt_max / DT_F;
  int tt = TOTAL_TIME;
  int isr = ISR; 
  double tti = tt / dt;
  int out_loop_tt = ceil(tti/isr);
  printf("dt max for guaranteed convergence of time evolution : %fs\n", dt_max);

  for (int i = 0; i < out_loop_tt; i++) {
    if (running == false) {
      printf(ANSI_COLOR_GREEN "\nSIGTERM detected, closing gnuplot pipe safely, terminate program\n" ANSI_COLOR_RESET);
      goto stop_heat_loop;
    }

    for (int j = 0; j < isr-1; j++) {
      // iterations without displaying on gnuplot
      matvec_primme(uk, vk, &blockSize, NULL);
      temperature_iterate(uk, vk, p.n, dt, &t);
    }

    matvec_primme(uk, vk, &blockSize, NULL);
    temperature_iterate(uk, vk, p.n, dt, &t);

    /* generating title */
    sprintf(title, "time : %g s", t); 

    hp.open(&hp, &p, title);
    hp.write(&hp, &p, uk);
  }

  stop_heat_loop:

  hp.close(&hp);

  free(vk); free(uk);

  vspace;
  #endif /*SHOW_TEMPERATURE_EVOL*/


  /* freeing memory */
  free(min_evals); free(max_evals);
  free(min_evecs); free(max_evecs);

  #if SOLVING_WITH_SLEPC
  free(slepc_evals); free(slepc_evecs);
  #endif

  p.close(&p);

  printf("program ended, press ENTER to exit\n");

  return EXIT_SUCCESS;
}