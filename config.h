#ifndef CONFIG_MEMBRANE_H
#define CONFIG_MEMBRANE_H

#define M_UNIT_STEPS 20
// number of points for the lenght of a unit square

#define EXTRACT_MAT 1
/* will extract the matrix to ./compare_mat
/!\ to verify with reference matrix, use M_UNIT_STEPS=3 
and then do ./diff_all.sh inside folder */

#define PRIMME_CONFIG_PRINT 1
// print solver config before solving

#define SLEPC_CONFIG_PRINT 1
#define SOLVING_WITH_SLEPC 1

#define CALC_RESIDUAL 1
// calculates residual from min eigen value of primme

#define DISPLAY_EIGENVEC 1
// display the eigen vector for min eigen value from primme

#define SHOW_TEMPERATURE_EVOL 1
// DT max is the limit for the progressive euler method to converge
#define DT_F 10 // fraction of DT max, dt = dt_max / DT_F
#define TOTAL_TIME 10000 // in seconds, not exact, targeted
#define ISR 100 // inverse sampling rate, for instance, 10 would mean that 1/10 iteration will be displayed
#define INITIAL_TEMP 10 // initial temperature
#define DIFFUSIVITY 9.7e-5 // diffusivity

#define broadcast(msg)\
  printf(ANSI_COLOR_YELLOW "-----");\
  printf(msg);\
  printf("-----" ANSI_COLOR_RESET);\
  printf("\n");

#define vspace printf("\n")

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#endif // !CONFIG_MEMBRANE_H