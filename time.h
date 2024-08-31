double mytimer_cpu();
double mytimer_wall();

/*
    Set of macros to make time measurements troughout the code cleaner

    /!\ if they use the same ti,tf : timinig macros should not be nested
*/

#define tictac(function, action, timing_function, ti, tf)\
    ti = timing_function();\
    function;\
    tf = timing_function();\
    printf("-> time taken for ");\
    printf("%s was %e seconds\n", action, tf-ti);

#define tic(timing_function, ti)\
    ti = timing_function();

#define tac(timing_function, tf, action)\
    tf = timing_function();\
    printf("-> time taken for ");\
    printf("%s was %e seconds\n", action, tf-ti);
