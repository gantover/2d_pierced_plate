# Installation instructions

the **install** file will 
- automatically download all the required librairies (including primme) and configure them. 
- configure the paths in libsources for the makefile.
- download the matrix (ia,ja,a) from http://metronu.ulb.ac.be 

> /!\ wget should be installed on your system beforehand : `apt install wget`

1. give execution permission to the file with `chmox +x install`
2. execute `./install`

> There should not be anything else to do before compiling the program

# Compilation

- `make`

# Execution instructions

> the program has to be run inside a wrapper, 
this wrapper will seek for the **ENTER** key to be pressed and send
a **SIGTERM** signal to the c program for it to handle the *closing of the
gnuplot pipe safely*

to execute program : `./wrapper`

> you may need to do `chmod +x wrapper`

# Change program settings

config.h has a lot of parameters to play with. Those options act as toogles for bits of code, 1 being on and 0 off.

# Note for debugging with valgrind

> Please disable the slepc solver part with `SOLVE_WITH_SLEPC 0` 

Even the exemples provided on the slepc website make memory errors upon running inside valgrind : https://petsc-users.mcs.anl.narkive.com/hufTDZv6/tons-of-valgrind-errors-with-simply-petscinitialize-and-petscfinalize

# Verify project matrix 

To verify the matrix of this project with the one obtained from the metronu website, 
use `define M 3`, `define EXTRACT_MAT 1`, execute the program and then go into compare_mat and do `./diff_all.sh`