# PRIMME
LIBP = -L./primme/ -lprimme
INCP = -I./primme/PRIMMESRC/COMMONSRC/

# ALL
LIB = $(LIBP) -lm -lblas -llapack

objects = prob.o gnuplot.o temperature.o time.o interface_primme.o interface_slepc.o
headers = $(objects:.c=.h)

COPT = -O2

default: executable_to_wrap 

# defines PETSC_DIR, PETSC_ARCH, SLEPC_DIR
include libsources 

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

# this extends the clean defined in slepc_common
clean:: 
	rm executable_to_wrap 

executable_to_wrap: main.c $(objects) $(headers) config.h
	$(LINK.C) $(COPT) $^ -o $@ ${SLEPC_EPS_LIB} $(LIB) 

%.o: %.c config.h
	$(LINK.C) $(COPT) -c $< -o $@ ${SLEPC_EPS_LIB} $(INCP) 