# Note: This Makefile requires GMP, MPFR, OpenMPI and fpLLL to be pre-installed
#       and all libraries, include files and binaries to be available via the
#       path and related environment variables in the shell.

cc  = mpicc
CC  = mpiCC
LD  = $(CC)

MKDIR = mkdir
RM = rm
DOXYGEN = doxygen

# For debug tracing: # -D DEBUG_TRACE_RNG -D DEBUG_TRACE_SAMPLING
cc_FLAGS = -m64 -O4 -Wall -Wextra
CC_FLAGS = $(cc_FLAGS) -std=c++11 -D OMPI_SKIP_MPICXX
LD_FLAGS  = -lgmp -lmpfr -lfplll

include common.mk
