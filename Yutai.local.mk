#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/usr/local/opt/openblas/lib -lopenblas -L/usr/local/opt/lapack/lib -llapack -L/usr/local/lib -lm 
INCLUDEBLASLOCAL=-I/usr/local/include -I/usr/local/opt/openblas/include -I/usr/local/opt/lapack/include
OPTCLOCAL=-fPIC -march=native
