OBJS = parse_combo.o parse_csv.o io.o dose_response_function.o

CC = g++
cc = gcc

PROFILE = #-pg -g
OPENMP = -fopenmp

FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++0x

# This Makefile assumes that the GNU Scientific Library is installed
# in your home directory in a subdirectory called "gsl". If this is
# not the case, you will need to edit the "GSL_LIB_DIR" and
# "GSL_INCLUDE_DIR" variables below to refer to the directories that
# contain the GSL library files and the GSL include files respectively.
GSL_LIB_DIR = $(HOME)/gsl/lib
GSL_INCLUDE_DIR = $(HOME)/gsl/include

INC = -I. -I$(GSL_INCLUDE_DIR)
LIBS = -lm -lz $(GSL_LIB_DIR)/libgsl.a $(GSL_LIB_DIR)/libgslcblas.a 

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(PROFILE) $(INC) -c $<
.c.o:
	$(cc) $(FLAGS) $(PROFILE) $(INC) -c $<

all: gemini_prep synergy_search

gemini_prep : $(OBJS) gemini_prep.o
	$(CC) $(PROFILE) -o gemini_prep gemini_prep.o $(OBJS) $(LIBS) $(OPENMP)


synergy_search : $(OBJS) synergy_search.o
	$(CC) $(PROFILE) -o synergy_search synergy_search.o $(OBJS) $(LIBS) $(OPENMP)
	
clean:
	-rm -f *.o


