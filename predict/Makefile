OBJS = random_forest.o io.o parse_csv.o cluster.o options.o \
	mpi_util.o one_hot.o

CC = mpic++
cc = mpicc

PROFILE = #-pg -g
OPENMP = -fopenmp

FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++0x

INC = -I. 
LIBS = -lm -lz

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(PROFILE) $(INC) -c $<
.c.o:
	$(cc) $(FLAGS) $(PROFILE) $(INC) -c $<

all: predict_synergy

predict_synergy : $(OBJS) main.o
	$(CC) $(PROFILE) -o predict_synergy main.o $(OBJS) $(LIBS) $(OPENMP)

clean:
	-rm -f *.o


