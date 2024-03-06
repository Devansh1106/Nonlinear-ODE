# This is a simple makefile just to build PROG
#COMP_FLAGS = -O0 -ggdb3 -fsanitize=address

#COMP_FLAGS = -O3

DEF_FLAGS_CPP  = $(COMP_FLAGS) -mcmodel=small	-Wall	-g -Wextra

DEF_FLAGS_CC = $(COMP_FLAGS) -g

CPP := mpic++ $(DEF_FLAGS_CPP)

CC := mpic++ $(DEF_FLAGS_CC)

DIR_SRC_CPP := $(wildcard src/*.cpp) 

DIR_SRC_CC:= src/Solver.c

CFLAGS = -I C/MUMPS_5.5.1/include -I include

PROG_CPP = main_cpp

####################################################################
# Lines from here to down should not be changed. They are the actual
# rules which make uses to build PROG
.KEEP_STATE:

all:$(PROG_CPP)

OBJS_CPP = $(DIR_SRC_CPP:%.cpp=%.o)
OBJS_MPI = $(DIR_SRC_CC:%.c=%.o)
$(PROG_CPP):$(OBJS_CPP) $(OBJS_MPI)
	mpic++ -o $(PROG_CPP) $(OBJS_CPP) $(OBJS_MPI) $(CFLAGS) -ldmumps -lmumps_common -lpord -lmetis -lscotch -lscotcherr -lmpi -lstdc++
	rm -f *.o $(OBJS_CPP) $(OBJS_MPI)

run: all
	mpirun -np 8 ./$(PROG_CPP)
	# mpirun -np 8 sh -c 'if [ $$OMPI_COMM_WORLD_RANK -eq 0 ]; then valgrind --leak-check=full ./$(PROG_CPP) > valgrind_output.txt 2>&1; else ./$(PROG_CPP); fi'

	#./$(PROG_MPI)