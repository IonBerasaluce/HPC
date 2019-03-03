CC=mpicxx
CXXFLAGS=-O3 -std=c++11

#loading the required libraries 
LDLIBS= -lscalapack-openmpi -lblacs-openmpi -lblacsCinit-openmpi -llapack -lblas -lboost_program_options
OBJS = HPC_cw.o task2funct.o task3funct.o task4funct.o task5funct.o basefunct.o

# # Source Files and function header files
# SRC = task2funct.cpp task3funct.cpp task4funct.cpp task5funct.cpp basefunct.cpp
# HDRS = headers.h

PROG = HPC_cw
SPECS = --length 10 --no_elements 24 --Area 0.012 --MomInertia 0.0000144 --Emodulus 210000000000 --density 7850 --Max_Time 3 --T1 1

#all target will compile all of the cpp files ready for execution
all: $(PROG)

$(PROG): $(OBJS)
	$(CC) -o $(PROG) $(OBJS) $(LDLIBS) 

#Individual target definition for each task including 2f and 3c, after each task the temporary .txt files are deleted 
task1: $(PROG)
	./HPC_cw $(SPECS) --Static 
	python PostProcessTask1.py $(SPECS)
	make -s clean

task2: $(PROG)
	./HPC_cw $(SPECS) --Explicit --Serial --Interval 20000
	python PostProcess.py
	make -s clean

task2f: $(PROG)
	./HPC_cw $(SPECS) --Task2-f --Interval 20000 
	python AmpProcess.py
	make -s clean

task3: $(PROG)
	./HPC_cw $(SPECS) --Implicit --Serial --Interval 20000
	python PostProcess.py
	make -s clean

task3c: $(PROG)
	./HPC_cw $(SPECS) --Task2-f --Interval 10000 
	python AmpProcess.py
	make -s clean
	
task4: $(PROG)
	mpiexec -np 2 ./HPC_cw $(SPECS) --Parallel --Explicit --Interval 100000
	python PostProcess.py
	make -s clean

task5: $(PROG)
	mpiexec -np 4 ./HPC_cw $(SPECS) --Parallel --Implicit --Interval 200000
	python PostProcess.py
	make -s clean

%.o: %.cpp
	mpicxx -c $(CXXFLAGS) $< -o $@

.PHONY: clean

#defining the clean (no object files are cleaned here as we dont want to recompile after each task is run)
clean:
	rm -f *.txt

#After complete execution of the program, the user can delete all of teh object files and temporary .txt files that may remain
fullclean:
	rm -f *.o *.txt

	