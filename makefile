# Executable
EXE = program.exe

# Compilers
CC = gcc
FC = gfortran
CXX = g++

# Compiler Flags
CFLG = -Wall -Wextra -pedantic -std=c17
CXXFLG = -Wall -Wextra -pedantic -std=c++14
F77FLG = 
F90FLG = -Wall -Wextra -pedantic

# Optimization 
OPT = -O2 -march=native
CFLG += $(OPT)
CXX += $(OPT)
F77FLG += $(OPT)
F90FLG += $(OPT)

# Compiler Library
CLIB = -lm -lgfortran -lquadmath

# Directories
BDIR = bin
HDIR = include
SDIR = src
ODIR = .obj
DDIR = .dep
MESH = mesh
DATA = data

# Dependencies Flag
DFLG = -MMD -MF $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $@)

# Degbug
DBG = -g

# Veriables
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))
CDEP = $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $(COBJ))

CXXSRC = $(wildcard $(SDIR)/*.cpp)
CXXOBJ = $(patsubst $(SDIR)/%.cpp, $(ODIR)/%.o, $(CXXSRC))
CXXDEP = $(patsubst $(ODIR)/%.o, $(DDIR)/%.d, $(CXXOBJ))

F77SRC = $(wildcard $(SDIR)/*.f)
F77OBJ = $(patsubst $(SDIR)/%.f, $(ODIR)/%.o, $(F77SRC))

F90SRC = $(wildcard $(SDIR)/*.f90)
F90OBJ = $(patsubst $(SDIR)/%.f90, $(ODIR)/%.o, $(F90SRC))

# Targets
.PHONY: all run
all: $(BDIR) $(ODIR) $(DDIR) $(BDIR)/$(EXE)

run: all
	@echo Executing $(EXE):
	@./$(BDIR)/$(EXE)

$(BDIR):
	@mkdir -pv $@

$(ODIR):
	@mkdir -pv $@

$(DDIR):
	@mkdir -pv $@

$(BDIR)/$(EXE): $(COBJ) $(CXXOBJ) $(F77OBJ) $(F90OBJ)
	$(CC) -o $@ $^ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(DFLG) $(DBG) $(CFLG) -o $@ -c $< -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CXX) $(DFLG) $(DBG) $(CXXFLG) -o $@ -c $< -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.f
	$(FC) $(F77FLG) -o $@ -c $<

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) $(F90FLG) $(DBG) -o $@ -c $<

# Clean
.PHONY: clean clean_all clean_data clean_data_all
clean: 
	@$(RM) -rv $(BDIR)/* $(ODIR)/* $(DDIR)/*

clean_all: clean
	@$(RM) -rv $(BDIR) $(ODIR) $(DDIR)

clean_data:
	@find ./$(DATA)/Test/ -type f -name '*.pdf' | xargs $(RM) -rv
	@find ./$(DATA)/Test/ -type f -name '*.txt' | xargs $(RM) -rv
	@find ./$(DATA)/Test/ -type f -name '*.log' | xargs $(RM) -rv
	@find ./$(DATA)/Test/ -type f -name '*.aux' | xargs $(RM) -rv
	@find ./$(DATA)/Test/ -type f -name '*.tex' | xargs $(RM) -rv

clean_data_all: clean_all clean_data
	@find ./$(DATA)/Test/ -type f -name '*.dat' | xargs $(RM) -rv

# Debug
.PHONY: valgrind
valgrind:
	valgrind --leak-check=full ./$(BDIR)/$(EXE) 

# Include Dependencies
-include $(CDEP) $(CXXDEP)