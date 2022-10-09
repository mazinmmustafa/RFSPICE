# Compilers
CC = gcc

# Executable
EXE = program.exe

# Compiler Libraries
CLIB = -lm

# Flags
CFLG = -Wall -Wextra -Ofast

# Directories
BDIR = bin
HDIR = include
ODIR = obj
SDIR = src

# Variables
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))

# Targets
all: $(BDIR)/$(EXE)

$(BDIR)/$(EXE): $(COBJ) $(FOBJ) $(F77OBJ)
	$(CC) -o $@ $^ -I $(HDIR) -L./ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLG) -c -o $@ $^
	
# Clean
.PHONY: clean

clean:
	del $(BDIR)\$(EXE) $(ODIR)\*.o
	cls