# Makefile for CUDA Event Generator
# Written by Siddharth Sule, 2023

# Directories
BINDIR = bin
INCDIR = include
OBJDIR = obj
SRCDIR = src

# Files
SRCS = $(wildcard $(SRCDIR)/*.cpp $(SRCDIR)/*.cu)
INCS = $(wildcard $(INCDIR)/*.h $(INCDIR)/*.cuh)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(patsubst $(SRCDIR)/%.cu, $(OBJDIR)/%.o, $(SRCS)))
PROG = $(BINDIR)/main

# Compiler
NVCC = nvcc

# Flags
NVCCFLAGS = -std=c++17 -I$(INCDIR)

# Rules
all: directories $(PROG)

directories:
        mkdir -p $(BINDIR)
        mkdir -p $(OBJDIR)

$(PROG): $(OBJS)
        $(NVCC) $(NVCCFLAGS) -o $(PROG) $(OBJS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCS)
        $(NVCC) $(NVCCFLAGS) -c -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cu $(INCS)
        $(NVCC) $(NVCCFLAGS) -c -o $@ $<

clean:
        rm -f $(PROG) $(OBJS)
        rm -rf $(BINDIR)
        rm -rf $(OBJDIR)

.PHONY: all clean