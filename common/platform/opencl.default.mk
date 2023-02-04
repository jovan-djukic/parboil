# (c) 2007 The Board of Trustees of the University of Illinois.

# Rules common to all makefiles

# Commands to build objects from source file using C compiler
# with gcc

# gcc (default)
# Variable provided by parboil
# CC = gcc

# Variable provided by parboil
# CXX = g++
  
# LINKER = g++

# Varaibles necessary to compile OpenCL code, OPENCL_PATH variable is defined in Makefile.conf
CC = $(OPENCL_PATH)/bin/nvcc
PLATFORM_CFLAGS = -O3

CXX = $(OPENCL_PATH)/bin/nvcc
PLATFORM_CXXFLAGS = -O3

LINKER = $(OPENCL_PATH)/bin/nvcc
PLATFORM_LDFLAGS = -lm -lpthread

