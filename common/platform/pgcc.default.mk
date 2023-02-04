# (c) 2007 The Board of Trustees of the University of Illinois.

# Rules common to all makefiles

# Commands to build objects from source file using C compiler

# pgcc (default)
CC = pgcc
PLATFORM_CFLAGS = -O3 -fast
  
CXX = pgc++
PLATFORM_CXXFLAGS = -O3 -fast
  
LINKER = pgc++
PLATFORM_LDFLAGS = -lm -lpthread