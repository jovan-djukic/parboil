# (c) 2007 The Board of Trustees of the University of Illinois.

# Rules common to all makefiles

# Commands to build objects from source file using C compiler
# with gcc

# gcc (default)
CC = $(OPENACC_PATH)/bin/pgcc
PLATFORM_CFLAGS =-O3 -fast -acc -cuda -Minfo=accel -I$(OPENACC_PATH)/include -L$(OPENACC_PATH)/lib
  
CXX = $(OPENACC_PATH)/bin/pgc++
PLATFORM_CXXFLAGS =-O3 -fast -acc -cuda -Minfo=accel -I$(OPENACC_PATH)/include -L$(OPENACC_PATH)/lib
  
LINKER = $(OPENACC_PATH)/bin/pgc++
PLATFORM_LDFLAGS = -lm -lpthread -O3 -fast -acc -cuda -Minfo=accel

