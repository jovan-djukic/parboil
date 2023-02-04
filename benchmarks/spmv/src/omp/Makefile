# (c) 2010 The Board of Trustees of the University of Illinois

LANGUAGE=pgcc
TOOLS_SRC=common_src/convert-dataset
SRCDIR_OBJS=main.o file.o 
# APP_CFLAGS=-I$(TOOLS_SRC) -fopenmp
APP_CFLAGS=-I$(TOOLS_SRC) -mp -Minfo=mp
APP_LDFLAGS=-lgomp

include $(TOOLS_SRC)/commontools.mk
