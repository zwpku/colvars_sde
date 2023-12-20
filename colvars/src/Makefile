# -*- makefile -*-

# This Makefile is useful for compilation tests only.  To use the library,
# please use instead the build recipes of the packages that include it.

COLVARS_LIB = libcolvars.a
COLVARS_OBJ_DIR =

.PHONY: default clean

default: $(COLVARS_LIB)

.SUFFIXES:
.SUFFIXES: .cpp .o

ifeq ($(CXX),)
CXX := g++
endif
ifeq ($(CXXFLAGS),)
CXXFLAGS := -O2 -fPIC
ifneq ($(CXX),"CC")
CXXFLAGS += -Wall
endif
endif

ifeq ($(AR),)
AR := ar
ARFLAGS := -rscv
endif
ifeq ($(SHELL),)
SHELL := /bin/sh
endif


# Detect debug settings
ifeq ($(COLVARS_DEBUG),)
COLVARS_DEBUG_INCFLAGS :=
else
COLVARS_DEBUG_INCFLAGS := -DCOLVARS_DEBUG
endif

COLVARS_INCFLAGS := $(COLVARS_DEBUG_INCFLAGS)

COLVARS_SRCS := $(wildcard *.cpp)

COLVARS_OBJS := $(COLVARS_SRCS:.cpp=.o) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(EXTRACOLVARSFLAGS) $(COLVARS_INCFLAGS) -c -o $@ $< 


$(COLVARS_LIB):	Makefile.deps $(COLVARS_OBJS)
	$(AR) $(ARFLAGS) $(COLVARS_LIB) $(COLVARS_OBJS)  

Makefile.deps: $(COLVARS_SRCS)
	@echo > $@
	@for src in $^ ; do \
	  obj=`basename $$src .cpp`.o ; \
	  $(CXX) -MM $(COLVARS_INCFLAGS)  \
	    -MT '$$(COLVARS_OBJ_DIR)'$$obj $$src >> $@ ; \
	  done

include Makefile.deps

clean:
	-rm -f $(COLVARS_OBJS) $(COLVARS_LIB) Makefile.deps

