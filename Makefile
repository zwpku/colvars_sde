COLVARS_LIB = ./src/colvars/libcolvars.a
COLVARS_SRC_DIR = ./src/colvars
MAIN_LIB = ./src/backend/main.o
POTENTIAL_LIB = ./src/backend/potentials.o

MAIN_NAME = sde_colvars

PREFIX = ~/local/bin/

CXXFLAGS := -std=c++11 -pedantic -g -O2 -fPIC

#TORCHDIR = /home/numerik/bzfzhang/local/libtorch
TORCHDIR = /home/wei/local/libtorch
TORCHINCFLAGS = -I$(TORCHDIR)/include -I$(TORCHDIR)/include/torch/csrc/api/include
EXTRACOLVARSFLAGS = -std=c++17 -DTORCH $(TORCHINCFLAGS)
EXTRALINKLIBS = -Wl,-rpath,$(TORCHDIR)/lib -L$(TORCHDIR)/lib -ltorch -ltorch_cpu -lc10
COLVARS_DEBUG=

#.PHONY: $(COLVARS_LIB)

default: $(MAIN_NAME)

$(MAIN_NAME): $(MAIN_LIB) $(POTENTIAL_LIB) $(COLVARS_LIB) 
	g++ -o $@ $< $(POTENTIAL_LIB) $(COLVARS_LIB) $(EXTRALINKLIBS) 

$(MAIN_LIB): ./src/backend/main.cpp
	make -C ./src/backend main.o

$(POTENTIAL_LIB): ./src/backend/potentials.cpp
	make -C ./src/backend potentials.o

$(COLVARS_LIB):
	COLVARS_DEBUG=$(COLVARS_DEBUG) EXTRACOLVARSFLAGS="$(EXTRACOLVARSFLAGS)" make -C $(COLVARS_SRC_DIR) libcolvars.a -j4

install: $(MAIN_NAME)
	cp $< $(PREFIX)

clean:
	rm -f $(MAIN_NAME) 
	make -C ./src/backend clean
	make -C $(COLVARS_SRC_DIR) clean

