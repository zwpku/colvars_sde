COLVARS_LIB = ./colvars/src/libcolvars.a
COLVARS_SRC_DIR = ./colvars/src
MAIN_LIB = ./src/main.o

CXXFLAGS := -std=c++11 -pedantic -g -O2 -fPIC

TORCHDIR = /home/numerik/bzfzhang/local/libtorch
#TORCHDIR = /home/wei/local/libtorch
TORCHINCFLAGS = -I$(TORCHDIR)/include -I$(TORCHDIR)/include/torch/csrc/api/include
EXTRACOLVARSFLAGS = -std=c++17 -DTORCH $(TORCHINCFLAGS)
EXTRALINKLIBS = -Wl,-rpath,$(TORCHDIR)/lib -L$(TORCHDIR)/lib -ltorch -ltorch_cpu -lc10

.PHONY: $(COLVARS_LIB)

default: main

main: 	$(MAIN_LIB) $(COLVARS_LIB) 
	g++ -o $@ $<  $(COLVARS_LIB) $(EXTRALINKLIBS)

$(MAIN_LIB): ./src/main.cpp
	make -C ./src main.o

$(COLVARS_LIB):
	EXTRACOLVARSFLAGS="$(EXTRACOLVARSFLAGS)" make -C $(COLVARS_SRC_DIR) libcolvars.a -j4

