COLVARS_LIB = ./colvars/src/libcolvars.a
COLVARS_SRC_DIR = ./colvars/src

CXXFLAGS := -std=c++11 -pedantic -g -O2 -fPIC

TORCHDIR = /home/numerik/bzfzhang/local/libtorch
TORCHINCFLAGS = -I$(TORCHDIR)/include -I$(TORCHDIR)/include/torch/csrc/api/include
EXTRACOLVARSFLAGS = -std=c++17 -DTORCH $(TORCHINCFLAGS)
EXTRALINKLIBS = -Wl,-rpath,$(TORCHDIR)/lib -L$(TORCHDIR)/lib -ltorch -ltorch_cpu -lc10

default: main

main: 	main.o $(COLVARS_LIB) 
	g++ -o $@ $<  $(COLVARS_LIB) $(EXTRALINKLIBS)

main.o: src/main.cpp
	g++ -c $< -pedantic -I$(COLVARS_SRC_DIR)

$(COLVARS_LIB):
	EXTRACOLVARSFLAGS="$(EXTRACOLVARSFLAGS)" make -C $(COLVARS_SRC_DIR) libcolvars.a -j4

