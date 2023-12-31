COLVARS_SRC_DIR = ../colvars

main.o: main.cpp
	g++ -c $< -pedantic -I$(COLVARS_SRC_DIR) -I./

potentials.o: potentials.cpp
	g++ -c $< -pedantic 

clean: 
	rm -f *.o


