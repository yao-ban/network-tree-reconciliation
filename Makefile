all: perecon simulate_species simulate_gene

perecon: perecon.cpp edgeset.h
	g++ perecon.cpp -o perecon

simulate_species: simulate_species.cpp edgeset.h
	g++ simulate_species.cpp -o simulate_species -std=c++11

simulate_gene: simulate_gene.cpp edgeset.h
	g++ simulate_gene.cpp -o simulate_gene -std=c++11
