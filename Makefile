all: garisom

garisom: gainriskSOModelCPP.cpp
	g++ -std=c++11 -O3 -ffast-math gainriskSOModelCPP.cpp -o garisom
