# Make file for GARISOM
#
# Author: German Vargas
# Date: June 2023.
#
# This files loads and set the instructions for compiling GARISOM.
# Each line is an compilation instruction or the location of a module:
# gainriskSOM_2.0.cpp
# TEST.cpp

all: garisom
 
garisom: gainriskSOM_2.0.cpp
	g++ -std=c++11 -O3 -ffast-math gainriskSOM_2.0.cpp \
	modules/00MainProgram.cpp \
	modules/01IOHandler.cpp \
	modules/02Soils.cpp \
	modules/03Hydraulics.cpp \
	modules/04Morphology.cpp \
	modules/05CAssimilation.cpp	\
	modules/06MiscFunctions.cpp \
	-o garisom