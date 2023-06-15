/* 
# Function to calculate the initial conditions in the model

## This module calls other modules in subdirectories

*/
#ifndef MISCFUNCTIONS
#define MISCFUNCTIONS

// Necessary cpp libraries
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// Some miscellaneous functions
class MiscFunctions
{
    public:
    double get_patm(double alt);
    double get_resist(double conductance);
};

#endif