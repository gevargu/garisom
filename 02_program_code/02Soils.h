#ifndef SOILS
#define SOILS

#include <iostream>
#include <string>
#include <cmath> // math utility functions
#include <vector>

using namespace std;
class soils
{
public:
    // Van Genuchten functions
    void get_vgparams(std::string &texture, long &layers, std::string (&soillayersTable)[2001][101], long &rowLR, long &colLR);
};

#endif