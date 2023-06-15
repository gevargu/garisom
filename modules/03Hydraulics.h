#ifndef HYDRAULICS
#define HYDRAULICS

#include <stdio.h>
#include <iostream>
#include <string> // the C++ String Class, easier to deal with than char arrays for this application
#include <fstream>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <cmath> // math utility functions
#include <ctime> // timers for performance testing
#include <iterator>
#include <cstring>

using namespace std;
class hydraulics
{
public:
    double get_cweibull(double p12, double p50);
    double get_bweibull(double p12, double c);
    double get_weibullfit(double &x, double b, double c, double ksat);
    double get_rootpercent(double leafpercent);
    double get_rootkmax(double rootpercent, double rsatp);
    double get_weibullfitroot(double &x, double b, double c, double ksatr[6], int z);
    void trapzdwbr(double &p1, double &p2, double &s, long &t);
    void qtrapwbr(double &p1, double &p2, double &s);
    double get_kmaxrh(double rhizor, double vgterm);
    void trapzdwbs(double &p1, double &p2, double &s, long &t);
    void qtrapwbs(double &p1, double &p2, double &s);
    double get_stempercent(double leafpercent);
    double get_stemkmax(double stempercent, double rsatp);
    double get_gmaxl(double gmax, double laperba);
    double get_leafkmax(double ksatp, double leafpercent);
    double get_LSC(double ksatl, double laperba);
};

#endif