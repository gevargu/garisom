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
    // genaral
    double get_cweibull(const double& p12, const double& p50);
    double get_bweibull(const double& p12, const double& c);
    double get_weibullfit(const double& x,const double& b, const double& c, const double& ksat);
    // root block
    double get_rootpercent(const double& leafpercent);
    double get_rootkmax(const double& rootpercent,const double& rsatp);
    double get_weibullfitroot(const double& x, const double& b, const double& c, double* ksatr, const int& z);
    void trapzdwbr(const long& t, double& p1, double& p2, const double& b, const double& c,double* ksatr, const int& z, double& s, long& it, long& tnm, double& del, double& x, double& sum);
    void qtrapwbr(double& olds, long& t, double& p1, double& p2, const double& b, const double& c,double* ksatr, const int& z, double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx);
    void get_rootPcrit(double er[][100001],double kr[][100001],long& z, const int& layers,double* ksatr, double& p1, double& p2, const double& pinc,long& k,double& e,double& olds, long& t, const double& b, const double& c, double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx, const double& kmin,double* pcritr);
    void get_rootPcrit_v(double er_v[][100001],double kr_v[][100001],long& z, const int& layers,double* ksatr, double& p1, double& p2, const double& pinc,long& k,double& e,double& olds, long& t, const double& b, const double& c, double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx, const double& kmin,double* pcritr);
    // stem block
    double get_stempercent(double leafpercent);
    double get_stemkmax(double stempercent, double rsatp);
    void trapzdwbs(const long &t, double &p1, double &p2, const double& b, const double& c, const double& ksats,double &s, long& it, long& tnm, double& del, double& x, double& sum);
    void qtrapwbs(double& olds, long &t, const long& f, double &p1, double &p2, const double& b, const double& c, const double& ksats,double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx);
    
    // leaves block
    double get_gmaxl(double gmax, double laperba);
    double get_leafkmax(double ksatp, double leafpercent);
    double get_LSC(double ksatl, double laperba);
};

#endif