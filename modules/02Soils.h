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
    // soil structure
    double get_depths(const long& k,const double& layers,const double& beta);
    double get_halfdepths(const long& k,const double& layers,const double& beta);
    double get_layvolume(const double& depth,const double& pi,const double& depthmax,const double& aspect);
    double get_radius(const double& vol,const double& depth,const double& pi);
    // rhizosphere resistance
    double get_rhizor(const double& rplant,const double& rhizotarg);
    double get_kmaxrh(double& rhizor, double& vgterm);
    // Van Genuchten functions
    std::vector<double> get_vgparams(const std::string& texture);
    double get_vgp(double* a, double* n,const double& x,const long& z);
    double get_vgterm(double* n,const double& vp,const long& z);
    double get_vg(double* a, double* n, const double& x,double* kmaxrh,const long& z);
    // Soil water status
    double get_rvg(double* a, double* n, const double& x, const long& z);
    double get_swc(double* a, double* n, const double& x, const long& z);
    // Van Genuchten Integration
    void trapzdvg(const long& t,double* a, double* n, double* kmaxrh,const long& z, double& p1, double& p2, double& s, long& it, long& tnm, double& del, double& x, double& sum);
    void qtrapvg(double& olds, long& t,double* a, double* n, double* kmaxrh,const long& z, double& p1, double& p2, double& s, long& it, long& tnm, double& del, double& x, double& sum);
    // Soil E(P) global curve
    void rhizocurves(long& z, int& layers, double& p1, double& p2,long& k, double& e, double erh[][100001], double krh[][100001], double& pinc, double& kmin, double& olds, long& t,double* a, double* n, double* kmaxrh, double& s, long& it, long& tnm, double& del, double& x, double& sum, double* pcritrh);
};

#endif