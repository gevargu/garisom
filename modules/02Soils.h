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
    double get_depths(long k,double layers, double beta);
    double get_halfdepths(long k, double layers, double beta);
    double get_layvolume(double depth, double pi, double depthmax, double aspect);
    double get_radius(double vol, double depth, double pi);
    // rhizosphere resistance
    double get_rhizor(double rplant, double rhizotarg);
    // Van Genuchten functions
    std::vector <double> get_vgparams(std::string texture);
    double get_vgp(double a[6], double n[6], double &x, long z);
    double get_vgterm(double n[6], double vp, long z);
    double get_vg(double a[], double n[], double &x,double kmaxrh[], long z);
    // Soil water status
    double get_rvg(double a[6], double n[6], double x, long z);
    double get_swc(double a[6], double n[6], double x, long z);
    void rhizor_change(double* rhizor, double* kmaxrh[6]);
    // Van Genuchten Integration
    //void trapzdvg(double* p1, double* p2, double* s, long* t);
    //void qtrapvg(double &p1, double &p2, double &s);
    // Soil E(P) global curve
    //void rhizocurves(long &z, int &layers, double &p1, double &p2, long &k, double &e,double kmaxrh[6], double pinc, double &s, double &x, double a[6], double n[6], double kmin);
};

#endif