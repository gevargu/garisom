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
    void get_soilwetness(double& drainage, double& runoff, long& dd, bool& isNewYear, bool& useGSData, long& gs_yearIndex, 
                            double& waterold, long& z, int& layers, double& x, double* thetafracres, double* thetafracfc, double* a, 
                            double* n, double& fieldcapfrac, double* thetafc, double* thetasat, double* water, double* depth,
                            double* fc, double& ffc, double dataCells[2000001][101], long& rowD, long& colD, long& dColF_End_watercontent,
                            double* gs_ar_input, double* gs_ar_waterInitial_OFF, double* gs_ar_waterInitial, std::string& night,
                            double& layerflow, double& baperga, double& lai,double elayer[6][100001], long& halt, double& laisl, long& haltsh,
                            double& laish, double& timestep, double* soilredist, double& deficit, double* swclimit, long& iter_Counter,
                            double& tod, long& stage_ID, bool& mode_predawns, long& dColRain, double& rain, bool& rainEnabled, bool& ground,
                            long* layer, double& waternew, double& waterchange, double* gs_ar_waterFinal, long& dColF_End_waterchange,
                            long& dColF_End_rain, long& dColF_End_gwater, double& gwflow, long& dColF_End_drainage, long& dColF_End_input,
                            long& dColF_End_runoff);
    // Van Genuchten Integration
    void trapzdvg(const long& t,double* a, double* n, double* kmaxrh,const long& z, double& p1, double& p2, double& s, long& it, long& tnm, double& del, double& x, double& sum);
    void qtrapvg(double& olds, long& t,double* a, double* n, double* kmaxrh,const long& z, double& p1, double& p2, double& s, long& it, long& tnm, double& del, double& x, double& sum);
    // Soil E(P) global curve
    void get_rhizoPcrit(long& z, int& layers, double& p1, double& p2,long& k, double& e, double erh[][100001], double krh[][100001], double& pinc, double& kmin, double& olds, long& t,double* a, double* n, double* kmaxrh, double& s, long& it, long& tnm, double& del, double& x, double& sum, double* pcritrh);
};

#endif