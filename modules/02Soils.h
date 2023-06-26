#ifndef SOILS
#define SOILS

#include <iostream>
#include <string>
#include <cmath> // math utility functions
#include <vector>
#include "01Macros.h"

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
    // Soil water status & water potential
    double get_rvg(double* a, double* n, const double& x, const long& z);
    double get_swc(double* a, double* n, const double& x, const long& z);
    void get_soilwetness(double& drainage, double& runoff, double& waterold, double& x, double* thetafracres, double* thetafracfc, double* a,
                            double* n, double& fieldcapfrac, double* thetafc, double* thetasat, double* water, double* depth, double* fc,
                            double& ffc, double dataCells[2000001][101], double* gs_ar_input, double* gs_ar_waterInitial_OFF, double* gs_ar_waterInitial, 
                            double& layerflow, double& baperga, double& lai,double elayer[6][100001], double& laisl, double& laish, double& timestep,
                            double* soilredist, double& deficit, double* swclimit, double& tod, double& rain, double& waternew, double& waterchange,
                            double* gs_ar_waterFinal, double& gwflow, double& transpirationtree, double& laperba, double& soilevap, double* gs_ar_E,
                            double* gs_ar_ET, double& atree, double* gs_ar_Anet, double* gs_ar_cica, double& cinc, double& ca, long* gs_ar_cica_N,
                            double* gs_ar_Aci, double* gs_ar_AnetDay, double& sum, double& kpday1, double& kxday1, double* gs_ar_waterInitial_GS,
                            double* gs_ar_waterFinal_OFF, double& iter_refK, double* gs_ar_PLCSum, double* gs_ar_PLCSum_N, double* gs_ar_PLCp,
                            double* gs_ar_PLC85, double* gs_ar_PLCx, double* gs_ar_waterInput_GS, double* gs_ar_waterFinal_GS, double* gs_ar_waterInput_OFF,
                            int& layers, std::string& night,
                            long& dd, long& gs_yearIndex, long& z, long& rowD, long& colD, long& dColF_End_watercontent, long& halt,  long& haltsh,
                            long& iter_Counter, long& stage_ID, long& dColRain, long* layer, long& dColF_End_waterchange, long& dColF_End_rain,
                            long& dColF_End_gwater, long& dColF_End_drainage, long& dColF_End_input, long& dColF_End_runoff, long& dColF_End_E,
                            long& o, long& dColF_End_soilEvap, long& dColF_End_ET, long& dColF_End_ANet, long& dColF_CP_kplant, long& dColF_CP_kxylem,
                            long& dColF_End_PLCplant, long& dColF_End_PLCxylem,
                            bool& isNewYear, bool& useGSData, bool& mode_predawns, bool& rainEnabled, bool& ground, bool& gs_inGrowSeason,
                            bool& gs_doneFirstDay);
    void get_predawns(long& k, long* layer, long& z, long& rowD, long& colD, long& dd, long& dColRain, long& t, long& failure, long& dColF_p1, long& o,
                            int& layers, std::string* layerfailure, std::string& failspot,bool& mode_predawns,
                            double* kminroot, double& theta, double* water, double& x, double* depth, double& pgrav,
                            double* thetasat, double dataCells[2000001][101], double* pd, double* a, double* n,
                            double* prh, double* pcritrh, double* pcritr, double& sum, double& pr, double& prinitial);
    // Van Genuchten Integration
    void trapzdvg(const long& t,double* a, double* n, double* kmaxrh,const long& z, double& p1, double& p2, double& s, long& it, long& tnm, double& del, double& x, double& sum);
    void qtrapvg(double& olds, long& t,double* a, double* n, double* kmaxrh,const long& z, double& p1, double& p2, double& s, long& it, long& tnm, double& del, double& x, double& sum);
    // Soil E(P) global curve
    void get_rhizoPcrit(long& z, int& layers, double& p1, double& p2,long& k, double& e, double erh[][100001], double krh[][100001], double& pinc, double& kmin, double& olds, long& t,double* a, double* n, double* kmaxrh, double& s, long& it, long& tnm, double& del, double& x, double& sum, double* pcritrh);
};

#endif