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
    // Van Genuchten functions
    void get_vgparams(std::string &texture, long &layers, std::string (&soillayersTable)[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS], long &rowLR, long &colLR);
    // Soil wetness
    void get_soilwetness(double &drainage, double &runoff, double &waterold, double &x, double (&thetafracres)[6], double (&a)[6],
        double (&n)[6], double (&thetafracfc)[6], double (&thetafc)[6], double (&depth)[6], double &fieldcapfrac, double (&thetasat)[6], double (&water)[6],
        double (&fc)[6], double &ffc, double (&dataCells)[DATAFILE_MAXROWS][DATAFILE_MAXCOLS], double (&gs_ar_input)[100], double (&gs_ar_waterInitial_OFF)[100],
        double (&gs_ar_waterInitial)[100], double &layerflow, double (&elayer)[6][100001], double &laisl, double &lai, double &laish,
        double &baperga, double &timestep, double (&soilredist)[6], double &deficit, double (&swclimit)[6], double &tod, double &rain,
        double &waternew, double &waterchange, double (&gs_ar_waterFinal)[100], double &gwflow, double &transpirationtree, double &laperba, 
        double (&gs_ar_E)[100], double &soilevap, double (&gs_ar_ET)[100], double &atree, double (&gs_ar_Anet)[100], double &cinc, double &ca,
        double (&gs_ar_cica)[100], long (&gs_ar_cica_N)[100], double (&gs_ar_Aci)[100], double (&gs_ar_AnetDay)[100], double &sum, double &kpday1,
        double &kxday1, double (&gs_ar_waterInitial_GS)[100], double (&gs_ar_waterFinal_OFF)[100], double &iter_refK, double (&gs_ar_PLCSum)[100],
        double (&gs_ar_PLCSum_N)[100], double (&gs_ar_PLCp)[100], double (&gs_ar_PLC85)[100], double (&gs_ar_PLCx)[100], double (&gs_ar_waterInput_GS)[100],
        double (&gs_ar_waterFinal_GS)[100], double (&gs_ar_waterInput_OFF)[100], long &dd, long &gs_yearIndex, long &z, long &layers, long &rowD,
        long &colD, long &dColF_End_watercontent, long &halt, long &haltsh, long &iter_Counter, long &stage_ID, long &dColRain, long &j, long (&layer)[6],
        long &dColF_End_waterchange, long &dColF_End_rain, long &dColF_End_gwater, long &dColF_End_drainage, long &dColF_End_input, long &dColF_End_runoff,
        long &o, long &dColF_End_E, long &dColF_End_soilEvap, long &dColF_End_ET, long &dColF_End_ANet, long &dColF_CP_kplant, long &dColF_CP_kxylem,
        long &dColF_End_PLCplant, long &dColF_End_PLCxylem, std::string &night, bool &isNewYear, bool &useGSData, bool &mode_predawns, bool &rainEnabled,
        bool &ground, bool &gs_inGrowSeason, bool &gs_doneFirstDay);
    
    // Soil predawn water potential
    void get_predawns(double (&kminroot)[6], double &theta, double (&water)[6], double (&depth)[6], double (&thetasat)[6],
        double &x, double (&pd)[6], double (&dataCells)[DATAFILE_MAXROWS][DATAFILE_MAXCOLS], double &pgrav, double (&a)[6], 
        double (&n)[6], double (&prh)[6], double (&pcritrh)[6], double (&pcritr)[6], double &sum, double &pr, double &prinitial,
        long &k, long &layers, long (&layer)[6], long &z, long &rowD, long &dd, long &colD, long &dColRain, long &t, long &failure,
        long &dColF_p1, long &o, std::string (&layerfailure)[6], std::string &failspot, bool &mode_predawns);
    // Soil water flow 
    void get_soilflow(double &store, double (&kmaxrh)[6], double (&kkmax)[6], double (&depth)[6], double (&pd)[6],
        double (&soilredist)[6], double &p1, double &pend, double &e, double &p2, double &pinc, double &olds, double (&a)[6],
        double (&n)[6], double &x, double &s, double &del, double &sum, double &eps, double (&soilf)[6], long &z, long &layers,
        long &tmax, long &t, long &it, long &tnm, long &j);
    // Deepsoil water flow
    void get_deepflow(double &store, double (&kmaxrh)[6], double (&kkmax)[6], double &grounddistance, double &pground,
        double (&pd)[6], double &p1, double &pend, double &e, double &p2, double &pinc, double &olds, double (&a)[6],
        double (&n)[6], double &x, double &s,double &del, double &sum, double &eps, double &groundflow, double &gwflow, double &timestep,
        double &drainage, double (&soilredist)[6],long &z, long &layers, long &tmax, long &t, long &it, long &tnm, long &j);
    // Soil Evaporation
    void get_soilevaporation(double &emission, double &emiss, const double &sbc, double& soiltemp, double& lc,
        double (&leaftempsh)[100001], double &rabssoil, double &soilabssol, double &sbottom, double &mdensair,
        double &patm, double &airtemp, double &gha, double &us, double &xh, double &zdispl, double &rough, double &zh,
        double &soilep, const double &sha, double &grad, double &lambda, double &rha, double &vpd, double &maxvpd, 
        double &rhs, double (&pd)[6], const double &gas, double &soilevap, double &baperga, double (&soilredist)[6],
        long &haltsh);
};

#endif