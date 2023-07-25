#ifndef MAINPROGRAM
#define MAINPROGRAM

#define FIO_PRECISION 12

// Necessary cpp libraries
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

// GARISOM modules
#include "01Macros.h"
#include "01CSVRow.h" // to read rows from H. Todd
#include "01IOHandler.h" // to manage inputs and outputs from H. Todd

// Modules yet to be added
#include "02Soils.h" // module to perform soil calculations
//#include "03Hydraulics.h"// module to fit weibull functions for xylem resistance
//#include "04Morphology.h" // module to calculate root morphological traits
//#include "05CAssimilation.h"// module with carbon assimilation and irradiance functions
//#include "06MiscFunctions.h" // set of functions used to calculate parameter values at the beginning

// Define classes
class MainProgram
{
public:
    // Define classess for each module
    IOHandler dSheet;
    CSVRow dataHeaderRow;
    CSVRow summaryHeaderRow;
    soils soilcalculator;

    // Variable declarations
    /* For reading data */
    double dummyDouble = 0.0;
    long rowD, colD, rowLR, colLR;
    // Configuration file
    std::string configCells[CONFIGFILE_MAXROWS][CONFIGFILE_MAXCOLS];
    // Parameter file
    std::string paramCells[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS];
    std::string nameTable[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS]; // needs to be same dimensions as the parameter file "cells"
    const long maxsp = 90;
    std::string soillayersTable[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS];
    
    // Data file
    const long maxYears = 90;
    const long bufferYears = 5; // throw out this many years when calculating optimal BA:GA to avoid the field capacity starting condition influencing the results
    double dataCells[DATAFILE_MAXROWS][DATAFILE_MAXCOLS]; // up to 2000k rows and 100 columns of input data
    double finalOutCells[MAX_SUMMARY_ROWS][MAX_SUMMARY_COLS];
    double obssolar, vpd, airtemp, maxvpd, wind, us, soiltemp;
    // Growing season data
    std::string GSCells[101][11]; // contains the growing season start/end days
    // a much simpler class for the full C++ version which just wraps access to the dataSet[][] array in the dsheet.Cells command
    // exists to preserve the code legacy of the model, which is written in VBA for Excel
    std::string GSFileName = "";

    /* For staging the model */
    long stage_ID;
    std::string stageNames[STAGE_ID_FUT_STRESS_NOACCLIM + 1];
    double stage_OptHistBAGA = 0.0;
    double stage_OptFutBAGA = 0.0;
    double stage_KmaxFut = 0.0;
    double stage_LeafResistPerFut = 0.0;
    double stage_CO2Fut = 0.0;
    double stage_Hist_refK = 0.0;
    double stage_Fut_refK = 0.0;
    
    // the minimum info necessary to identify any historical run
    std::string name_region = "";
    std::string name_site = "";
    std::string name_species = "";
    // and future
    std::string name_scen = "";
    std::string name_model = "";

    /* Model constants */
    const double pi = 3.14159;
    const double sbc = 0.0000000567; //'stefan boltzman constant in W m-2 K-4
    const double sha = 29.3; //'specific heat of air in J mol-1C-1
    const double gas = 8.3144598; //'universal gas constant J mol-1K-1
    const double oa = 0.21; //'mole fraction of o2
    const double solar = 1362; //'solar constant W m-2
                             //'Const tau = 0.65 'clear sky transmissivity, CN p. 173
    const double absolar = 0.5; //'absorptivity of solar for leaves
    const double abspar = 0.8; //'absorptivity of par for leaves
    const double absnir = 0.2; //'absorptivity of near infrared for leaves

    /* Model configuration */
    bool hysteresis, refilling, ground, soilred, sevap, raining, rainEnabled, bucket, reset_kmax,nonGS_water;
    bool nonGS_evaporation, GS_mode,useGSData, mode_predawns, iter_runSupplyCurve,iter_useAreaTable;
    bool iter_gwEnable,iter_ffcEnable, iter_bagaEnable, iter_yearsAsCount;
    int species_no;
    long iter_code; // used to override the normal iteration behavior if we want to walk back a half step
    long iter_Counter; //how many iterations of this data set have we run? For finding the supply curve
    double iter_refK;

    /* Model parameter file */
    // Site Identifiers & Parameters
    std::string region, siteID, species;
    double lat, longitude, alt, slope, slopeasp;
    // Sky Parameters
    double tau, tsncorr, emiss, ca;
    // Soils Parameters
    std::string texture;
    long layers;
    double rockfrac, rhizopercent, ffc, fieldcapfrac,soilabssol, rough, zdispl, zh, pground,grounddistance;
    // Stand-level Parameters
    double baperga, lai, xh, height, pgrav;
    // Tree-level Parameters
    double aspect, root_depth, xang, laperba, leafwidth;
    // Hydraulics
    double leafpercent , ksatp, root_p12, root_p50, stem_p12, stem_p50, leaf_p12, leaf_p50, einc;
    // Carbon Assimilation
    double vmax25, jmax25, lightcomp, qmax, kc25, ko25, comp25, thetac, havmax, hdvmax, svvmax, hajmax, hdjmax, svjmax, lightcurv;
    // BAGA Optimization
    double iter_gwDist; //we//ll keep track of the gw dist separately when running iterations, then override the grounddistance variable with this value
    long iter_ddOutMod; //output offset so we can keep multiple iterations of data
    double iter_gwInc, iter_gwStart, iter_gwEnd, iter_gwRunning;
    double iter_ffc, iter_ffcStart, iter_ffcEnd, iter_ffcInc;
    double iter_baga, iter_bagaStart, iter_bagaEnd, iter_bagaInc, iter_bagaRef, iter_bagaCutoff;

    /* Initial conditions */
    long runmean, f, k, z, i, unknowns, t;
    double cutoff, minwind, epsx, sthresh, gmax, gmaxl, patm, pinc, sum;
    // For soil parameterization
    double n[6], a[6], layerDepths[20], depthmax, vertdistance[11], depth[6], vol, radius[6];
    double length[6], shallow,kkmax[6], thetasat[6], rhizotarg, vp, vgterm, coef;
    std::vector <double> vgparams;
    // For hydraulics parameterization
    double rootpercent, stempercent, root_c, root_b, stem_c, stem_b, leaf_c, leaf_b, kmin;
    double x, rootr, stemr, leafr, plantr, ksatl, lsc, rsatp, ksats, ksatr[6], ksatroot;
    double rstem, rleaf, rplant, rhizor, kinc, kmaxr[6], kmaxrh[6], rrhizofrac;
    // Carbon Assimilation
    // Morphological traits
    double beta;

    /* Defining model variables used or to be calculated by the model */
    // Used for counting years
    bool isNewYear;
    long gs_yearIndex; /* this is a counter from 0 (for the first year) indicating how many years have passed
    - It gets the actual year from gs_ar_years(gs_yearIndex)
    - this avoids having to make year an input column in the model 
    - will just count how many years have passed when running*/
    long gs_prevDay, yearCounter;
    bool gs_inGrowSeason;
    bool gs_doneFirstDay; //done all first day of grow season calculations?
    long jd, year_cur, year_start; //not to be confused with the year array index, which is year_cur - year_start
    long yearVal; // the temporary variable where we hold the year read from the sheet
    long gs_ar_years[100], gs_ar_starts[100], gs_ar_ends[100], growSeasonCount;
    double gs_ar_Ca[100];
    // used in Pcrit calculations
    std::string failspot, layerfailure[6], tlayerfailure[6];
    long failure, layer[6], tlayer[6];
    double pcritrh[6], erh[6][100001], krh[6][100001];// rhizosphere
    double pcritr[6], kr[6][100001], er[6][100001],er_v[6][100001], kr_v[6][100001], tkr[6][100001], ter[6][100001];// roots
    double pcrits, ksh, es[100001], es_v[100001], tes[100001]; // stems
    double pcritl, el[100001], el_v[100001], tel[100001]; // leaves
    double pcritsystem;
    // E(P) curves Integration variables
    long it, tnm, j, tmax;
    double del, eps, olds, p1, p2, e, s;
    // hydraulics
    double kminroot[6], kminstem, kminleaf,kminplant,kpday1, kxday1, kplantold;
    // soil variables
    double drainage, gwflow, thetafc[6], thetafracres[6], thetafracfc[6],fc[6], theta, pd[6],prh[6], pr,
            runoff, waterold, water[6], layerflow,elayer[6][100001],laisl,laish, soilredist[6], deficit,
            swclimit[6], rain, waternew, waterchange, soilevap, transpirationtree, cinc, atree, prinitial;
    // time-step counter
    long dd, halt,haltsh, o, chalk, p, skip, check, d;
    double tod, timestep, psynmax, psynmaxsh, psynmaxmd, psynmaxshmd;
    // solar calculations
    std::string night;
    double fet, et, sm, tsn, sindec, dec, cosdec, tim, coszen, zen, cosaz,
        az, m, sp, sb, sd, st, cloud, fcd, kbe, kbezero, mleafang,
        rad, k1, t1, told, t2, kd, qd, qds, qdt, qb, qbt, qsc, qsh, qsl,
        parsh, parsl, parbottom, nirsh, nirsl, sshade, ssun, sbottom, ssunb,
        ssund, sref, par, ppfd, ea, eac, la, lg;
    // rhizosphere, root, stem and leaves flow/pressures
    double plow, elow, ehigh, klow, khigh, estart, klower, efinish, kupper, flow,
        aamax, jmatrix[7][7],vv[7], dum, indx[7], func[7], frt, dfrdpr, dfrhdprh[7],
        dfrdprh[7], dfrhdpr[7], threshold, initialthreshold, ps, pl;
    long imax, ii, ll, weird, ticks, test;
    // composite curves
    double kroot[6][100001], kleaf[100001], kstem[100001], kplant[100001],
        pstem[100001], proot[100001], pleaf[100001], eplant[100001], dedp[100001], dedpf[1000001],
        ecritsystem, prhizo[6][100001], gcmdsh;
    long total, totalv;
    // canopy pressure
    double dedplmin,pleafv[1000001], predawn, plold, dedplzero, dedpl, klossv[1000001],
        leaftshmd, lavpdshmd, maxkloss, dpmax, dpamax, amaxmax, dpamin,amaxfrac[100001],
        dpa[100001], rmean, md, dpasun, mdsh, amaxfracsh[100001], transpiration,
        psynact,gcmd, transpirationsh,psynactsh,lavpdmdsh, cincsh;
    // leaf energy balance functions
    double eplantl[1000001],numerator,denominator,leaftemp[100001],lavpd[100001],rabs,lambda,
        grad,gha,leaftempsh[100001],lavpdsh[100001], leaftmd, lavpdmd;
    // photosynthesis functions
    double gcanw[100001], gcanc[100001], comp, rday25, rday[100001], ci, jact, je, jc, var,
        psyn[100001], cin[1000001], marker, vmax, kc, ko, jmax, gcanwsh[100001], gcancsh[100001],
        rdaysh[100001], psynsh[100001], cinsh[100001], gcanwmd, emd, gcancmd, rdaymd, psynmd[100001],
        cinmd, gcanwshmd, gcancshmd,psynshmd[100001], cinshmd, rdayshmd;
    // e(p) curve resetting
    long phigh;
    // soil ground water flow and evaporation
    double store, pend, soilf[6], groundflow, emission, lc, rabssoil, mdensair, soilep, rha,
        rhs;

    // BA optimization
    double treeToPhotoLAI;

    /* Model outputs*/
    // time-step file columns
    long dColYear, dColDay, dColTime, dColSolar, dColWind, dColRain, dColTAir, dColTSoil, dColD;
    //Data column positions - NOTE THAT THEY ARE OFFSET FROM colD, the starting data column (to avoid wasting io array space)
    long dColF_p1, dColF_p2, dColF_p3, dColF_p4, dColF_p5, dColF_predawn, dColF_P, dColF_E, dColF_Gw, dColF_laVPD,
        dColF_leaftemp, dColF_ANet, dColF_s1m2, dColF_ci, dColF_PPFD, dColF_S_P, dColF_S_E, dColF_S_Gw, dColF_S_laVPD, 
        dColF_S_leaftemp, dColF_S_Anet, dColF_S_s1m2, dColF_S_ci, dColF_S_PPFD, dColF_T_E, dColF_T_ANet, dColF_T_s1m2,
        dColF_T_pcrit, dColF_T_Ecrit, dColF_CP_Pstem, dColF_CP_Proot, dColF_CP_kstem, dColF_CP_kleaf, dColF_CP_kplant,
        dColF_CP_kxylem, dColF_CP_kroot1, dColF_CP_kroot2, dColF_CP_kroot3, dColF_CP_kroot4, dColF_CP_kroot5, dColF_CP_krootAll,
        dColF_CP_Eroot1, dColF_CP_Eroot2, dColF_CP_Eroot3, dColF_CP_Eroot4, dColF_CP_Eroot5, dColF_CP_Empty1, dColF_CP_Empty2,
        dColF_End_watercontent, dColF_End_waterchange, dColF_End_rain, dColF_End_gwater, dColF_End_E, dColF_End_drainage,
        dColF_End_soilEvap, dColF_End_ET, dColF_End_ANet, dColF_End_input, dColF_End_PLCplant, dColF_End_PLCxylem,
        dColF_End_runoff;
    
    // summary file columns
    long dColF_GS_year, dColF_GS_input, dColF_GS_Anet, dColF_GS_E, dColF_GS_PLCp, dColF_GS_PLCx, dColF_GS_kPlant, dColF_GS_kXylem, dColF_GS_ET;

    // Output variables 
    /* all of these are arrays of size 100 - but this should never matter as long as 
    the gs_yearIndex never exceeds 99*/
    double gs_ar_input[100], gs_ar_Anet[100], gs_ar_E[100], gs_ar_PLCp[100], gs_ar_PLCx[100], gs_ar_kPlant[100];
    double gs_ar_kXylem[100], gs_ar_ET[100], gs_ar_PLC85[100], gs_ar_PLCSum[100], gs_ar_PLCSum_N[100];
    double gs_ar_kPlantMean[100], gs_ar_waterInitial[100], gs_ar_waterFinal[100], gs_ar_waterInitial_GS[100];
    double gs_ar_waterFinal_GS[100], gs_ar_waterInput_GS[100], gs_ar_waterInitial_OFF[100], gs_ar_waterFinal_OFF[100];
    double gs_ar_waterInput_OFF[100], gs_ar_nrFailConverge_Water[100], gs_ar_nrFailConverge_WaterMax[100], gs_ar_cica[100];
    double gs_ar_Aci[100], gs_ar_AnetDay[100];
    long gs_ar_kPlantMean_N[100], gs_ar_nrFailConverge[100], gs_ar_nrFailThreshold[100], gs_ar_cica_N[100];

    /* Functions used in the model */
    
    // For reading files
    
    void readPARSheet();                    // read the parameter file data
    void readGSSheet();                     // read growing season data and atm CO2 
    void readGrowSeasonData();              // extract cell values from growing season data and CO2
    double getCarbonByYear(const long& yearNum, std::string GSCells[101][11], const long& maxYears);//
    void readDataSheet();                   // read climate forcing data
    void setConfig();                       // sets up model configuration
    
    // For starting the model
    
    void InitialConditions();               // sets up initial conditions (substitutes ReadIn() in V 1.0.0)
    short InitModelVars();                  // initialize model variables
    short setIterationParams(long &iter);   // parameters used in the iterator
    bool isInGrowSeasonSimple();            // sets whether the model is running within the growing season
    short cleanModelVars();                 // clear values
    bool locateRanges();                    // might delete this function
    void readSiteAreaValues();              // TO-DO: incorporate this function into the model or delete
    
    // Soils module

    double get_vg(double &x);                                   // gets soil k
    void trapzdvg(double &p1, double &p2, double &s, long &t);  // integrates VG
    void qtrapvg(double &p1, double &p2, double &s);            // evaluates accuracy of van genuchten integration
    double rvg(double &x);                                      // gives soil Y in MPa
    double swc(double &x);                                      // gives soil water content
    void get_soilwetness();                                      // gets predawns accounting for rain, groundwater flow, transpiration, and redistribution
    void get_rhizoPcrit();                                      // generate soil E(P) global curve
    void get_rhizoflow();                                       // calculates rhizosphere flow using global E(P) curve
    
    // Hydraulics module

    double get_cweibull(double &p12, double &p50);                          // gets parameter c from p12 and p50
    double get_bweibull(double &p12, double &c);                            // gets parameter b from par c and p12
    double get_weibullfit(double &x, double &b, double &c, double &ksat);   // fits weibull function, used for stems and leaves
    double get_weibullfitroot(double &x);                                   // fits root-specific weibull function
    void trapzdwbr(double &p1, double &p2, double &s, long &t);             // integrates weibull function across root element z
    void qtrapwbr(double &p1, double &p2, double &s);                       // evaluates integration
    void get_rootPcrit();                                                   // generates fresh global E(P) curve for the element(erases history)
    void get_rootPcrit_v();                                                 // generates fresh global E(P) virgin curve for the element(erases history)
    void get_rootflow();                                                    // gets flow through root using global E(P) curve
    void ludcmp();                                                          // does LU decomposition on the jacobian prior to solution by lubksb
    void lubksb();                                                          // solves the decomposed jacobian delta p's
    void newtonrhapson();                                                   // returns rhizosphere pressures and root pressure, pr, as function of pd's and e
    void trapzdwbs(double &p1, double &p2, double &s, long &t);             // integrates stem element weibull
    void qtrapwbs(double &p1, double &p2, double &s);                       // evaluates stem integration
    void get_stemPcrit(bool vCurve = false);                                // generates stem E(P) curve
    void get_stem();                                                        // gets stem pressure and conductance from pr and e.
    void trapzdwbl(double &p1, double &p2, double &s, long &t);             // integrates leaf element weibull
    void qtrapwbl(double &p1, double &p2, double &s);                       // evaluates leaf integration
    void get_leafPcrit(bool vCurve = false);                                // generates fresh global E(P) curve for the leaf element(erases history)
    void get_leaf();                                                        // gets leaf pressure from stem pressure and e
    void componentpcrits();                                                 // gets critical pressures (Pcrits)
    void compositecurve();                                                  // stores composite E(P)curve and the element conductances
    void storehistory();                                                    // stores historical element e(P) curves
    void gethistory();                                                      // restores historical element e(P) curves
    void updatecurves();                                                    // resets element E(P) curves
    void canopypressure();                                                  // computes carbon-based middays

    // Photosynthesis module
    
    // For time-step simulations
    
    long modelTimestepIter(long& VBA_dd);   // time step iterator
    void modelProgramNewYear();             // new year parameterization
    void saveOutputSheet(std::string filename, std::string sheetType);
    long Rungarisom();                      // function to run the program it replaces "modelProgramMain()" (H. Todd)
};

#endif