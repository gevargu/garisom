#ifndef MAINPROGRAM
#define MAINPROGRAM

#define FIO_PRECISION 12

// For data inputs
#define MAX_SUMMARY_COLS 121
#define MAX_SUMMARY_ROWS 2001
#define PARAMFILE_MAXROWS 101
#define PARAMFILE_MAXCOLS 101
#define DATAFILE_MAXROWS 2000001
#define DATAFILE_MAXCOLS 101
#define CONFIGFILE_MAXROWS 3
#define CONFIGFILE_MAXCOLS 20

// For staging the model
#define STAGE_ID_NONE 0 // use this mode to skip all of the stage code and just run based on param sheet settings (a "normal" run, the default config)
#define STAGE_ID_HIST_OPT 1
#define STAGE_ID_HIST_STRESS 2
#define STAGE_ID_FUT_OPT 3
#define STAGE_ID_FUT_STRESS 4
#define STAGE_ID_FUT_STRESS_NOACCLIM 5

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
#include "01CSVRow.h" // to read rows from H. Todd
#include "01IOHandler.h" // to manage inputs and outputs from H. Todd
#include "02Soils.h" // module to perform soil calculations
#include "03Hydraulics.h"// module to fit weibull functions for xylem resistance
#include "04Morphology.h" // module to calculate root morphological traits
#include "06MiscFunctions.h" // set of functions used to calculate parameter values at the beginning

// Define global variables

// Define classes
class MainProgram
{
private:
    IOHandler dSheet;
    CSVRow dataHeaderRow;
    CSVRow summaryHeaderRow;
    soils soilcalculator;
    hydraulics hydraulicscalculator;
    morphology morphologycalculator;
    MiscFunctions misc_functions;

public:
    // Variable declarations
    // For reading data
    double dummyDouble = 0.0;
    long rowD, colD, rowLR, colLR;
    // Configuration file
    std::string configCells[CONFIGFILE_MAXROWS][CONFIGFILE_MAXCOLS];
    // Parameter file
    std::string paramCells[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS];
    std::string nameTable[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS]; // needs to be same dimensions as the parameter file "cells"
    const long maxsp = 90;
    // Data file
    const long maxYears = 90;
    const long bufferYears = 5; // throw out this many years when calculating optimal BA:GA to avoid the field capacity starting condition influencing the results
    double dataCells[DATAFILE_MAXROWS][DATAFILE_MAXCOLS]; // up to 2000k rows and 100 columns of input data
    double finalOutCells[MAX_SUMMARY_ROWS][MAX_SUMMARY_COLS];
    // Growing season data
    long GSCells[101][11]; // contains the growing season start/end days
    // a much simpler class for the full C++ version which just wraps access to the dataSet[][] array in the dsheet.Cells command
    // exists to preserve the code legacy of the model, which is written in VBA for Excel
    std::string GSFileName = "";

    // For staging the model
    long stage_ID;
    std::string stageNames[STAGE_ID_FUT_STRESS_NOACCLIM + 1];

    // Model constants
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

    // Model configuration
    bool hysteresis, refilling, ground, soilred, sevap, raining, bucket;
    bool reset_kmax, nonGS_water, nonGS_evaporation, GS_mode,useGSData, mode_predawns;
    int species_no;

    // Defining model variables from parameter_data.csv file
    // Site Identifiers & Parameters
    std::string region, siteID, species;
    double lat, lon, alt, slope, slope_asp;
    // Sky Parameters
    double tau, tsncorr, emiss;
    // Soils Parameters
    std::string texture;
    int layers;
    double rockfrac, rhizopercent, ffc, fieldcapfrac,soilabssol, rough, zdispl, zh, pground,grounddistance;
    // Stand-level Parameters
    double baperga, lai, xh, height;
    // Tree-level Parameters
    double aspect, root_depth, xang, laperba, leafwidth;
    // Hydraulics
    double leafpercent , ksatp, root_p12, root_p50, stem_p12, stem_p50, leaf_p12, leaf_p50;
    // Carbon Assimilation
    double vmax25, jmax25, lightcomp, qmax, kc25, ko25, comp25, thetac, havmax, hdvmax, svvmax, hajmax, hdjmax, svjmax, lightcurv;
    
    // Used in initial conditions
    long runmean, f, k, z, i, unknowns, t;
    double cutoff, minwind, epsx, sthresh, gmax, gmaxl, patm, pinc, sum;
    // For soil parameterization
    double n[6], a[6], layerDepths[20], depthmax, vertdistance[11], depth[6], vol, radius[6];
    double length[6], shallow,kkmax[6], thetasat[6], rhizotarg, vp, vgterm, coef;
    std::vector <double> vgparams;
    // For hydraulics parameterization
    double rootpercent, stempercent, root_c, root_b, stem_c, stem_b, leaf_c, leaf_b, kmin;
    double x, rootr, stemr, leafr, plantr, ksatl, lsc, rsatp, ksats, ksatr[6], ksatroot, test;
    double rstem, rleaf, rplant, rhizor, kinc, kmaxr[6], kmaxrh[6], rrhizofrac;
    // Carbon Assimilation
    // Morphological traits
    double beta;

    // Defining model variables used or to be calculated by the model
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

    // Iterator variables
    long iter_Counter; //how many iterations of this data set have we run? For finding the supply curve
    // Integration variables
    long it, tnm, j, tmax;
    double del, eps, olds;

    // Output columns
    // time-step file
    long dColYear, dColDay, dColTime, dColSolar, dColWind, dColRain, dColTAir, dColTSoil, dColD;
    //Data column positions - NOTE THAT THEY ARE OFFSET FROM colD, the starting data column (to avoid wasting io array space)
    long dColF_p1, dColF_p2, dColF_p3, dColF_p4, dColF_p5, dColF_predawn, dColF_P, dColF_E, dColF_Gw, dColF_laVPD, dColF_leaftemp, dColF_ANet,
    dColF_s1m2, dColF_ci, dColF_PPFD, dColF_S_P, dColF_S_E, dColF_S_Gw, dColF_S_laVPD, dColF_S_leaftemp,
    dColF_S_Anet, dColF_S_s1m2, dColF_S_ci, dColF_S_PPFD;
    long dColF_T_E, dColF_T_ANet, dColF_T_s1m2,
    dColF_T_pcrit, dColF_T_Ecrit, dColF_CP_Pstem, dColF_CP_Proot, dColF_CP_kstem, dColF_CP_kleaf, dColF_CP_kplant,
    dColF_CP_kxylem, dColF_CP_kroot1, dColF_CP_kroot2, dColF_CP_kroot3, dColF_CP_kroot4, dColF_CP_kroot5, dColF_CP_krootAll,
    dColF_CP_Eroot1, dColF_CP_Eroot2, dColF_CP_Eroot3, dColF_CP_Eroot4, dColF_CP_Eroot5, dColF_CP_Empty1, dColF_CP_Empty2,
    dColF_End_watercontent, dColF_End_waterchange, dColF_End_rain, dColF_End_gwater, dColF_End_E, dColF_End_drainage,
    dColF_End_soilEvap, dColF_End_ET, dColF_End_ANet, dColF_End_input, dColF_End_PLCplant, dColF_End_PLCxylem,
    dColF_End_runoff;
    // summary file
    long dColF_GS_year, dColF_GS_input, dColF_GS_Anet, dColF_GS_E, dColF_GS_PLCp, dColF_GS_PLCx, dColF_GS_kPlant, dColF_GS_kXylem, dColF_GS_ET;

    // Functions used in the model
    void readPARSheet();                    // read the parameter file data
    void readGSSheet();                     // read growing season data and atm cO2 
    void readDataSheet();                   // read climate forcing data
    void setConfig();                       // sets up model configuration
    void InitialConditions();               // sets up initial conditions (substitutes ReadIn() in V 1.0.0)
    short InitModelVars();                  // initialize model variables
    short setIterationParams(long &iter);   // parameters used in the iterator
    bool isInGrowSeasonSimple();            // sets whether the model is running within the growing season
    short cleanModelVars();                 // clear values
    bool locateRanges();                    // might delete this function
    long Rungarisom();                      // function to run the program it replaces "modelProgramMain()" (H. Todd)

};

#endif