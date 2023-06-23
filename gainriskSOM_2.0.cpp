/* 
# Carbon Gain vs. Hydraulic Risk Stomatal Optimization Model V.2.0

## Author: German Vargas G.
## Contact: german.vargas@utah.edu
## Date: Jun-2023

### Model structure (file extensions)
    Garisom
        a) Input files (csv files):
            model_config
            parameter_data
            seasonlimits
            dataheader
            dataset

        b) Compilation instructions:
            Makefile
        
        c) Program code (cpp):
            gainriskSOM_2.0.0.cpp
        
        d) Modules:
            00: Main Program module
                00MainProgram (h/cpp)
            01: Input management functions 
                01CSVROW (h)
                01IOHandler (h/cpp)
            02: Soil hydrodynamics module
                02Soils (h/cpp)
            03: Plant hydraulics module
                03Hydraulics (h/cpp)
            04: Plant Morphology module
                04Morphology (h/cpp)
            05: Carbon assimilation module
                05CAssimilation (h/cpp)
            06: Miscellanious functions
                06MiscFunctions (h/cpp)

### Model outputs:

### Notes on recent updates:
    1. This version is a modular program that facilitates modifications, troubleshooting, and additions to the model.
    2. Estimation of weibull parameters via p12 and p50, which facilitates data acquisition.
    3. Also includes a cavitation fatigue module to re-calculate the weibull function while accounting for previous year damage to the xylem.
    4. It accounts for increases in atmospheric CO2 due to climate change.
    5. Unique parameter file which facilitates the inclusion of multiple species/PFTs per site
    6. Addition of unique file with model configuration
    
### Notes on future updates:
    1. Addition of a data assimilation routine to facilitate iterative forecasting. Will do using PEcAn.
    2. Addition of a hydrological competition module to quantify the effects of multiple species on the water budget.

*/ 
//
// Load the necessary libraries
//
// Cpp libraries
//
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

/* GARISOM modules
    - Each module has a set of functions.
    - All functions should be declared in the main header hile: "MainProgram.h".
    - If editing please study the documentation associated with each module.
    - Add new functions as modules.
*/ 
#include "modules/00MainProgram.h" // Main functions to run garisom
#include "modules/01CSVRow.h" // to read rows (H. Todd)
#include "modules/01IOHandler.h"// to reference array cells (H. Todd)
#include "modules/02Soils.h" // module to perform soil calculations
#include "modules/03Hydraulics.h"// module to fit weibull functions for xylem resistance
#include "modules/04Morphology.h" // module to calculate root morphological traits
#include "modules/05CAssimilation.h"// module with carbon assimilation and irradiance functions
#include "modules/06MiscFunctions.h" // set of functions used to calculate parameter values at the beginning

// Define classes
MainProgram garisom;

int main() {
    // seed the random number generator with something crazy
    srand((unsigned)(time(0) * time(0)));
    // set the cout decimal precision
    std::cout.precision(12);

    garisom.stageNames[STAGE_ID_NONE] = "standard"; //stage IDs and names are related to various modes the model was designed to run in for specific projects -- that code has been removed, but a few vestiges remain

    // model initialization

   long result = 0;   
   // to do a normal run that's only based on local folder parameter sheet settings, set the stage_ID to zero
   garisom.stage_ID = STAGE_ID_NONE;
   result = garisom.Rungarisom(); // for returning failure error codes... not really used in this version

   // Model finalization
	if (!result) {
        std::cout << "Model Failure! Stage " << garisom.stage_ID << std::endl;
        return 0;
    } else {
        std::cout << "Model Success! Stage " << garisom.stage_ID << std::endl;
    }
    return 1;
}