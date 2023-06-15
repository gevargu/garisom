/* 
# Carbon Gain vs. Hydraulic Risk Stomatal Optimization Model V.1.0

## Author: German Vargas G.
## Contact: german.vargas@utah.edu
## Date: Jun-2023

### Original model citations:
##### Original VSB code
Authors: Sperry, John S., Martin D. Venturas, William R. L. Anderegg, Maurizio Mencuccini, D. Scott Mackay, Yujie Wang, and David M. Love. 
Title: “Predicting Stomatal Responses to the Environment from the Optimization of Photosynthetic Gain and Hydraulic Cost.” 
Publication Plant, Cell & Environment 40, no. 6 (2017): 816–30. https://doi.org/10.1111/pce.12852.

##### Original cpp code 
Authors: Venturas, Martin D., John S. Sperry, David M. Love, Ethan H. Frehner, Michael G. Allred, Yujie Wang, and William R. L. Anderegg. 
Title: “A Stomatal Control Model Based on Optimization of Carbon Gain versus Hydraulic Risk Predicts Aspen Sapling Responses to Drought.” 
Publication: New Phytologist 220, no. 3 (2018): 836–50. https://doi.org/10.1111/nph.15333.

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
#include "modules/06MiscFunctions.h" // set of functions used to calculate parameter values at the beginning

// Define classes
MainProgram garisom;

int main() {
    // 1. Reading model parameters
    /* 
    - The following loop goes through each line of the parameter file and reads the values in each column. 
    - TO-DO: vectorize parameters for multi-species calculations.
    */
    garisom.cleanModelVars();
    std::cout << " ---------------------------------------------" << endl;
    std::cout << "|          READING MODEL INPUT FILES          |" << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << endl;
    std::cout << "- Configuration and parameter files: " << endl;
    garisom.readPARSheet();
    std::cout << " ---------------------------------------------" << endl;
    std::cout << "|             MODEL CONFIGURATION             |" << endl;
    std::cout << " ---------------------------------------------" << endl;
    garisom.setConfig();
    std::cout << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << "|        READING ClIMATE FORCING FILES        |" << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << endl;
    std::cout << "- Growing Season Data: " << endl;
    garisom.readGSSheet();
    std::cout << endl;
    std::cout << "- Climate forcing data: " << endl;
    garisom.readDataSheet();
    std::cout << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << "|     INITIALIZING C GAIN VS H RISK MODEL     |" << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << endl;
    std::cout << " Testing soil texture: " << garisom.paramCells[2][12] << endl;
    garisom.InitialConditions();
    std::cout << endl;
    std::cout << " Testing vg function: " << garisom.rhizor << endl;
    std::cout << endl;

    /*// seed the random number generator with something crazy
   //srand((unsigned)(time(0) * time(0)));
   // set the cout decimal precision
   std::cout.precision(12);

   // Model initialization
   long result = 0;
   std::cout << "Model Startup." << std::endl;

   // to do a normal run that's only based on local folder parameter sheet settings, set the stage_ID to zero
   mainProg.stage_ID = STAGE_ID_NONE;
   result = mainProg.modelProgramMain(); // for returning failure error codes... not really used in this version

   // Model finalization
	if (!result)
   {
      std::cout << "Model Failure! Stage " << mainProg.stage_ID << std::endl;
      return 0;
   }
   else
   {
      std::cout << "Model Success! Stage " << mainProg.stage_ID << std::endl;
   }*/
   return 1;
}