#include "00MainProgram.h"
// The functions here go in order of implementation

/* Reading parameter data sheet*/
std::istream& operator>>(std::istream& str, CSVRow& data)
    {
    data.readNextRow(str);
    return str;
    }

void MainProgram::readPARSheet ()
{   // This function reads the model configuration and the parameter datasheet
    // This is done 1 time at the beginning of the model
    // set the model control file name
    std::string mconfigFileName = "model_config.csv";

    // reading model controls
    std::cout << std::endl;
    std::cout << "Reading model configuration from " << mconfigFileName << std::endl;

    std::ifstream configFile(mconfigFileName);
    if (!configFile.is_open()){
        std::cout << "UNRECOVERABLE: Failed to open model configuration file: " << mconfigFileName << std::endl;
        std::cout << "Model stops " << std::endl;
        std::cout << std::endl;
        std::exit(1);
    }

    long rowCount = -1;
    CSVRow row;
    while(configFile >> row){
        rowCount = rowCount +1;
        if (rowCount >= CONFIGFILE_MAXROWS){
            break;
        }

        for (int rC = 0; rC < row.size(); rC++){
            if (rowCount < CONFIGFILE_MAXROWS && rowCount >= 0 && rC < CONFIGFILE_MAXCOLS && rC >= 0){
                configCells[rowCount+1][rC + 1] = row[rC];
            }
        }
    }

    // set the parameter and nametable filenames
    std::string paramFileName = "parameters.csv";
    
    std::cout << std::endl;
    std::cout << "Reading parameters from " << paramFileName << std::endl;

    std::ifstream paramFile(paramFileName);
    if (!paramFile.is_open()) {
        std::cout << "UNRECOVERABLE: Failed to open parameter file " << paramFileName << std::endl;
        std::cout << "Model stops " << std::endl;
        std::cout << std::endl;
        std::exit(1);
    }

    // load the string values of all the parameters into an array
    // in most cases these need to be converted to doubles when called on:
    // double x = atof(stringVar.c_str());
    rowCount = -1;
    while (paramFile >> row){

        rowCount = rowCount + 1;
        if (rowCount >= PARAMFILE_MAXROWS){
            break;
        }

        for (int rC = 0; rC < row.size(); rC++){
            if (rowCount < PARAMFILE_MAXROWS && rowCount >= 0 && rC < PARAMFILE_MAXCOLS && rC >= 0){
                paramCells[rowCount + 1][rC + 1] = row[rC];
                //if (row[rC] != "" && row[rC] != "\r")
                //std::cout << "Reading parameters. " << row[rC] << std::endl;
            }
        }
        //std::cout << "4th Element(" << row[3] << ")\n";
    }

}

/* Set model configuration */
// I placed the  instructions from readIterationSettings() function in v 0.1 here
void MainProgram::setConfig(){
    // Model Configurations
    std::cout << "              Xylem hysteresis: "; 
    // turns on/off xylem hysteresis from previous growing season. Values: off; on 
    if(configCells[2][1] == "on"){
        hysteresis = true;
        std::cout << "On" << endl;
    } else {
        hysteresis = false;
        std::cout << "Off" << endl;
    }
    std::cout << "               Xylem refilling: ";
    // turns on/off xylem refilling within a growing season. Values: off; on
    if(configCells[2][2] == "on"){
        refilling = true;
        std::cout << "On" << endl;
    } else {
        refilling = false;
        std::cout << "Off" << endl;
    }
    std::cout << "                Pre-dawns mode: ";
    // turns on/off if model should consider measured pre-dawn water potential values. Values: off; on
    if(configCells[2][3] == "on"){
        mode_predawns = true;
        std::cout << "On" << endl;
        std::cout << "MODE: Running from predawn inputs (in rain column, MPa). Soil sim disabled." << std::endl;
    } else {
        mode_predawns = false;
        std::cout << "Off" << endl;
        std::cout << "MODE: Running from soil simulation (calculated predawns)." << std::endl;
    }
    std::cout << "             Ground water flow: "; 
    // turns on/off groundwater flow. Values: off; on
    if(configCells[2][4] == "on"){
        ground = true;
        std::cout << "On" << endl;
    } else {
        ground = false;
        std::cout << "Off" << endl;
    }
    std::cout << "     Soil water redistribution: ";
    // turns on/off soil redistribution routine. Values: off; on
    if(configCells[2][5] == "on"){
        soilred = true;
        std::cout << "On" << endl;
    } else {
        soilred = false;
        std::cout << "Off" << endl;
    }
    std::cout << "        Soil water evaporation: ";
    // turns on/off soil evaporation routine. Values: off; on
    if(configCells[2][6] == "on"){
        sevap = true;
        std::cout << "On" << endl;
    } else {
        sevap = false;
        std::cout << "Off" << endl;
    }
    std::cout << "                 Soil drainage: ";
    // tunrs on/off bucket option that blocks drainage through the bottom and allows reaching saturation. Values: off; on
    if(configCells[2][7] == "on"){
        bucket = true;
        std::cout << "Off" << endl;
    } else {
        bucket = false;
        std::cout << "On" << endl;
    }
    std::cout << "               Rainfall inputs: ";
    // turns on/off rain inputs. Values: off; on
    if(configCells[2][8] == "on"){
        raining = true;
        std::cout << "On" << endl;
    } else {
        raining = false;
        std::cout << "Off" << endl;
    }
    std::cout << "        Multiple years enabled: ";
    // turns on/off if we are working with multiple growing seasons. Values: off; on
    if(configCells[2][9] == "on"){
        GS_mode = true;
        std::cout << "On" << endl;
    } else {
        GS_mode = false;
        std::cout << "Off" << endl;
    }
    std::cout << "      Resetting kmax each year: ";
    // tunrs on/off reseting kmax at the beginning of each growing season. Values: off (continues the run with the accumulated losses in conductivity); on (resets kmax)
    if(configCells[2][10] == "on"){
        reset_kmax = true;
        std::cout << "On" << endl;
    } else {
        reset_kmax = false;
        std::cout << "Off" << endl;
    }
    std::cout << "          GS soil water budget: ";
    // water budget during non growing season (nonGS). Values: off (no water budget and soil filled to fraction of field capacity at beginning of new GS); on (computes water budget during nonGS)
    if(configCells[2][11] == "on"){
        nonGS_water = true;
        std::cout << "On" << endl;
    } else {
        nonGS_water = false;
        std::cout << "Off" << endl;
    }
    std::cout << " Non-GS soil water evaporation: ";
    // soil evaporation during non growing season (nonGS)  if nonGS_water == on. Values: off; on
    if(configCells[2][12] == "on"){
        nonGS_evaporation = true;
        std::cout << "On" << endl;
    } else {
        nonGS_evaporation = false;
        std::cout << "Off" << endl;
    }
    std::cout << "            Multi-Species Mode: ";
    // soil evaporation during non growing season (nonGS)  if nonGS_water == on. Values: off; on
    if(configCells[2][13] == "on"){
        species_no = std::atoi(configCells[2][14].c_str());
        std::cout << "On" << endl;
    } else {
        species_no = 1;
        std::cout << "Off" << endl;
    }

    std::cout << "        Iteration Supply Curve: ";
    // Turns on the iterations in the BAGA optimization routine. Values: off; on
    if(configCells[2][15] == "on"){
        iter_runSupplyCurve = true;
        std::cout << "On" << endl;
    } else {
        iter_runSupplyCurve = false;
        std::cout << "Off" << endl;
    }

    std::cout << "                Use Area Table: ";
    // Are we using BA:GA values from another csv file. Values: off; on
    if(configCells[2][16] == "on"){
        iter_useAreaTable = true;
        std::cout << "On" << endl;
    } else {
        iter_useAreaTable = false;
        std::cout << "Off" << endl;
    }

    std::cout << "          Iterate Ground Water: ";
    // Iterating ground water in for each stand. Values: off; on
    if(configCells[2][17] == "on"){
        iter_gwEnable = true;
        std::cout << "On" << endl;
    } else {
        iter_gwEnable = false;
        std::cout << "Off" << endl;
    }
    // we use 1 value for all species, for now...
    iter_gwInc = std::atof(paramCells[2][59].c_str()); //how much should be increment the ground water every iteration? Smaller = higher resolution supply curve
    iter_gwStart = std::atof(paramCells[2][60].c_str()); //Some data sets may not respond to increasing ground water until a high threshold, so can start with a minimum level and increment from there
    iter_gwEnd = std::atof(paramCells[2][61].c_str());

    std::cout << "  Iterate Field Capacity (FFC): ";
    // Iterating ground water in for each stand. Values: off; on
    if(configCells[2][18] == "on"){
        iter_ffcEnable = true;
        std::cout << "On" << endl;
    } else {
        iter_ffcEnable = false;
        std::cout << "Off" << endl;
    }
    // we use 1 value for all species, for now...
    iter_ffcInc = std::atof(paramCells[2][64].c_str());
    iter_ffcStart = std::atof(paramCells[2][62].c_str());
    iter_ffcEnd = std::atof(paramCells[2][63].c_str());

    std::cout << "           Iterate Field BA:GA: ";
    // Iterating BA:GA for each stand. Values: off; on
    if(configCells[2][19] == "on"){
        iter_bagaEnable = true;
        std::cout << "On" << endl;
    } else {
        iter_bagaEnable = false;
        std::cout << "Off" << endl;
    }
    // we use 1 value for all species, for now... I should move this to initialization function
    iter_bagaInc = std::atof(paramCells[2][57].c_str());
    iter_bagaStart = std::atof(paramCells[2][55].c_str());
    iter_bagaEnd = std::atof(paramCells[2][56].c_str());
    iter_bagaRef = std::atof(paramCells[2][54].c_str());
    iter_bagaCutoff = std::atof(paramCells[2][58].c_str());
    iter_bagaRef = 1.0; //156.269; // 1.0; // TODO TEMP can remove after param sheets updated
    iter_bagaEnd = 500.0; // TODO TEMP " -> for bisection method, allow extreme range

    std::cout << "       Iterate Years as Counts: ";
    // Iterating ground water in for each stand. Values: off; on
    if(configCells[2][20] == "on"){
        iter_yearsAsCount = true;
        std::cout << "On" << endl;
    } else {
        iter_yearsAsCount = false;
        std::cout << "Off" << endl;
    }
}

/* Reading growing season data sheet*/
void MainProgram::readGSSheet()
{
    memset(GSCells, 0, sizeof(GSCells));

    if (stage_ID == STAGE_ID_NONE)
    {
        GSFileName = "seasonlimits.csv";
    }

    if (GS_mode == true && GSFileName != "") // only if we actually selected a file, and we're using the GS data files in the first place
    {
		//std::cout << "Reading growing season limits." << std::endl;

        std::ifstream dataFile(GSFileName);
        if (!dataFile.is_open()){
			std::cout << "FAILED to open GS RANGE file " << GSFileName << std::endl;
		}
			
		std::cout << "Opened GS RANGE file " << GSFileName << std::endl;
			
		long yearCount = 0; // how many unique years have we encountered? 
		long curYear = 0;
		long readYear = 0;
		int rowCount = -1;
        
        CSVRow row;
        while (dataFile >> row){
			rowCount = rowCount + 1;
            if (rowCount > 100){
				break;
			}

            if (rowCount + 1 == 1){
                ; // do nothing, we don't need the header row
            } else if (rowCount + 1 > 1) {
            	readYear = std::atol(row[0].c_str());
            	if (readYear > curYear){
                	curYear = readYear;
                  	yearCount++;
                  	
					if (yearCount > maxYears) {
                     	std::cout << "Read GS range for " << yearCount - 1 << " years, quitting read." << std::endl;
                     	break;
                  	}
                  	
					std::cout << "Reading GS range year " << curYear << " (" << yearCount << "/" << maxYears << ")" << std::endl;
               	}

               for (int rC = 0; rC < row.size(); rC++) {
                	if (rowCount < 100 && rowCount >= 0 && rC < 10 && rC >= 0) {
						GSCells[rowCount + 1][rC + 1] = row[rC]; // load array as string
                        // all data i/o is in double
					}
               	}
            }
            //std::cout << "4th Element(" << row[3] << ")\n";
        }
        
		std::cout << "Finished GS RANGE read." << std::endl;
    } else {
        std::cout << "Skipping GS RANGE load -- Years will be treated as independent." << std::endl;
    }
}

void MainProgram::readGrowSeasonData(){
    long startRow = 2; // skipping the header
    long startCol = 1;

    // assumed table layout:
    /* Year  |  Start  |  End  | Ca_ppm,  with "Year" label as the named cell above */
    long gsC = -1;//growSeasonCount = gsC
    do {
        gsC = gsC + 1;
        gs_ar_years[gsC] = std::atol(GSCells[startRow + gsC][1].c_str()); // Years
        gs_ar_starts[gsC] = std::atol(GSCells[startRow + gsC][2].c_str()); // Start 
        gs_ar_ends[gsC] = std::atol(GSCells[startRow + gsC][3].c_str()); // End
        gs_ar_Ca[gsC] = std::atof(GSCells[startRow + gsC][4].c_str()); // Atmospheric CO2 (Ca_ppm)
    } while (!(std::atol(GSCells[startRow + gsC + 1][startCol].c_str()) == 0)); // hopefully we get zeros in the blank space...
    //Loop Until IsEmpty(gsRange.Offset(gsC + 1, 0).Value2)
}

/* Reading climate forcing and output data sheets*/
void MainProgram::readDataSheet()
{
    // try to find a header file
    bool foundHeaderFile = false;
    std::string headerFileName = "dataheader.csv";

    std::ifstream headerFile(headerFileName);
    if (!headerFile.is_open()){
        std::cout << "FAILED to open DATA HEADER file " << headerFileName << std::endl;
    } else {
        std::cout << "Reading DATA HEADER file " << headerFileName << std::endl;
        
        if (headerFile >> dataHeaderRow) {// should be the only row in the file
            foundHeaderFile = true;
            //std::cout << headerRow << std::endl;
        }
    }

    bool foundSumHeaderFile = false;
    headerFileName = "sumheader.csv";

    //std::cout << "Reading summary header." << std::endl;

    std::ifstream sumHeaderFile(headerFileName);
    
    if (!sumHeaderFile.is_open()){
        std::cout << "FAILED to open SUMMARY HEADER file " << headerFileName << std::endl;
    } else {
        std::cout << "Reading SUMMARY HEADER file " << headerFileName << std::endl;
        
        if (sumHeaderFile >> summaryHeaderRow){ // should be the only row in the file
            foundSumHeaderFile = true;
            //std::cout << headerRow << std::endl;
        }
    }

    // otherwise, we'll load whatever head exists in the weather file later

    memset(dataCells, 0, sizeof(dataCells));

    std::string dataFileName = "dataset.csv";

    if (stage_ID == STAGE_ID_NONE){
        dataFileName = "dataset.csv";
    }

    //std::cout << "Reading data set." << std::endl;

    std::ifstream dataFile(dataFileName);
    if (!dataFile.is_open()){
        std::cout << "FAILED to open DATA file " << dataFileName << std::endl;
    }

    std::cout << "Reading DATA file " << dataFileName << std::endl;

    // only want to read 30 years of data
    // long maxYears = 30; // change this to increase limit // moved to global
    // for now we'll take the first 30 we encounter
    long yearCount = 0; // how many unique years have we encountered? 
    long curYear = 0;
    long readYear = 0;

    int rowCount = -1;
    CSVRow row;
    while (dataFile >> row){
        rowCount = rowCount + 1;
        if (rowCount >= DATAFILE_MAXROWS){ // ~maximum size of an Excel spreadsheet
            break;
        }

        if (rowCount + 2 == rowD && !foundHeaderFile){
            // header row
            dataHeaderRow = row;
        } else if (rowCount + 2 > rowD) {// check the years
            readYear = std::atol(row[0].c_str());
            
            if (readYear > curYear){
                curYear = readYear;
                yearCount++;
                if (yearCount > maxYears){
                    std::cout << "Read in " << yearCount - 1 << " years, quitting read." << std::endl;
                    break;
                }

                std::cout << "Reading year " << curYear << " (" << yearCount << "/" << maxYears << ")" << std::endl;
            }

            for (int rC = 0; rC < row.size(); rC++){
                if (rowCount < DATAFILE_MAXROWS && rowCount >= 0 && rC < DATAFILE_MAXCOLS && rC >= 0){
                    dataCells[rowCount + 1][rC + 1] = atof(row[rC].c_str()); // load array as double
                } else { // all data i/o is in double
                    long breakpoint = 1137;
                }
            }
        }
        //std::cout << "4th Element(" << row[3] << ")\n";
    }
    
    std::cout << "Finished DATA read." << std::endl;
}

/* Initialize parameter values*/
void MainProgram::InitialConditions()
{
    /* 
    - The following loop will allow to extract parameters for each species/PFT in a single 
    location or new parameters.
    - For 1 species/PFT there is only 1 line
    TO-DO: vectorize the variables or use an array for storing values
    */
    for (int sp = 0; sp <= species_no ; sp++) {
        if(sp + 1 <= 1){// + 1 helps as the loop starts in 0
            // do nothing, we skip the first row
        } else {
            // Site Identifiers & Parameters
            region = paramCells[sp+1][1];                                // name of the region where the simulations take place
            name_region = region;
            siteID = paramCells[sp+1][2];                                // site/simulation ID. Usually the plot/stand identifier  
            name_site = siteID;         
            lat = std::atof(paramCells[sp+1][3].c_str());                // latitude in degree fraction north
            lon = std::atof(paramCells[sp+1][4].c_str());                // longitude in degree fraction west
            species = paramCells[sp+1][5];                               // species or PFT represented in parameter data
            name_species = species;                             
            alt = std::atof(paramCells[sp+1][6].c_str());                // site elevation in m above sea level
            slope = std::atof(paramCells[sp+1][7].c_str());              // slope inclination; degrees from horizontal
            slope_asp = std::atof(paramCells[sp+1][8].c_str());          // slope aspect; counterclockwise degrees from south
            
            // Sky Parameters
            tau = std::atof(paramCells[sp+1][9].c_str());                // atmospheric transmittance from weather data (set to 0.65 as default if no data available)
            tsncorr = std::atof(paramCells[sp+1][10].c_str());           // solar noon correction from weather data in hours
            emiss = std::atof(paramCells[sp+1][11].c_str());             // long wave emissivity
            
            // Soil Parameters                      
            texture = paramCells[sp+1][12];                              // USDA soil texture category (equal for all layers but could also be determined per layer)
            if(mode_predawns == true){ // one soil layer for pre-dawn configuration
                layers = 1;
            } else {
                layers = std::atoi(paramCells[sp+1][13].c_str());        // number of soil layers (select 1-n)
            }
            rockfrac = std::atof(paramCells[sp+1][14].c_str());          // fraction of soil volume as rocks (0-1)
            rockfrac = 1 - rockfrac;                                     // fraction of volume with no rocks
            rhizopercent = std::atof(paramCells[sp+1][15].c_str());      // average percent of whole plant resistance in rhizosphere (maximum soil limitation)
            ffc = std::atof(paramCells[sp+1][16].c_str()) / 100.0;       // fraction of field capacity for starting the season
            fieldcapfrac = std::atof(paramCells[sp+1][17].c_str());      // fraction that field capacity is of saturation (minus residual)
            soilabssol = std::atof(paramCells[sp+1][18].c_str());        // absorptivity of soil surface for solar
            pground = std::atof(paramCells[sp+1][19].c_str());           // ground water pressure
            grounddistance = std::atof(paramCells[sp+1][20].c_str());    // distance to ground water source in m from the bottom of the rootsystem
            rough = 0.01;                                                // soil Zm, eqn 14.9, using table 5.1 for smooth surface, cm
            zdispl = 6.5 * rough;                                        // soil d, eqn 14.9, using d = 6.5 Zm, eq 5.2,5.3
            zh = 0.2 * rough;                                            // roughness for temperature, important when simulating rain
            
            // Stand-level Parameters
            baperga = std::atof(paramCells[sp+1][21].c_str()) * 0.0001;  // basal area per ground area m2 ha-1
            lai = std::atof(paramCells[sp+1][22].c_str());               // canopy lai
            xh = std::atof(paramCells[sp+1][23].c_str());                // height above soil surface for understory wind and gh in m
            height = std::atof(paramCells[sp+1][24].c_str());            // average tree height in m
            
            // Tree-level Parameters
            aspect = std::atof(paramCells[sp+1][25].c_str());            // max radius of root system per max depth
            root_depth = std::atof(paramCells[sp+1][26].c_str());        // maximum rooting depth in m
            xang = std::atof(paramCells[sp+1][27].c_str());              // leaf angle parameter; CN 15.4
            laperba = std::atof(paramCells[sp+1][28].c_str());           // initial leaf area per basal area per individual tree; m2 m-2
            leafwidth = std::atof(paramCells[sp+1][29].c_str()) * 0.72;  // leaf width in m x factor = characteristic dimension(campbell and norman)
            
            // Hydraulics
            leafpercent = std::atof(paramCells[sp+1][30].c_str());       // saturated % of tree resistance in leaves
            ksatp = std::atof(paramCells[sp+1][31].c_str());             // kmax of tree in kg hr-1 m-2 MPa-1 per basal area
            root_p12 = std::atof(paramCells[sp+1][32].c_str());          // root element p12
            root_p50 = std::atof(paramCells[sp+1][33].c_str());          // root element p50
            stem_p12 = std::atof(paramCells[sp+1][34].c_str());          // stem p12
            stem_p50 = std::atof(paramCells[sp+1][35].c_str());          // stem p50
            leaf_p12 = std::atof(paramCells[sp+1][36].c_str());          // leaf p12
            leaf_p50 = std::atof(paramCells[sp+1][37].c_str());          // stem p50
            pinc = std::atof(paramCells[sp+1][38].c_str());              // pressure increment for curve calculation
            if (pinc <= 0.0){
                pinc = 0.00075; 
            } else {
            // do nothing pinc already there
            }

            // Carbon Assimilation
            vmax25 = std::atof(paramCells[sp+1][39].c_str());            // umol m-2 s-1; maximum carboxylation rate (vmax) at 25C
            jmax25 = std::atof(paramCells[sp+1][40].c_str());            // umol m-2 s-1; maximum electron transport rate (jmax) at 25C
            lightcomp = std::atof(paramCells[sp+1][41].c_str());         // light compensation point in ppfd
            qmax = std::atof(paramCells[sp+1][42].c_str());              // quantum yield of electron transport; moles e per mols photons
            kc25 = std::atof(paramCells[sp+1][43].c_str());              // m-m constant for CO2 in mole fraction at 25C. Bernacchi T response
            ko25 = std::atof(paramCells[sp+1][44].c_str());              // m-m constant for O2 in mole fraction at 25C. Bernacchi T response
            comp25 = std::atof(paramCells[sp+1][45].c_str());            // photorespiratory compensation point in mole fraction at 25C. Bernacchi T response
            thetac = std::atof(paramCells[sp+1][46].c_str());            // shape factor for A-ci colimitation
            havmax = std::atof(paramCells[sp+1][47].c_str());            // J mol-1; temp-dependency parameters from Leunig 2002
            hdvmax = std::atof(paramCells[sp+1][48].c_str());            // J mol-1; temp-dependency parameters from Leunig 2002
            svvmax = std::atof(paramCells[sp+1][49].c_str());            // J mol-1 K-1; temp-dependency parameters from Leunig 2002
            hajmax = std::atof(paramCells[sp+1][50].c_str());            // J mol-1; temp-dependency parameters from Leunig 2002
            hdjmax = std::atof(paramCells[sp+1][51].c_str());            // J mol-1; temp-dependency parameters from Leunig 2002
            svjmax = std::atof(paramCells[sp+1][52].c_str());            // J mol-1 K-1; temp-dependency parameters from Leunig 2002
            lightcurv = std::atof(paramCells[sp+1][53].c_str());         // temp-dependency parameters from Leunig 2002
        
        if (stage_ID == STAGE_ID_NONE || stage_ID == STAGE_ID_HIST_STRESS || stage_ID == STAGE_ID_FUT_STRESS || stage_ID == STAGE_ID_FUT_STRESS_NOACCLIM) {
            if (GS_mode == 1) {
                useGSData = true;
            } else {
                useGSData = false;
            } // endif
        } else if (stage_ID == STAGE_ID_HIST_OPT || stage_ID == STAGE_ID_FUT_OPT) {
            if (GS_mode == 1) { // different parameter name
                useGSData = true;
            } else {
                useGSData = false;
            } // endif
        } else{ // impossible unknown stage failsafe
            useGSData = false;
        }
        
        // 2. Calculating initial conditions
        /* 
        - This section calculate the variables used in the first iteration. 
        - It is within the species loop as in the future can be easily modified to calculate everything as vector of species or pfts
        - TO-DO: vectorize values for multiple species
        */
        std::cout << "----------- Setting up model parameters. (year count = " << gs_yearIndex << ")" << std::endl;
        std::cout << endl;
   
        /* Site initial conditions */
        runmean = 1; //Cells(6, 9) //running mean for profit maximization
        cutoff = 1.1; //Cells(7, 9) //cutoff for stopping dpamax search
        minwind = 0.4515; //m s-1'minimum wind threshold
        f = 70; //itmax for trapzd routine used for xylem only
        //atmospheric pressure, T = 15 C, average sealevel patm; approximation
        //Call setValueFromName("o_atmP", patm) //Cells(1, 18) = patm //atmospheric pressure in kPa
        patm = misc_functions.get_patm(alt);

        /* Plant initial conditions */
        epsx = 0.0001; //fractional acceptable error for trapzd routine in xylem
        sthresh = 1.0001; //acceptable error for e integral for xylem
        gmax = 1000000; //Cells(2, 13) //maximum G, wet soil, vpd = 0, in kg m - 2 hr - 1 basal area
        rsatp = misc_functions.get_resist(ksatp); // convert to resistance
        gmaxl = hydraulicscalculator.get_gmaxl(gmax, laperba); //convert to gmax per leaf area in mmol m-2s-1
        rootpercent = hydraulicscalculator.get_rootpercent(leafpercent); // saturated % of tree resistance in roots
        stempercent = hydraulicscalculator.get_rootpercent(leafpercent); // saturated % of tree resistance in stems
        ksatl = hydraulicscalculator.get_leafkmax(ksatp, leafpercent); //leaf conductance per basal area
        lsc = hydraulicscalculator.get_LSC(ksatl, laperba); //lsc in kg hr-1m-2MPa-1
        ksatroot = hydraulicscalculator.get_rootkmax(rootpercent, rsatp); //kmax of root system; assumes zero % rhizosphere resistance in WET soil
        ksats = hydraulicscalculator.get_stemkmax(stempercent, rsatp); //kmax of stem system
        kmin = ksatp / 2000.0; //"instantaneous K" cutoff for global K(P) curves for each element
        beta = morphologycalculator.get_rootbeta(root_depth); // amount of root biomass above the max rooting depth
        
        /* Soil & Hydraulic initial conditions 
        Here we work with vectors because of the multi-layer characteristics of the soil
        They are calculated after the plant initial hydraulic conditions as they depend on root properties
        TO-DO: incorporate the layer loops into the functions
        */
        for (k = 1; k <= layers; k++)
        {// set layer depths and % root ksat
            layerDepths[k] = soilcalculator.get_depths(k, layers, beta);
        }
        depthmax = layerDepths[layers];// get maximum depth in m

        // calculate transport distances
        for (k = 1; k <= layers * 2.0; k++)
        {
            vertdistance[k] = soilcalculator.get_halfdepths(k, layers, beta); //get half depths
        }
        
        i = 0;
        for (k = 1; k <= layers * 2; k += 2) // To layers * 2 Step 2
        {
            i = i + 1;
            vertdistance[i] = vertdistance[k]; //take every other vertdistance
        }
        
        //now get radial distances
        for (k = 1; k <= layers; k++) //To layers //get thicknesses
        {
            if (k == 1)
            {
                depth[k] = layerDepths[k]; //Cells(8 + k, 11)
            } else {
                depth[k] = layerDepths[k] - layerDepths[k - 1]; //lSheet.Cells(rowLR + k, colLR + 11) - lSheet.Cells(rowLR + k - 1, colLR + 11)  //depth is layer thickness in meters
                } // endif
        }
        vol = soilcalculator.get_layvolume(depth[1],pi, depthmax, aspect);//volume of first, and hence all, layers
        
        //get radial widths of each layer and transport length
        for (k = 1; k <= layers; k++) //To layers
        {
            radius[k] = soilcalculator.get_radius(vol, depth[k], pi); //width in m
            length[k] = radius[k] + vertdistance[k]; //transport distance
            if (k == 1)
            {
                shallow = length[k];
                length[k] = length[k] / shallow; //normalize to shallowest layer
            }
        }
        
        unknowns = layers + 1; //number of unknowns to be solved for; also dimensions of matrix
        depth[0] = 0; //depth of surface
        // get Van Genuchten parameters and propagate for each layer
        // if different layers have different textures, include the function inside
        vgparams = soilcalculator.get_vgparams(texture);
        for (k = 1; k <= layers; k++)
        {
            a[k] = vgparams[0];//van genuchten alpha
            n[k] = vgparams[1];//van genuchten n
            kkmax[k] = vgparams[2];//saturated conductivity of soil in kg hr-1MPa-1 m-1
            thetasat[k] = vgparams[3]; //theta sat in volume/volume
            thetasat[k] = thetasat[k] * rockfrac; //reduce for actual rock-free fraction of soil
        }

        //now add toplayer (layer 0) of rootless soil 2 cm thick w. same properties as layer 1
        a[0] = a[1];
        n[0] = n[1];
        kkmax[0] = kkmax[1];
        thetasat[0] = thetasat[1];
        depth[0] = 0.02; //sets top layer to 2 cm
        //now solve for kmax rhizosphere that gives the desired ave % rhizosphere resistance
        rhizotarg = rhizopercent/100;
        z = 1; //use layer 1 as stand in for whole root system
        ksatr[1] = ksatroot; //set to whole root system
        x = 0.5; //start by finding kmaxrh at 0.5 MPa...a deliberate under-shoot
        // get resistances along plant components
        // root
        root_c = hydraulicscalculator.get_cweibull(root_p12,root_p50);
        root_b = hydraulicscalculator.get_bweibull(root_p12,root_c);
        rootr = 1.0/ hydraulicscalculator.get_weibullfitroot(x, root_b, root_c,ksatr,z);
        // stem
        stem_c = hydraulicscalculator.get_cweibull(stem_p12,stem_p50);
        stem_b = hydraulicscalculator.get_bweibull(stem_p12,stem_c);
        rstem = 1.0/ hydraulicscalculator.get_weibullfit(x, stem_b, stem_c, ksats);
        // leaves
        leaf_c = hydraulicscalculator.get_cweibull(leaf_p12,leaf_p50);
        leaf_b = hydraulicscalculator.get_bweibull(leaf_p12,leaf_c);
        rleaf = 1.0/ hydraulicscalculator.get_weibullfit(x, leaf_b, leaf_c, ksatl);
        // whole plant
        rplant = rootr + rstem + rleaf; //rplant here is just the xylem part
        // rhizosphere
        rhizor = soilcalculator.get_rhizor(rplant, rhizotarg); // Rhizosphere resistance
        vp = soilcalculator.get_vgp(a,n,x,z); //van genuchten terms 
        vgterm = soilcalculator.get_vgterm(n, vp, z); //van genuchten terms // vgterm = vp ^ ((n[z] - 1) / (2 * n[z])) * ((1 - vp) ^ ((n[z] - 1) / n[z]) - 1) ^ 2
        kmaxrh[1] = hydraulicscalculator.get_kmaxrh(rhizor, vgterm); //solve for kmaxrh[1]
        kinc = kmaxrh[1] * 0.1;

        do // loop through rhizosphere kmax 
        {
            kmaxrh[1] = kmaxrh[1] + kinc;
            x = 0;
            sum = 0;
            do // loop through pressures
            {
                x = x + 0.1;
                rootr = 1.0 / hydraulicscalculator.get_weibullfitroot(x, root_b, root_c,ksatr,z);
                rstem = 1.0 / hydraulicscalculator.get_weibullfit(x, stem_b, stem_c, ksats);
                rleaf = 1.0 / hydraulicscalculator.get_weibullfit(x, leaf_b, leaf_c, ksatl);
                rhizor = 1.0 / soilcalculator.get_vg(a,n,x,kmaxrh,z);
                rplant = rootr + rstem + rleaf + rhizor;
                rrhizofrac = rhizor / rplant; //fraction of resistance in rhizosphere
                sum = sum + rrhizofrac; //add up fractions
            } while (!(1.0 / rplant < kmin)); //Loop Until 1 / rplant < kmin //average over full range
            sum = sum / (x / 0.1); //average fraction
        } while (!(sum < rhizotarg)); // Until sum < rhizotarg //loop until desired soil limitation is reached

        kmaxrh[1] = kmaxrh[1] / layers; //divide whole root rhizokmax into equal portions for each layer
                                        //end of soil limitation adjustment
        //now set soil layer parameters based on aroot of entire root system
        for (z = 1; z <= layers; z++) //z = 1 To layers
        {
            kmaxrh[z] = kmaxrh[1]; //soil to root MAXIMUM conductance in kg hr-1 MPa-1//re - set for individual layers
            //paramCells[rowLR + z][colLR + 7] = std::to_string(kmaxrh[z]); //note: this is kMAX...at P=0; not at saturated PD
        }
        t = 0;
        // loop to find ksatroot for each layer
        coef = 0.0;
        do
        {
            coef = coef + 0.01;
            sum = 0.0;
            for (k = 1; k <= layers; k++)//k = 1 To layers //soil layers from top to bottom
            {
                ksatr[k] = coef / length[k]; //assumes ksatr proportional to biomass/length
                sum = sum + ksatr[k];
            }
        } while (!(sum > ksatroot)); //Loop Until sum > ksatroot //loop until each layer adds to total

        //for (k = 1; k <= layers; k++)//k = 1 To layers
        //{
            //paramCells[rowLR + k][colLR + 2] = std::to_string(100.0 * (ksatr[k] / ksatroot));
            //paramCells[rowLR + k][colLR + 9] = std::to_string(ksatr[k]);
        //}
        // ca = ca * 0.000001; //ambient co2 in moles per mole
        //Cells(1, 21) = ca * patm * 1000 //ambient co2 in Pa
        //Call setValueFromName("o_co2AmbPa", ca * patm * 1000) // this is a pure output from the Excel version, OK to disable here
        // the value in this cell will not be referenced again by the model
        // Will replace when loading the GS data
            
            if(mode_predawns == true){ // one soil layer for pre-dawn configuration
                layers = 1;
            } else {}
            if (stage_ID == STAGE_ID_NONE || stage_ID == STAGE_ID_HIST_STRESS || stage_ID == STAGE_ID_FUT_STRESS || stage_ID == STAGE_ID_FUT_STRESS_NOACCLIM) {
                if (GS_mode == true) {
                    useGSData = true;
                } else {
                    useGSData = false;
                } // endif
            } else if (stage_ID == STAGE_ID_HIST_OPT || stage_ID == STAGE_ID_FUT_OPT) {
                if (GS_mode == 1) { // different parameter name
                    useGSData = true;
                } else {
                    useGSData = false;
                } // endif
            } else{ // impossible unknown stage failsafe
                useGSData = false;
            }
        }
    }
}

/* Initial model variables*/
short MainProgram::InitModelVars(){
    // this can be called on the start of every iteration BUT NOT ON NEW YEARS
    isNewYear = true;
    // this is a good place for general initialization too
    gs_yearIndex = 0;
    gs_prevDay = 0;
    gs_inGrowSeason = false;
    gs_doneFirstDay = false;
    year_cur = 0;
    year_start = 0;
    yearVal = 0;

    return -1;
}

/* Cleaning parameter values*/
short MainProgram::cleanModelVars()
{
    // model program vars
    root_b = 0;
    root_p12 = 0;
}

/* Model iterator */
short MainProgram::setIterationParams(long &iter){
    if (iter < 1000 && iter > -1){
        //gs_yearIndex = iter;
        iter_Counter = iter;
        return -1; // true in short language
    } else {
        return 0;
    }
}

/* Iteration settings*/
// It seems this function was deprecated in version 0.1 
// parameters still used: i_iter_bagaInc 
// All these settings should come from the model control now

/* Growing season boolean */
bool MainProgram::isInGrowSeasonSimple(){
    if (useGSData){
        if (gs_yearIndex >= 0 && gs_yearIndex < 100){
            if (gs_ar_years[gs_yearIndex] > 0){
                if (jd >= gs_ar_starts[gs_yearIndex] && jd <= gs_ar_ends[gs_yearIndex]){
                  return true;
               }
            }
        }
        return false;
    } else {
        return true; // if the GS limits are disabled, we're always in the growing season
    }
}

/* Range locator */
bool MainProgram::locateRanges() {
      // an ugly hack to setup the column order
      dColYear = 1;
      dColDay = dColYear + 1;
      dColTime = dColDay + 1;
      dColSolar = dColTime + 1;
      dColRain = dColSolar + 1;
      dColWind = dColRain + 1;
      dColTAir = dColWind + 1;
      dColTSoil = dColTAir + 1;
      dColD = dColTSoil + 1;
      dColF_p1 = dColD + 1;
      dColF_p2 = dColF_p1 + 1;
      dColF_p3 = dColF_p2 + 1;
      dColF_p4 = dColF_p3 + 1;
      dColF_p5 = dColF_p4 + 1;
      dColF_predawn = dColF_p5 + 1;
      dColF_P = dColF_predawn + 1;
      dColF_E = dColF_P + 1;
      dColF_Gw = dColF_E + 1;
      dColF_laVPD = dColF_Gw + 1;
      dColF_leaftemp = dColF_laVPD + 1;
      dColF_ANet = dColF_leaftemp + 1;
      dColF_s1m2 = dColF_ANet + 1;
      dColF_ci = dColF_s1m2 + 1;
      dColF_PPFD = dColF_ci + 1;
      dColF_S_P = dColF_PPFD + 1;
      dColF_S_E = dColF_S_P + 1;
      dColF_S_Gw = dColF_S_E + 1;
      dColF_S_laVPD = dColF_S_Gw + 1;
      dColF_S_leaftemp = dColF_S_laVPD + 1;
      dColF_S_Anet = dColF_S_leaftemp + 1;
      dColF_S_s1m2 = dColF_S_Anet + 1;
      dColF_S_ci = dColF_S_s1m2 + 1;
      dColF_S_PPFD = dColF_S_ci + 1;
      dColF_T_E = dColF_S_PPFD + 1;
      dColF_T_ANet = dColF_T_E + 1;
      dColF_T_s1m2 = dColF_T_ANet + 1;
      dColF_T_pcrit = dColF_T_s1m2 + 1;
      dColF_T_Ecrit = dColF_T_pcrit + 1;
      dColF_CP_Pstem = dColF_T_Ecrit + 1;
      dColF_CP_Proot = dColF_CP_Pstem + 1;
      dColF_CP_kstem = dColF_CP_Proot + 1;
      dColF_CP_kleaf = dColF_CP_kstem + 1;
      dColF_CP_kplant = dColF_CP_kleaf + 1;
      dColF_CP_kxylem = dColF_CP_kplant + 1;
      dColF_CP_kroot1 = dColF_CP_kxylem + 1;
      dColF_CP_kroot2 = dColF_CP_kroot1 + 1;
      dColF_CP_kroot3 = dColF_CP_kroot2 + 1;
      dColF_CP_kroot4 = dColF_CP_kroot3 + 1;
      dColF_CP_kroot5 = dColF_CP_kroot4 + 1;
      dColF_CP_krootAll = dColF_CP_kroot5 + 1;
      dColF_CP_Eroot1 = dColF_CP_krootAll + 1;
      dColF_CP_Eroot2 = dColF_CP_Eroot1 + 1;
      dColF_CP_Eroot3 = dColF_CP_Eroot2 + 1;
      dColF_CP_Eroot4 = dColF_CP_Eroot3 + 1;
      dColF_CP_Eroot5 = dColF_CP_Eroot4 + 1;
      dColF_CP_Empty1 = dColF_CP_Eroot5 + 1;
      dColF_CP_Empty2 = dColF_CP_Empty1 + 1;
      dColF_End_watercontent = dColF_CP_Empty2 + 1;
      dColF_End_waterchange = dColF_End_watercontent + 1;
      dColF_End_rain = dColF_End_waterchange + 1;
      dColF_End_gwater = dColF_End_rain + 1;
      dColF_End_E = dColF_End_gwater + 1;
      dColF_End_drainage = dColF_End_E + 1;
      dColF_End_soilEvap = dColF_End_drainage + 1;
      dColF_End_ET = dColF_End_soilEvap + 1;
      dColF_End_ANet = dColF_End_ET + 1;
      dColF_End_input = dColF_End_ANet + 1;
      dColF_End_PLCplant = dColF_End_input + 1;
      dColF_End_PLCxylem = dColF_End_PLCplant + 1;
      dColF_End_runoff = dColF_End_PLCxylem + 1;

      dColF_GS_year = 1;//dColF_End_runoff + 1;
      dColF_GS_input = dColF_GS_year + 1;
      dColF_GS_Anet = dColF_GS_input + 1;
      dColF_GS_E = dColF_GS_Anet + 1;
      dColF_GS_PLCp = dColF_GS_E + 1;
      dColF_GS_PLCx = dColF_GS_PLCp + 1;
      dColF_GS_kPlant = dColF_GS_PLCx + 1;
      dColF_GS_kXylem = dColF_GS_kPlant + 1;
      dColF_GS_ET = dColF_GS_kXylem + 1;

      // TODO HACK get the region, site, spec names now -- these aren't in the current nameTable files so we'll do it manually
      // fix this to use the nameTable when we have time to re-run all the acclimations... or write a script to update all the nametables
      // note: this is specific to the BA/GA optimization

      //std::cout << "Read site and scenario names (region | site | species | scenario | model): " << name_region << " | " << name_site << " | " << name_species << " | " << name_scen << " | " << name_model << std::endl;

    //
    return true;
}

/* Read site area values
TO-DO: incorporate this function into the model or delete
*/
void MainProgram::readSiteAreaValues(){
   }

// /* Get pcrits */
// void MainProgram::componentpcrits(){ //'gets pcrits
//     soilcalculator.rhizocurves(z,layers,p1,p2,k,e,kmaxrh,pinc,s,x,a,n,kmin);
//     //rootcurves(); //'gets root element curves
//     //stemcurve(); //'gets stem element curve
//     //pcrits;
//     //leafcurve(); //'gets leaf element curve
//     //pcritl;
    
//     //rootcurves_v(); //erases history for md solution
//     //stemcurve(true);
//     //leafcurve(true);

//     //memset(ter, 0, sizeof(ter));
//     //memset(tkr, 0, sizeof(tkr));
//     //memset(tes, 0, sizeof(tes));
//     //memset(tel, 0, sizeof(tel));
// }

/* Running the model */ 
long MainProgram::Rungarisom(){
    std::cout << " ---------------------------------------------" << endl;
    std::cout << "|  CARBON GAIN VS HYDRAULIC RISK MODEL V 1.0  |" << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << endl;

    //Dim ddOutMod As Long // moved to module global
    memset(finalOutCells, 0, sizeof(finalOutCells)); // clear the final outputs data. Normal data sheets get cleared on each iteration, but this one only per-run
    
    bool lrSuccess = MainProgram::locateRanges(); //Finds the output sections across the workbook
    if (!lrSuccess) {
        std::cout << "Unrecoverable model failure!" << std::endl;
        std::cout << "Model stops " << std::endl;
        std::cout << std::endl;
        std::exit(1);
        //return 0; // failure, unrecoverable
    }

    MainProgram::cleanModelVars();
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << "|          READING MODEL INPUT FILES          |" << std::endl;
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "- Configuration and parameter files: " << std::endl;
    MainProgram::readPARSheet();
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << "|             MODEL CONFIGURATION             |" << std::endl;
    std::cout << " ---------------------------------------------" << std::endl;
    MainProgram::setConfig();
    std::cout << std::endl;
    iter_ddOutMod = 0;
    iter_Counter = 0;
    iter_code = 0;

    MainProgram::InitModelVars();
    MainProgram::InitialConditions();
    
    if (stage_ID == STAGE_ID_FUT_STRESS_NOACCLIM){ // override some if we're doing the odd "no acclimation stress profile"
        std::cout << "Stage " << stage_ID << "; NoAcclim Stress Profile, overriding historical ca " << ca << " -> " << stage_CO2Fut << " and ksatp " << ksatp << " -> " << stage_KmaxFut << std::endl;
    }
    std::cout << " ---------------------------------------------" << endl;
    std::cout << "|        READING ClIMATE FORCING FILES        |" << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << endl;
    std::cout << "- Growing Season Data: " << endl;
    MainProgram::readGSSheet();
    MainProgram::readGrowSeasonData();
    std::cout << endl;  
    std::cout << "- Climate forcing data: " << endl;
    MainProgram::readDataSheet();
    std::cout << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << "|             INITIALIZING MODEL              |" << endl;
    std::cout << " ---------------------------------------------" << endl;
    std::cout << endl;

    if((iter_useAreaTable)){
        MainProgram::readSiteAreaValues();
    }

    if (iter_Counter == 0){ //we//re on the first iteration (or we//re not using iterations)
        gs_yearIndex = 0; // multiyear
        gs_prevDay = 0;
        gs_inGrowSeason = false;
    } else {
        //things to do ONLY if we//re NOT on the first iteration
    } // End If
    
    std::memset(gs_ar_input, 0, sizeof(gs_ar_input));
    std::memset(gs_ar_Anet, 0, sizeof(gs_ar_Anet));
    std::memset(gs_ar_E, 0, sizeof(gs_ar_E));
    std::memset(gs_ar_PLCp, 0, sizeof(gs_ar_PLCp));
    std::memset(gs_ar_PLCx, 0, sizeof(gs_ar_PLCx));
    std::memset(gs_ar_kPlant, 0, sizeof(gs_ar_kPlant));
    std::memset(gs_ar_kXylem, 0, sizeof(gs_ar_kXylem));
    std::memset(gs_ar_ET, 0, sizeof(gs_ar_ET));
    std::memset(gs_ar_PLC85, 0, sizeof(gs_ar_PLC85));
    std::memset(gs_ar_PLCSum, 0, sizeof(gs_ar_PLCSum));
    std::memset(gs_ar_PLCSum_N, 0, sizeof(gs_ar_PLCSum_N));

    std::memset(gs_ar_kPlantMean, 0, sizeof(gs_ar_kPlantMean));
    std::memset(gs_ar_kPlantMean_N, 0, sizeof(gs_ar_kPlantMean_N));
    std::memset(gs_ar_waterInitial, 0, sizeof(gs_ar_waterInitial));
    std::memset(gs_ar_waterFinal, 0, sizeof(gs_ar_waterFinal));

    std::memset(gs_ar_waterInitial_GS, 0, sizeof(gs_ar_waterInitial_GS));
    std::memset(gs_ar_waterFinal_GS, 0, sizeof(gs_ar_waterFinal_GS));
    std::memset(gs_ar_waterInput_GS, 0, sizeof(gs_ar_waterInput_GS));

    std::memset(gs_ar_waterInitial_OFF, 0, sizeof(gs_ar_waterInitial_OFF));
    std::memset(gs_ar_waterFinal_OFF, 0, sizeof(gs_ar_waterFinal_OFF));
    std::memset(gs_ar_waterInput_OFF, 0, sizeof(gs_ar_waterInput_OFF));

    std::memset(gs_ar_nrFailConverge, 0, sizeof(gs_ar_nrFailConverge));
    std::memset(gs_ar_nrFailConverge_Water, 0, sizeof(gs_ar_nrFailConverge_Water));
    std::memset(gs_ar_nrFailThreshold, 0, sizeof(gs_ar_nrFailThreshold));

    std::memset(gs_ar_cica, 0, sizeof(gs_ar_cica));
    std::memset(gs_ar_cica_N, 0, sizeof(gs_ar_cica_N));

    std::memset(gs_ar_Aci, 0, sizeof(gs_ar_Aci));
    std::memset(gs_ar_AnetDay, 0, sizeof(gs_ar_AnetDay));

    for (k = 0; k <= layers; k++){ // k = 0 To layers //assign source pressures, set layer participation
        layerfailure[k] = "ok";
        layer[k] = 0; //1 if out of function
    } // Next k

    failure = 0; //=1 for system failure at midday...terminates run
    failspot = "no failure";
    //componentpcrits(); //gets pcrits for each component
    failspot = "no failure";
    std::cout << std::endl;
    std::cout << "Rhizosphere Pcrit for layer 1: " << pcritrh[1] << std::endl;
    std::cout << std::endl;
    std::cout << "------------------ TESTING AREA ------------------"<< std::endl; 
    std::cout << "Rhizosphere Resistance: " << rhizor << std::endl;
   // soilcalculator.rhizor_change(&rhizor, kmaxrh[6]);
    std::cout << std::endl;
    std::cout << "Rhizosphere Resistance: " << rhizor << std::endl;
    std::cout << std::endl;
    for (k = 1; k <= layers; k++){ // k = 1 To layers //exclude the top layer
        //kminroot[k] = ksatr[k];
    } // Next k

    //kminstem = ksats;
    //kminleaf = ksatl;
    //kminplant = ksatp;

    //gwflow = 0; //inflow to bottom of root zone
    //drainage = 0; //drainage from bottom of root zone

    //dd = 0;
    long ddMod = 0;
    long successCode = 0;

}