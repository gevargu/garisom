#include "00MainProgram.h"
// The functions here go in order of implementation

/* Reading parameter data sheet*/
std::istream& operator>>(std::istream& str, CSVRow& data){
    data.readNextRow(str);
    return str;
    }

void MainProgram::readPARSheet (){   
    // This function reads the model configuration and the parameter datasheet
    // This is done 1 time at the beginning of the model
    // set the model control file name
    std::string mconfigFileName = "configuration_2.0.0.csv";

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
    std::string paramFileName = "parameters_2.0.0.csv";
    
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
// I placed the  instructions from readIterationSettings() function in v 1.0.24 here
void MainProgram::setConfig(){
    // // Model Configurations
    // Soil
    std::cout << "             Ground water flow: "; //i_gWaterEnable
    // turns on/off groundwater flow. Values: off; on
    if (configCells[4][1] == "on"){
        ground = true;
        std::cout << "On" << endl;
    } else {
        ground = false;
        std::cout << "Off" << endl;
    }
    std::cout << "     Soil water redistribution: "; //i_soilRedEnable
    // turns on/off soil redistribution routine. Values: off; on
    if(configCells[4][2] == "on"){
        soilred = true;
        std::cout << "On" << endl;
    } else {
        soilred = false;
        std::cout << "Off" << endl;
    }
    std::cout << "        Soil water evaporation: "; //i_soilEvapEnable
    // turns on/off soil evaporation routine. Values: off; on
    if(configCells[4][3] == "on"){
        sevap = true;
        std::cout << "On" << endl;
    } else {
        sevap = false;
        std::cout << "Off" << endl;
    }
    // Climate
    std::cout << "               Rainfall inputs: "; //i_rainEnable
    // turns on/off rain inputs. Values: off; on
    if(configCells[4][4] == "off"){
        rainEnabled = false;
        std::cout << "On" << endl;
    } else {
        rainEnabled = true; //enabled is made the "else" case becasue the "default" state is to process rain, so any input OTHER than //off// enables
        std::cout << "Off" << endl;
    }
    std::cout << "        Multiple years enabled: "; 
    // turns on/off if we are working with multiple growing seasons. Values: off; on
    if (stage_ID == STAGE_ID_NONE || stage_ID == STAGE_ID_HIST_STRESS || stage_ID == STAGE_ID_FUT_STRESS || stage_ID == STAGE_ID_FUT_STRESS_NOACCLIM){
        if (configCells[4][5] == "on") {//"i_useGSDataStress"
            useGSData = true;
            std::cout << "On" << endl;
        } else {
            useGSData = false;
            std::cout << "Off" << endl;
        } // endif
    } else if (stage_ID == STAGE_ID_HIST_OPT || stage_ID == STAGE_ID_FUT_OPT){
        if (configCells[4][14] == "on") { //"i_useGSDataOpt" different parameter name
            useGSData = true;
            std::cout << "On" << endl;
        } else {
            useGSData = false;
            std::cout << "Off" << endl;
        } // endif
    } else {// impossible unknown stage failsafe
        useGSData = false;
        std::cout << "Off" << endl;
    }

    // hydraulics
    std::cout << "               Xylem refilling: ";
    // turns on/off xylem refilling within a growing season. Values: off; on
    if(configCells[4][6] == "on"){ //"i_refilling"
        refilling = true;
        std::cout << "On" << endl;
    } else {
        refilling = false;
        std::cout << "Off" << endl;// default
    }
    std::cout << "                Pre-dawns mode: ";
    // turns on/off if model should consider measured pre-dawn water potential values. Values: off; on
    if(configCells[4][7] == "on"){//"i_predawnsMode"
        mode_predawns = true;
        std::cout << "On" << endl;
        std::cout << "MODE: Running from predawn inputs (in rain column, MPa). Soil sim disabled." << std::endl;
    } else {
        mode_predawns = false;// default
        std::cout << "Off" << endl;
        std::cout << "MODE: Running from soil simulation (calculated predawns)." << std::endl;
    }
    std::cout << "              Xylem hysteresis: "; 
    // turns on/off xylem hysteresis from previous growing season. Values: off; on 
    if(configCells[4][15] == "on"){// "i_xylemHysteresis"
        hysteresis = true;
        std::cout << "On" << endl;
    } else {
        hysteresis = false;// default
        std::cout << "Off" << endl;
    }

    // BA:GA optimization routine
    std::cout << "        Iteration Supply Curve: ";
    // Turns on the iterations in the BAGA optimization routine. Values: off; on
    if(configCells[4][8] == "on"){ // "i_iter_runSupplyCurve"
        iter_runSupplyCurve = true;
        std::cout << "On" << endl;
    } else {
        iter_runSupplyCurve = false; // default
        std::cout << "Off" << endl;
    }
    std::cout << "                Use Area Table: ";
    // Are we using BA:GA values from another csv file. Values: off; on
    if(configCells[4][9] == "on"){// "i_iter_useAreaTable"
        iter_useAreaTable = true;
        std::cout << "On" << endl;
    } else {
        iter_useAreaTable = false; // default
        std::cout << "Off" << endl;
    }
    std::cout << "          Iterate Ground Water: ";
    // Iterating ground water in for each stand. Values: off; on
    if(configCells[4][10] == "on"){// "i_iter_gwEnable"
        iter_gwEnable = true;
        std::cout << "On" << endl;
    } else {
        iter_gwEnable = false;// default
        std::cout << "Off" << endl;
    }
    // we use 1 value for all species, for now...
    iter_gwInc = std::atof(paramCells[4][62].c_str()); //how much should be increment the ground water every iteration? Smaller = higher resolution supply curve
    iter_gwStart = std::atof(paramCells[4][63].c_str()); //Some data sets may not respond to increasing ground water until a high threshold, so can start with a minimum level and increment from there
    iter_gwEnd = std::atof(paramCells[4][64].c_str());

    std::cout << "  Iterate Field Capacity (FFC): ";
    // Iterating ground water in for each stand. Values: off; on
    if(configCells[4][11] == "on"){ // "i_iter_ffcEnable"
        iter_ffcEnable = true;
        std::cout << "On" << endl;
    } else {
        iter_ffcEnable = false;// default
        std::cout << "Off" << endl;
    }
    // we use 1 value for all species, for now...
    iter_ffcInc = std::atof(paramCells[4][67].c_str());
    iter_ffcStart = std::atof(paramCells[4][65].c_str());
    iter_ffcEnd = std::atof(paramCells[4][66].c_str());
    
    std::cout << "           Iterate Field BA:GA: ";
    // Iterating BA:GA for each stand. Values: off; on
    if(configCells[4][12] == "on"){// "i_iter_bagaEnable"
        iter_bagaEnable = true;
        std::cout << "On" << endl;
    } else {
        iter_bagaEnable = false;//default
        std::cout << "Off" << endl;
    }
    // we use 1 value for all species, for now... I should move this to initialization function
    iter_bagaInc = std::atof(paramCells[4][60].c_str());
    iter_bagaStart = std::atof(paramCells[4][58].c_str());
    iter_bagaEnd = std::atof(paramCells[4][59].c_str());
    iter_bagaRef = std::atof(paramCells[4][57].c_str());
    iter_bagaCutoff = std::atof(paramCells[4][61].c_str());
    iter_bagaRef = 1.0; //156.269; // 1.0; // TODO TEMP can remove after param sheets updated
    iter_bagaEnd = 500.0; // TODO TEMP " -> for bisection method, allow extreme range

    std::cout << "       Iterate Years as Counts: ";
    // Iterating ground water in for each stand. Values: off; on
    if(configCells[4][13] == "on"){// "i_iter_"
        iter_yearsAsCount = true;
        std::cout << "On" << endl;
    } else {
        iter_yearsAsCount = false;// default
        std::cout << "Off" << endl;
    }

    // plant community version
    std::cout << "            Multi-Species Mode: ";
    // soil evaporation during non growing season (nonGS)  if nonGS_water == on. Values: off; on
    if(configCells[4][16] == "on"){ // "i_multipleSP"
        species_no = std::atoi(configCells[3][17].c_str());
        std::cout << "On" << endl;
    } else {
        species_no = 1; // default
        std::cout << "Off" << endl;
        std::cout << "MODE: Setting species number to: " << species_no << std::endl;

    }
}

/* Reading growing season data sheet*/
void MainProgram::readGSSheet(){
    memset(GSCells, 0, sizeof(GSCells));

    if (stage_ID == STAGE_ID_NONE)
    {
        GSFileName = "seasonlimits.csv";
    }

    if (useGSData == true && GSFileName != "") // only if we actually selected a file, and we're using the GS data files in the first place
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
    long startRow = 1; // skipping the header
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

/* Extract CO2 from growing season data*/
double MainProgram::getCarbonByYear(const long& yearNum, std::string GSCells[101][11], const long& maxYears) {
    long startRow = 2; // skipping the header
    double ca_year;
    for (long gsC = 0; gsC <= maxYears; gsC++){ // note this list had no header so we don't skip the first row  
        gsC = gsC;
        if (yearNum == std::atol(GSCells[gsC+startRow][1].c_str()) && std::atof(GSCells[gsC+startRow][4].c_str()) > 0) {
            //return carbonListCells[carbCount][2];
            ca_year = std::atof(GSCells[gsC+startRow][4].c_str()); // get current year atmospheric CO2 in ppm
        }
    }
    return ca_year;
}

/* Reading climate forcing and output data sheets*/
void MainProgram::readDataSheet(){
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
void MainProgram::InitialConditions(){
    /* 
    - The following loop will allow to extract parameters for each species/PFT in a single 
    location or new parameters.
    - For 1 species/PFT there is only 1 line
    TO-DO: vectorize the variables or use an array for storing values
    */
    long row_limit = species_no + 2;
    for (long sp = 0; sp <= row_limit; sp++) {
        if(sp + 1 <= 3){// + 1 helps as the loop starts in 0
            // do nothing, we skip the first 3 rows
        } else {
            std::cout << "INIT: Setting up model parameters. (year count = " << gs_yearIndex << ")" << std::endl;
            
            // site identifiers & parameters
            species = paramCells[sp+1][1];                               //"i_sp": species or PFT represented in parameter data
            name_species = species;
            region = paramCells[sp+1][2];                                //"i_region": name of the region where the simulations take place
            name_region = region;
            siteID = paramCells[sp+1][3];                                //"i_site": site/simulation ID. Usually the plot/stand identifier  
            name_site = siteID;
            lat = std::atof(paramCells[sp+1][4].c_str());                //"i_latitude": latitude in degree fraction north
            lat = pi / 180 * lat; //converted degrees to radians
            longitude = std::atof(paramCells[sp+1][5].c_str());          //"i_longitude": longitude in degree fraction west
            alt = std::atof(paramCells[sp+1][6].c_str());                //"i_elevation": elevation in m
            slope = std::atof(paramCells[sp+1][7].c_str());              //"i_slopei": slope inclination, degrees from horizontal
            slopeasp = std::atof(paramCells[sp+1][8].c_str());           //"i_slopeA": slope aspect, counterclockwise degrees from south
            pground = std::atof(paramCells[sp+1][9].c_str());            //"i_gWaterP": ground water pressure
            grounddistance = std::atof(paramCells[sp+1][10].c_str());    //"i_gWaterDist": distance to ground water source in m
            
            // atmosphere
            tau = std::atof(paramCells[sp+1][11].c_str());               //"i_atmTrans": atmospheric transmittance from weather data
            tsncorr = std::atof(paramCells[sp+1][12].c_str());           //"i_solarNoon": solar noon correction from weather data
            emiss = std::atof(paramCells[sp+1][13].c_str());             //"i_emiss": long wave emissivity
            ca = std::atof(paramCells[sp+1][14].c_str());                //"i_co2AmbPPM": ambient co2 in ppm
            ca = ca * 0.000001; //ambient co2 in moles per mole
            //the value in this cell will not be referenced again by the model

            // soil
            layers = std::atol(paramCells[sp+1][15].c_str());            //"i_layers": number of soil layers (1-5)
            if (mode_predawns){// set to 1 if predawns mode
                layers = 1;
            }
            fieldcapfrac = std::atol(paramCells[sp+1][16].c_str());      //"i_fieldCapFrac": fraction that field capacity is of saturation(minus residual)
            ffc =  std::atol(paramCells[sp+1][17].c_str())/ 100.0;       //"i_fieldCapPercInit": fraction of field capacity for starting the season
            rockfrac = std::atof(paramCells[sp+1][18].c_str());          //"i_rockFrac": fraction of soil volume as rocks
            rockfrac = 1 - rockfrac; //fraction of volume with no rocks
            soilabssol = std::atof(paramCells[sp+1][19].c_str());        //"i_soilAbsSol": absorptivity of soil surface for solar
            rhizotarg = std::atof(paramCells[sp+1][20].c_str()) / 100.0; //"i_rhizoPer": average fraction of whole plant resistance in rhizosphere(maximum soil limitation)
            texture = paramCells[sp+1][21];                              //"i_texture": texture based on USDA soil classification
            std::cout << "texture:" << texture << std::endl;

            // stand
            baperga = std::atof(paramCells[sp+1][22].c_str()) * 0.0001;  //"i_baperga": basal area per ground area converted from m2 / Ha to m2 / m2
            lai = std::atof(paramCells[sp+1][23].c_str());               //"i_leafAreaIndex": canopy lai
            xh = std::atof(paramCells[sp+1][24].c_str()) * 100;          //"i_soilXHeight": height above soil surface for understory wind and gh in cm
            height = std::atof(paramCells[sp+1][25].c_str());            //"i_height": average tree height in m
            treeToPhotoLAI = std::atof(paramCells[sp+1][26].c_str());    //"i_treeToPhotoLAI"
            laperba = std::atof(paramCells[sp+1][27].c_str());           // i_leafPerBasal: initial leaf area per basal area, m2 m - 2
            
            // tree
            leafwidth = std::atof(paramCells[sp+1][28].c_str()) * 0.72;  //"i_leafWidth": leaf width x factor = characteristic dimension(campbell and norman)
            xang = std::atof(paramCells[sp+1][29].c_str());              //"i_leafAngleParam": leaf angle parameter, CN 15.4
            aspect = std::atof(paramCells[sp+1][30].c_str());            // "i_aspect": max radius of root system per max depth
            //root_depth = std::atof(paramCells[sp+1][31].c_str());        //"i_rootDepth": maximum rooting depth
            // [TO-DO: replace by function here and rooting depth]
            beta = std::atof(paramCells[sp+1][32].c_str()); //"i_rootBeta": root beta for Y = 1 - B ^ d // amount of root biomass above the max rooting depth 

            // hydraulics
            // TO-DO replace here for p12 and p50
            leafpercent = std::atof(paramCells[sp+1][33].c_str());      //"i_leafPercRes": saturated % of tree R in leaves
            ksatp = std::atof(paramCells[sp+1][34].c_str());            //"i_kmaxTree": kmax of tree in kg hr - 1m - 2MPa - 1 per basal area
            pinc = std::atof(paramCells[sp+1][35].c_str());             // "i_pinc" MPa increment for global K(P) curves...has to be small enough for NR convergence.
            if (pinc <= 0.0){
                pinc = 0.00075;
            } // override if it's zero
            root_b = std::atof(paramCells[sp+1][36].c_str()); // = br //weibull b for each root element
            root_c = std::atof(paramCells[sp+1][37].c_str()); // = cr //weibull c for each root element
            stem_b = std::atof(paramCells[sp+1][38].c_str()); // = bs //weibull parameters for stem
            stem_c = std::atof(paramCells[sp+1][39].c_str()); // = cs
            leaf_b = std::atof(paramCells[sp+1][40].c_str()); // = bl //weibull parameters for leaf
            leaf_c = std::atof(paramCells[sp+1][41].c_str()); // = cl

            // photosynthesis
            qmax = std::atof(paramCells[sp+1][42].c_str());              //"i_qMax": quantum yield of electron transport, moles e per mols photons
            vmax25 = std::atof(paramCells[sp+1][43].c_str());            //"i_vmax25": vmax at 25C
            jmax25 = std::atof(paramCells[sp+1][44].c_str());            //"i_jmax25": jmax at 25C
            //rday25 = Cells(8, 21) //day respiration at 25C, = 0.015vmax cf.collatz et al.
            kc25 = std::atof(paramCells[sp+1][45].c_str());              //"i_kc25": m - m constant for CO2 in mole fraction at 25C
            ko25 = std::atof(paramCells[sp+1][46].c_str());              //"i_ko25": m - m constant for O2 in mole fraction at 25C
            comp25 = std::atof(paramCells[sp+1][47].c_str());            //"i_comp25": photorespiratory compensation point in mole fraction at 25C
            thetac = std::atof(paramCells[sp+1][48].c_str());            //"i_thetaC": shape factor for A - ci colimitation, 0.98
            havmax = std::atof(paramCells[sp+1][49].c_str());            //"i_havmax": these are all temp - dependency parameters from Leunig 2002
            hdvmax = std::atof(paramCells[sp+1][50].c_str());            //"i_hdvmax"
            svvmax = std::atof(paramCells[sp+1][51].c_str());            //"i_svvmax"
            lightcurv = std::atof(paramCells[sp+1][52].c_str());         //"i_lightCurv"
            lightcomp = std::atof(paramCells[sp+1][53].c_str());         //"i_lightComp"
            hajmax = std::atof(paramCells[sp+1][54].c_str());            //"i_hajmax"
            hdjmax = std::atof(paramCells[sp+1][55].c_str());            //"i_hdjmax"
            svjmax = std::atof(paramCells[sp+1][56].c_str());            //"i_svjmax"
            
            // 2. Calculating initial conditions
            /* 
            - This section calculate the variables used in the first iteration. 
            - It is within the species loop as in the future can be easily modified to calculate everything as vector of species or pfts
            - TO-DO: vectorize values for multiple species
            */
            /* Site initial conditions */
            runmean = 1; //running mean for profit maximization
            cutoff = 1.1; //cutoff for stopping dpamax search
            minwind = 0.4515; //m s-1'minimum wind threshold
            f = 70; //itmax for trapzd routine used for xylem only
            epsx = 0.0001; //fractional acceptable error for trapzd routine in xylem
            sthresh = 1.0001; //acceptable error for e integral for xylem
            
            // get atmospheric pressure
            patm = 101.325 * pow((1 - 0.0065 * alt / (288.15 + 0.0065 * alt)), 5.257); //atmospheric pressure, T = 15 C, average sealevel patm; approximation
            //Call setValueFromName("o_atmP", patm) //Cells(1, 18) = patm //atmospheric pressure in kPa

            // hydraulics
            gmax = 1000000; //maximum G, wet soil, vpd = 0, in kg m - 2 hr - 1 basal area
            gmaxl = gmax * (1.0 / laperba) * (1 / 3600.0) * 55.56 * 1000.0; //convert to gmax per leaf area in mmol m-2s-1 //'Cells(3, 13) = gmaxl
            //[HNT] calculate and output the stem and root percent resistances @ ksat -- originally was done by formula in sheet
            // added stemPercent and rootPercent globals at start
            stempercent = 1.0 / 3.0 * (100.0 - leafpercent);
            rootpercent = 2.0 / 3.0 * (100.0 - leafpercent);
            ksatl = ksatp * (100.0 / leafpercent); //leaf conductance per basal area
            lsc = ksatl * 1.0 / laperba; //lsc in kg hr-1m-2MPa-1//lsc per leaf area in kg hr - 1
            rsatp = 1.0 / ksatp; //convert to resistance
            ksatroot = 1.0 / ((rootpercent / 100.0) * rsatp); //kmax of root system; assumes zero % rhizosphere resistance in WET soil
            ksats = 1.0 / ((stempercent / 100.0) * rsatp); //kmax of stem system
            //Call setValueFromName("o_ksatLeaf", ksatl) //Cells(2, 4) = ksatl //parallel kmax for leaves
            //dbh = ((ksatp * 0.785) / Cells(6, 9)) ^ (1 / (Cells(7, 9) - 2)) //basal diameter in m from k by D allometry
            pgrav = height * 0.01; //pressure drop from gravity in MPa
            einc = ksatp / 500.0; //e increment in kg hr-1 m-2 basal area for composite curve
            kmin = ksatp / 2000.0; //"instantaneous K" cutoff for global K(P) curves for each element

            // soil parameters
            for (k = 1; k <= layers; k++){//set layer depths and % root ksat
                layerDepths[k] = 0.01 * log(1.0 - k * 0.995 / layers) / log(beta); //lower depth of each layer converted to m
                soillayersTable[rowLR+k][colLR+1] = std::to_string(k);
                soillayersTable[rowLR+k][colLR+11] = std::to_string(layerDepths[k]);
                //paramCells[rowLR + k][colLR + 11] = std::to_string(layerDepths[k]); // now calculating soil initial conditions using 2D array format
            }
            depthmax = layerDepths[layers];
            //calculate transport distances
            //first get vertical distance to biomass center of each layer
            for (k = 1; k <= layers * 2.0; k++){
                vertdistance[k] = 0.01 * log(1 - k * 0.995 / (layers * 2.0)) / log(beta); //get half depths
            }
            i = 0;
            for (k = 1; k <= layers * 2; k += 2){// To layers * 2 Step 2
                i = i + 1;
                vertdistance[i] = vertdistance[k]; //take every other vertdistance
            }
            //now get radial distances
            for (k = 1; k <= layers; k++){//To layers //get thicknesses
                if (k == 1){
                    depth[k] = layerDepths[k]; //Cells(8 + k, 11)
                } else {
                    depth[k] = layerDepths[k] - layerDepths[k - 1]; //lSheet.Cells(rowLR + k, colLR + 11) - lSheet.Cells(rowLR + k - 1, colLR + 11)  //depth is layer thickness in meters
                } // endif
            }
            vol = depth[1] * pi * pow((depthmax * aspect), 2.0); //volume of first, and hence all, layers
            //get radial widths of each layer and transport length
            for (k = 1; k <= layers; k++){//To layers
                radius[k] = pow((vol / (depth[k] * pi)), 0.5); //width in m
                length[k] = radius[k] + vertdistance[k]; //transport distance
                if (k == 1){
                    shallow = length[k];
                }
                length[k] = length[k] / shallow; //normalize to shallowest layer
            }
            unknowns = layers + 1; //number of unknowns to be solved for; also dimensions of matrix
            depth[0] = 0; //depth of surface
            
            // get Van Genuchten parameters and propagate for each layer
            // if different layers have different textures, include the function inside
            soilcalculator.get_vgparams(texture,layers,soillayersTable,rowLR,colLR);
            for (k = 1; k <= layers; k++){//To layers //read in soil
                //kkmax(k) = Cells(8 + k, 3) //saturated conductivity of soil in kg hr - 1MPa - 1 m - 1
                //a[k] = atof(paramCells[rowLR + k][colLR + 3].c_str());
                a[k] = atof(soillayersTable[rowLR + k][colLR + 3].c_str()); //van genuchten alpha
                n[k] = atof(soillayersTable[rowLR + k][colLR + 4].c_str()); //lSheet.Cells(rowLR + k, colLR + 4) //vg n
                kkmax[k] = atof(soillayersTable[rowLR + k][colLR + 5].c_str());//saturated conductivity of soil in kg hr-1MPa-1 m-1
                thetasat[k] = atof(soillayersTable[rowLR + k][colLR + 6].c_str()); //theta sat in volume/volume
                thetasat[k] = thetasat[k] * rockfrac; //reduce for actual rock-free fraction of soil
            }
            //now add toplayer (layer 0) of rootless soil 2 cm thick w. same properties as layer 1
            a[0] = a[1];
            n[0] = n[1];
            kkmax[0] = kkmax[1];
            thetasat[0] = thetasat[1];
            depth[0] = 0.02; //sets top layer to 2 cm
            //now solve for kmax rhizosphere that gives the desired ave % rhizosphere resistance
            z = 1; //use layer 1 as stand in for whole root system
            ksatr[1] = ksatroot; //set to whole root system
            x = 0.5; //start by finding kmaxrh at 0.5 MPa...a deliberate under-shoot
            rootr = 1.0 / get_weibullfitroot(x);
            rstem = 1.0 / get_weibullfit(x,stem_b,stem_c,ksats);
            rleaf = 1.0 / get_weibullfit(x,leaf_b,leaf_c,ksatl);
            rplant = rootr + rstem + rleaf; //rplant here is just the xylem part
            rhizor = rplant * (rhizotarg / (1.0 - rhizotarg)); //solve for what rhizor has to be at the target
            
            // rhizosphere
            vp = 1.0 / (pow((a[z] * x), n[z]) + 1); //van genuchten terms // vp = 1 / ((a(z) * x) ^ n(z) + 1) 
            vgterm = pow(vp, ((n[z] - 1) / (2.0 * n[z]))) * pow((pow((1 - vp), ((n[z] - 1) / n[z])) - 1), 2.0); //van genuchten terms // vgterm = vp ^ ((n[z] - 1) / (2 * n[z])) * ((1 - vp) ^ ((n[z] - 1) / n[z]) - 1) ^ 2
            kmaxrh[1] = (1.0 / rhizor) / vgterm; //solve for kmaxrh[1]
            kinc = kmaxrh[1] * 0.1;
            do{//loop through rhizosphere kmax
                kmaxrh[1] = kmaxrh[1] + kinc; //increase from deliberate undershoot
                x = 0; //
                sum = 0;
                do{//loop through pressures
                    x = x + 0.1;
                    rootr = 1.0 / get_weibullfitroot(x);
                    rstem = 1.0 / get_weibullfit(x,stem_b,stem_c,ksats);
                    rleaf = 1.0 / get_weibullfit(x,leaf_b,leaf_c,ksatl);
                    rhizor = 1.0 / get_vg(x);
                    rplant = rootr + rstem + rleaf + rhizor;
                    rrhizofrac = rhizor / rplant; //fraction of resistance in rhizosphere
                    sum = sum + rrhizofrac; //add up fractions
                } while (!(1.0 / rplant < kmin)); //Loop Until 1 / rplant < kmin //average over full range
                sum = sum / (x / 0.1); //average fraction
            } while (!(sum < rhizotarg)); // Until sum < rhizotarg //loop until desired soil limitation is reached

            kmaxrh[1] = kmaxrh[1] / layers; //divide whole root rhizokmax into equal portions for each layer
            //end of soil limitation adjustment, now set soil layer parameters based on aroot of entire root system
            for (z = 1; z <= layers; z++){//z = 1 To layers
                kmaxrh[z] = kmaxrh[1]; //soil to root MAXIMUM conductance in kg hr-1 MPa-1//re - set for individual layers
                //lSheet.Cells(rowLR + z, colLR + 7) = kmaxrh(z) //note: this is kMAX...at P=0; not at saturated PD
                soillayersTable[rowLR + z][colLR + 7] = std::to_string(kmaxrh[z]); //note: this is kMAX...at P=0; not at saturated PD
            }
            t = 0;
            //loop to find ksatroot for each layer
            coef = 0.0;
            do{
                coef = coef + 0.01;
                sum = 0.0;
                for (k = 1; k <= layers; k++){//k = 1 To layers //soil layers from top to bottom
                    ksatr[k] = coef / length[k]; //assumes ksatr proportional to biomass/length
                    sum = sum + ksatr[k];
                }
            } while (!(sum > ksatroot));//Loop Until sum > ksatroot //loop until each layer adds to total

            for (k = 1; k <= layers; k++){//k = 1 To layers
                //lSheet.Cells(rowLR + k, colLR + 2) = 100 * (ksatr(k) / ksatroot);
                //lSheet.Cells(rowLR + k, colLR + 9) = ksatr(k);
                soillayersTable[rowLR + k][colLR + 2] = std::to_string(100.0 * (ksatr[k] / ksatroot));
                soillayersTable[rowLR + k][colLR + 9] = std::to_string(ksatr[k]);
            }
            
            rough = 0.01; //soil Zm, eqn 14.9, using table 5.1 for smooth surface, cm
            zdispl = 6.5 * rough; //soil d, eqn 14.9, using d = 6.5 Zm, eq 5.2,5.3
            zh = 0.2 * rough; //roughness for temperature
            //rainsim = getValueFromName("i_rainEnable"); //Cells(11, 27) //y/n...turns on simulated rain
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
short MainProgram::cleanModelVars(){
    // Model arrays
    long iii = 0;
    long qqq = 0;
    // 6 row array or 100K
    for (iii = 0; iii < 6; iii++){
        //6 arrays 1-d
        n[iii] = 0;
        a[iii] = 0;
        //ksatrh[iii] = 0; // not used in the model
        kmaxrh[iii] = 0;
        ksatr[iii] = 0;
        kmaxr[iii] = 0;
        depth[iii] = 0;
        pcritr[iii] = 0;
        prh[iii] = 0;
        kminroot[iii] = 0;
        length[iii] = 0;
        radius[iii] = 0;
        layer[iii] = 0;
        tlayer[iii] = 0;
        layerfailure[iii] = "ok";
        tlayerfailure[iii] = "ok";
        soilf[iii] = 0;
        kkmax[iii] = 0;
        thetasat[iii] = 0;
        swclimit[iii] = 0;
        
        if (!(useGSData && gs_yearIndex > 0)){
            thetafc[iii] = 0;
            thetafracres[iii] = 0;
            thetafracfc[iii] = 0;
            water[iii] = 0;
            fc[iii] = 0;
        }
        pcritrh[iii] = 0;
        pd[iii] = 0;
        soilredist[iii] = 0;

        for (qqq = 0; qqq < 1000001; qqq++){
            // 1-D 100k and 1mil arrays
            if (iii == 0){
               // 100k arrays
                if (qqq < 100001){
                    psynshmd[qqq] = 0;
                    psynmd[qqq] = 0;
                    psynsh[qqq] = 0;
                    amaxfracsh[qqq] = 0;
                    leaftempsh[qqq] = 0;
                    lavpdsh[qqq] = 0;
                    rdaysh[qqq] = 0;
                    gcanwsh[qqq] = 0;
                    gcancsh[qqq] = 0;
                    cinsh[qqq] = 0;
                    rday[qqq] = 0;
                    leaftemp[qqq] = 0;
                    lavpd[qqq] = 0;
                    amaxfrac[qqq] = 0;
                    dpa[qqq] = 0;
                    psyn[qqq] = 0;
                    gcanw[qqq] = 0;
                    gcanc[qqq] = 0;
                    dedp[qqq] = 0;
                    pleaf[qqq] = 0;
                    kplant[qqq] = 0;
                    eplant[qqq] = 0;
                    kstem[qqq] = 0;
                    pstem[qqq] = 0;
                    proot[qqq] = 0;
                    kleaf[qqq] = 0;
                    tel[qqq] = 0;
                    tes[qqq] = 0;
                    el[qqq] = 0;
                    es[qqq] = 0;
                }// 1 mil arrays
                //kloss[qqq] = 0; // not really used in the model
                klossv[qqq] = 0;
                pleafv[qqq] = 0;
                eplantl[qqq] = 0;
                dedpf[qqq] = 0;
                cin[qqq] = 0;
            }
            // 6,100k 2d arrays 
            if (qqq < 100001){//don't mix up the 100k and 1mil arrays
                er[iii][qqq] = 0;
                kr[iii][qqq] = 0;
                erh[iii][qqq] = 0;
                krh[iii][qqq] = 0;
                tkr[iii][qqq] = 0;
                ter[iii][qqq] = 0;
                elayer[iii][qqq] = 0;
                prhizo[iii][qqq] = 0;
                kroot[iii][qqq] = 0;
            }// 6,1 mil arrays 2d
        }
    }

    // 7 row arrays
    for (iii = 0; iii < 7; iii++){
        indx[iii] = 0;
        func[iii] = 0;
        vv[iii] = 0;
        dfrhdprh[iii] = 0;
        dfrdprh[iii] = 0;
        dfrhdpr[iii] = 0;
    }

    // lazy memsets for anything that doesn't fit -- usually good enough
    std::memset(jmatrix, 0, sizeof(jmatrix));
    memset(vertdistance, 0, sizeof(vertdistance));
    //memset(inp, 0, sizeof(inp)); // not used in the model
    //jmatrix[7][7] = 0;
    //vertdistance[11] = 0;
    //inp[101];
    
    // Model variables
    e = 0;
    x = 0;
    ksh = 0;
    root_b = 0;
    root_c = 0;
    stem_b = 0;
    stem_c = 0;
    leaf_b = 0;
    leaf_c = 0;
    ksats = 0;
    ksatl = 0;
    pcrits = 0;
    pcritl = 0;
    ps = 0;
    pinc = 0;
    einc = 0;
    ksatp = 0;
    rsatp = 0;
    ksatroot = 0;
    vp = 0;
    s = 0;
    del = 0;
    eps = 0;
    epsx = 0;
    sthresh = 0;
    sum = 0;
    //sums = 0; // not used by the model
    olds = 0;
    kmin = 0;
    elow = 0;
    ehigh = 0;
    plow = 0;
    estart = 0;
    efinish = 0;
    flow = 0;
    dfrdpr = 0;
    frt = 0;
    klow = 0;
    khigh = 0;
    kupper = 0;
    klower = 0;
    pr = 0;
    threshold = 0;
    //kfactor; // not used by the model
    aamax = 0;
    dum = 0;
    kminleaf = 0;
    kminstem = 0;
    md = 0;
    ecritsystem = 0;
    //pcritsystem;
    phigh = 0;
    failspot = "";
    //setting = ""; // not used by the model
    //refilling = ""; // set as boolean in v 2.0
    //ground = ""; // set as boolean in v 2.0
    //soilred = ""; // set as boolean in v 2.0
    //dp = 0; // not used by the model
    kminplant = 0;
    //toplayer; // not used by the model
    p1 = 0;
    p2 = 0;
    pl = 0;
    plold = 0;
    predawn = 0;
    dedpl = 0;
    //prtarget = 0; // not used by the model
    //kmaxs = 0; // not used by the model
    //pstarget = 0; // not used by the model
    //kmaxl = 0; // not used by the model
    //pltarget = 0; // not used by the model
    //rhizsat = 0; // not used by the model
    //rootsat = 0; // not used by the model
    //sum2 = 0; // not used by the model
    beta = 0;
    //rhizok = 0; // not used by the model
    patm = 0;
    vpd = 0;
    //dedplzero;
    gcmd = 0;
    dpmax = 0;
    gmax = 0;
    //vpdsat = 0; // not used by the model
    //pstop = 0; // not used by the model
    halt = 0;
    failure = 0;
    k = 0;
    check = 0;
    t = 0;
    j = 0;
    it = 0;
    f = 0;
    //dd = 0;
    tmax = 0;
    tnm = 0;
    test = 0;
    layers = 0;
    d = 0;
    z = 0;
    p = 0;
    //pmax = 0; // not used by the model
    unknowns = 0;
    i = 0;
    imax = 0;
    ii = 0;
    ll = 0;
    total = 0;
    kplantold = 0;
    //frac = 0; // not used by the model
    //eincdef = 0; // not used by the model
    aspect = 0;
    vol = 0;
    depthmax = 0;
    shallow = 0;
    coef = 0;
    rhizotarg = 0;
    rplant = 0;
    rstem = 0;
    rootr = 0;
    rleaf = 0;
    rhizor = 0;
    rrhizofrac = 0;
    kinc = 0;
    //kplantmax = 0; // not used by the model
    vgterm = 0;
    ci = 0;
    ca = 0;
    marker = 0;
    var = 0;
    transpiration = 0;
    //dedpmax = 0; // not used by the model
    lavpdmd = 0;
    psynact = 0;
    psynmax = 0;
    //dpamax;
    grad = 0;
    gha = 0;
    numerator = 0;
    denominator = 0;
    lambda = 0;
    emiss = 0;
    airtemp = 0;
    rabs = 0;
    laperba = 0;
    par = 0;
    qmax = 0;
    vmax25 = 0;
    kc25 = 0;
    ko25 = 0;
    comp25 = 0;
    theta = 0;
    wind = 0;
    leafwidth = 0;
    comp = 0;
    vmax = 0;
    kc = 0;
    ko = 0;
    je = 0;
    jc = 0;
    //lavpdc1 = 0; // not used by the model
    //lavpdsum = 0; // not used by the model
    //lavpdh = 0; // not used by the model
    gmaxl = 0;
    cinc = 0;
    jmax = 0;
    jmax25 = 0;
    havmax = 0;
    hdvmax = 0;
    svvmax = 0;
    hajmax = 0;
    hdjmax = 0;
    svjmax = 0;
    jact = 0;
    maxvpd = 0;
    maxkloss = 0;
    chalk = 0;
    skip = 0;
    //psynstop = 0; // not used by the model
    totalv = 0;
    lsc = 0;
    timestep = 0;
    leafpercent = 0;
    //dbh = 0; // not used by the model
    height = 0;
    pgrav = 0;
    ffc = 0;
    pground = 0;
    grounddistance = 0;
    baperga = 0;
    layerflow = 0;
    // groundwater = 0; // not used by the model
    deficit = 0;
    rain = 0;
    //sumrain = 0; // not used by the model
    store = 0;
    pend = 0;
    groundflow = 0;
    rday25 = 0;
    if (!(useGSData && gs_yearIndex > 0)){
        waterold = 0;
        waternew = 0;
        waterchange = 0;
    }
    drainage = 0;
    gwflow = 0;
    thetac = 0;
    lat = 0;
    longitude = 0;
    slope = 0;
    slopeasp = 0;
    lai = 0;
    xang = 0;
    fet = 0;
    et = 0;
    sm = 0;
    lc = 0;
    tsn = 0;
    sindec = 0;
    dec = 0;
    cosdec = 0;
    tim = 0;
    tod = 0;
    coszen = 0;
    zen = 0;
    cosaz = 0;
    az = 0;
    m = 0;
    sp = 0;
    sb = 0;
    sd = 0;
    st = 0;
    cloud = 0;
    obssolar = 0;
    fcd = 0;
    kbe = 0;
    kbezero = 0;
    mleafang = 0;
    rad = 0;
    k1 = 0;
    t1 = 0;
    told = 0;
    t2 = 0;
    kd = 0;
    qd = 0;
    qds = 0;
    qdt = 0;
    qb = 0;
    qbt = 0;
    qsc = 0;
    qsh = 0;
    qsl = 0;
    laisl = 0;
    laish = 0;
    parsh = 0;
    parsl = 0;
    parbottom = 0;
    nirsh = 0;
    nirsl = 0;
    sshade = 0;
    ssun = 0;
    sbottom = 0;
    ssunb = 0;
    ssund = 0;
    sref = 0;
    ppfd = 0;
    ea = 0;
    eac = 0;
    la = 0;
    lg = 0;
    jd = 0;
    o = 0;
    mdsh = 0;
    cincsh = 0;
    lavpdmdsh = 0;
    gcmdsh = 0;
    psynactsh = 0;
    transpirationsh = 0;
    haltsh = 0;
    psynmaxsh = 0;
    transpirationtree = 0;
    //anetsh = 0; // not used by the model
    //anettree = 0; // not used by the model
    atree = 0;
    //anet = 0; // not used by the model
    rabssoil = 0;
    alt = 0;
    tsncorr = 0;
    rha = 0;
    rhs = 0;
    soilevap = 0;
    soilabssol = 0;
    rough = 0;
    zdispl = 0;
    xh = 0;
    zh = 0;
    us = 0;
    mdensair = 0;
    soilep = 0;
    emission = 0;
    soiltemp = 0;
    night = "";
    //sevap = ""; // set as boolean in v 2.0
    //pet = ""; // not used by the model
    //rainsim = ""; // not used by the model
    //nightBool = 0; // not used by the model
    tau = 0;
    minwind = 0;
    runmean = 0;
    weird = 0;
    //sign = 0; // not used by the model
    //xx = 0; // not used by the model
    //runs = 0; // not used by the model
    ticks = 0;
    rmean = 0;
    cutoff = 0;
    dedplmin = 0;
    amaxmax = 0;
    lightcurv = 0;
    emd = 0;
    leaftmd = 0;
    leaftshmd = 0;
    lavpdshmd = 0;
    gcanwmd = 0;
    gcancmd = 0;
    rdaymd = 0;
    psynmaxmd = 0;
    cinmd = 0;
    gcanwshmd = 0;
    gcancshmd = 0;
    rdayshmd = 0;
    cinshmd = 0;
    psynmaxshmd = 0;
    prinitial = 0;
    dpasun = 0;
    lightcomp = 0;
    fieldcapfrac = 0;
    rockfrac = 0;
    initialthreshold = 0;
    //kref = 0; // not used by the model
    //lscref = 0; // not used by the model
    //pdref = 0; // not used by the model
    //leafpercentref = 0; // not used by the model
    kpday1 = 0;
    kxday1 = 0;
    runoff = 0;
    //kmaxits = 0; // not used by the model
    //bst = 0; // not used by the model
    //startbst = 0; // not used by the model
    //scenario = 0; // not used by the model
    //kmaxset = 0; // not used by the model
    //leafpercentset = 0; // not used by the model
    return -1; // true in VBA langauge
}

/* Model iterator */
short MainProgram::setIterationParams(long& iter){
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
    if (useGSData) {
        if (gs_yearIndex >= 0 && gs_yearIndex < 100) {
            if (gs_ar_years[gs_yearIndex] > 0) {
                if (jd >= gs_ar_starts[gs_yearIndex] && jd <= gs_ar_ends[gs_yearIndex]) {
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
    colD = 1 - 1;
    rowD = 1 + 1;
    rowLR = 60;
    colLR = 1 - 1;

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

    std::cout << "Read site and species names (region | site | species): " << name_region << " | " << name_site << " | " << name_species << " | " << std::endl;
    return true;
}

/* Read site area values
TO-DO: incorporate this function into the model or delete
*/
void MainProgram::readSiteAreaValues(){
}

// Soils module
// TO-DO: finish pasting these functions to a separate file as its own class

/*Soil hydraulic conductivity*/
/* Van Genuchten function
gets soil k
*/
double MainProgram::get_vg(double &x) {
    //vp = 1 / ((a(z) * x) ^ n(z) + 1)
    //vg = kmaxrh(z) * vp ^ ((n(z) - 1) / (2 * n(z))) * ((1 - vp) ^ ((n(z) - 1) / n(z)) - 1) ^ 2
    vp = 1 / (pow((a[z] * x), n[z]) + 1);
    return kmaxrh[z] * pow(vp, ((n[z] - 1) / (2 * n[z]))) * pow((pow((1 - vp), ((n[z] - 1) / n[z])) - 1), 2);
}

/* Van Genuchten function integration
integrates VG
*/
void MainProgram::trapzdvg(double &p1, double &p2, double &s, long &t){
    if (t == 1){
        s = 0.5 * (p2 - p1) * (get_vg(p1) + get_vg(p2));
        it = 1;
    } else {
        tnm = it;
        del = (p2 - p1) / tnm;
        x = p1 + 0.5 * del;
        sum = 0;
        for (j = 1; j <= it; j++) {
            sum = sum + get_vg(x);
            x = x + del;
        }
        s = 0.5 * (s + (p2 - p1) * sum / tnm);
        it = 2 * it;
    }
}

/* VG integration accuracy
evaluates accuracy of van genuchten integration
*/
void MainProgram::qtrapvg(double &p1, double &p2, double &s) {
    eps = 0.001; //fractional refinement threshold for integral
    tmax = 70; //limit of calls to trapzdvg
    olds = -1; //starting point unlikely to satisfy if statement below
    for (t = 1; t <= tmax; t++){
        trapzdvg(p1, p2, s, t);
        if (std::abs(s - olds) < eps * std::abs(olds)){
            return;
        }
        olds = s;
    }
}

// Soil water status
/* Soil pressure
gives soil Y in MPa from soil theta/thetasat=x
*/
double MainProgram::rvg(double &x) {
    //rvg = (x ^ (1 / (1 - 1 / n(z))) + 1) ^ (1 / n(z)) / (x ^ (1 / (n(z) - 1)) * a(z))
    double aa = pow((pow(x, (1 / (1 - 1 / n[z]))) + 1), (1 / n[z]));
    double bb = (pow(x, (1 / (n[z] - 1))) * a[z]);
    return aa / bb;
}

/* Soil water content
gives soil water content, theta/thetasat from soil water potential in MPa=x
*/
double MainProgram::swc(double &x) {
    //swc = (1 / (1 + (a(z) * x) ^ n(z))) ^ (1 - 1 / n(z))
    return pow((1 / (1 + pow((a[z] * x), n[z]))), (1 - 1 / n[z]));
}

/* Soil wetness
gets predawns accounting for rain, groundwater flow, transpiration, and redistribution via soil and roots
*/
void MainProgram::get_soilwetness() {
    double tempDouble = 0.0;
    drainage = 0;
    runoff = 0;
    if (dd == 1 || isNewYear) { //if// //'every layer starts at initial % of field capacity
        if (!(useGSData && gs_yearIndex > 0)){
            waterold = 0; //'total root zone water content (mmol m-2 ground area)
            for (z = 0; z <= layers; z++){//z = 0 To layers //'
                x = 10; //'MPa water potential for getting residual thetafrac
                thetafracres[z] = swc(x); //'residual thetafrac
                thetafracfc[z] = (1 - thetafracres[z]) * fieldcapfrac + thetafracres[z]; //'thetafrac at field capacity
                thetafc[z] = thetasat[z] * thetafracfc[z]; //'water content at field capacity
            } //for//z //'
            for (z = 0; z <= layers; z++){//z = 0 To layers
                water[z] = thetafc[z] * depth[z]; //'field capacity estimated as 1/2 saturated capacity, water content of layer in m3 water per m2 ground area
                fc[z] = water[z]; //'records field capacity in m3 water volume per m2 ground area.
                water[z] = ffc * water[z]; //'start off with initial fraction of field capacity
                //'Cells(16 + dd, 42 + z) = water[z]
                waterold = waterold + water[z]; //'in m3/m2 ground area
            } //for//z
            dataCells[rowD + dd][colD + dColF_End_watercontent] = waterold * 1000; //'root zone water content in mm m-2 ground area
            //dSheet.Cells(rowD + dd, colD + dColF_End_watercontent) = waterold * 1000; //'root zone water content in mm m-2 ground area
            //'waterold = waterold * 55555556# //'convert m3 water per ground area to mmol water per ground area
            // [HNT] starting water now counts as an input
            gs_ar_input[gs_yearIndex] = gs_ar_input[gs_yearIndex] + waterold * 1000;
        }

        if (gs_yearIndex == 0){// if it's the first year, store this as the off-season starting water because we don't have a real value
            gs_ar_waterInitial_OFF[gs_yearIndex] = waterold * 1000;
        }
        // store the initial water content to check how much we consume at the end
        gs_ar_waterInitial[gs_yearIndex] = waterold * 1000;
        // [/HNT]
    } //End if// //'dd=1 if
    //'if pet = "y" Or pet = "n" { //if// //'do the original routine...doesn//'t run for PET scenario
    if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'get flows that happened during previous timestep
        for (z = 0; z < layers; z++){//z = 0 To layers - 1 //'transpiration, root and soil redistribution
            if (night == "n") { //if// //'it//'s day, must adjust elayer for sun vs. shade weighting
                layerflow = elayer[z][halt] * laisl / lai + elayer[z][haltsh] * laish / lai; //'weighted flow; NOTE: elayer = 0 for layer 0
            } else {
                layerflow = elayer[z][halt]; //'no adjustment necessary at night; NOTE: elayer = 0 for layer 0
            } //End if// //'night if
            layerflow = layerflow * baperga * 1.0 / 998.2 * timestep; //'rootflow into (= negative rootflow) or out (positive flow) of layer in m3/m2 ground area
            layerflow = layerflow + soilredist[z] * 1.0 / 998.2 * timestep; //'redistribution between layers (negative is inflow, positive is outflow). NOTE: soilredist(0) includes soil evaporation for layer 0
            water[z] = water[z] - layerflow; //'subtracts rootflow from layer on per ground area basis
        } //for//z
        //'now do the bottom layer and potential groundwater input
        if (night == "n") { //if// //'it//'s day, must adjust layerflow for sun vs. shade weighting
            layerflow = elayer[layers][halt] * laisl / lai + elayer[layers][haltsh] * laish / lai; //'weighted flow
        } else {
            layerflow = elayer[layers][halt]; //'no adjustment necessary at night
        } //End if// //'night if
        layerflow = layerflow * baperga * 1 / 998.2 * timestep; //'rootflow into (= negative rootflow) or out (positive flow) of layer in m3/m2 ground area
        layerflow = layerflow + soilredist[layers] * 1 / 998.2 * timestep; //'redistribution between layers (negative is inflow, positive is outflow)
        //'water(layers) = water(layers) - layerflow //'subtracts rootflow from layer on per ground area basis
        if (layerflow < 0) { //if// //'water is added
            for (z = layers; z >= 0; z--){//z = layers To 0 Step -1 //'start at bottom, go up
                deficit = thetasat[z] * depth[z] - water[z]; //'m of water required to wet up layer to SATURATION
                if (-1 * layerflow - deficit >= 0) { //if// //'there//'s enough to wet the layer...remember, negative flow is flow into the layer
                    water[z] = thetasat[z] * depth[z]; //'m water at saturation in layer
                    layerflow = layerflow + deficit; //'reduce what//'s left over for next layers
                } else { //'just soak up all the groundwater
                    water[z] = water[z] - layerflow; //'add to bottom layer
                    layerflow = 0; //' all gone
                } //End if// //'wetting if
            } //for//z
            runoff = runoff - layerflow; //'add what//'s left to runoff...
        } else { //'groundwater is positive...bottom layer is losing water
            water[layers] = water[layers] - layerflow; //'subtract from bottom layer
        } //End if// //'layerflow if

        //'now reset any exhausted layers to extraction limit
        if (water[0] <= 0){//'set lower limit to surface water content
            water[0] = 0.00001;
        }
        for (z = 1; z <= layers; z++){//z = 1 To layers
            if (water[z] < swclimit[z]){//'water at limit
                water[z] = swclimit[z];
            }
        } //for//z

        bool rainOverride = false;
        if (iter_Counter == 0 && tod == 23 && (stage_ID == STAGE_ID_HIST_OPT || stage_ID == STAGE_ID_FUT_OPT)){
            rainOverride = true;
        }
        //predawns mode
        if (mode_predawns){
            rainOverride = true;
        } // end predawns mode

        //'now check for rain during PREVIOUS TIME STEP
        if (dSheet.Cells(rowD + dd - 1, colD + dColRain) >= 0.0 || rainOverride == true) { //if// //'there//'s been rain!
            rain = dSheet.Cells(rowD + dd - 1, colD + dColRain) * 0.001; //'convert mm of rain to m depth...equivalent to m3/m2 water volume per ground area
            //'if raining(xx) = "n" { //if// rain = 0 //'rain turned off

            if (rainOverride){// just force everything to FC on every timestep
                for (z = 0; z <= layers; z++){//z = 0 To layers //'add in rain
                    water[z] = fc[z];
                }
            } // do this before the rain addition, so that we see all rain as runoff

            if (rainEnabled == false || mode_predawns) {
               rain = 0;
            } // predawns mode also turns off rain
            // [HNT] temp
            //if (rain > 0.1) // more than 1/10th meter per hour is a data anomoly, however this is a poor assumption and these should be checked in the weather curator
            //   rain = 0.0;
            // [/HNT]
            //'sumrain = rain * 1000 + sumrain //'total rain input in mm m-2
            for (z = 0; z <= layers; z++)//z = 0 To layers //'add in rain
            {
               if (rain <= 0)
                  break; //'rain//'s used up
               deficit = fc[z] - water[z]; //'m3 of water per m2 to wet up layer to fc
                                           //if (deficit >= -0.001) { //if// //'layer is at or below field capacity
               if (rain - deficit >= 0) { //if// //'there//'s enough to wet the layer
                  water[z] = fc[z];
                  rain = rain - deficit;
                  //'sumrain = sumrain + deficit //'absorbed rain
                  drainage = rain * 1000.0; //'any left over will drain out the bottom unless layer is rising above field capacity
               }
               else {
                  water[z] = water[z] + rain; //'it all goes into the first layer
                                              //'sumrain = sumrain + rain
                  rain = 0; //'rain used up
                  drainage = 0;
               } //End if// //'wetting up to field capacity "if"
                 //}
            } //for//z

              // If there's drainage, and the ground water is on, that should be used to fill up to saturation
            if (rain > 0.0 && ground == "y") // the remaining "drainage" is also still stored in the rain variable
            {
               // this is kind of inefficient, but the rain routine actually drained all the layers to FC even if GW was on and we should have been filling to sat
               // now we start at the bottom and fill the layers to saturation using the drainage

               for (j = layers; j >= 0; j--) //j = z - 1 To 0 Step -1 //'go back up to fill profile to saturation
               {
                  if (rain <= 0) //if// Exit for
                     break;
                  deficit = thetasat[j] * depth[j] - water[j];
                  if (deficit >= 0) { //if// //'got capacity
                     if (rain - deficit >= 0) { //if// //'enough rain to saturate the layer
                        water[j] = thetasat[j] * depth[j]; //'saturate the layer
                        rain = rain - deficit; //'reduce rain
                     }
                     else { //'rain absorbed by layer
                        water[j] = water[j] + rain;
                        rain = 0; //'use up rain
                     } //End if// //'deficit=>0 "if"
                  }
                  else { //'deficit<0...layer//'s saturated
                     rain = rain - deficit; //'increase rain by super-saturated amount (deficit is negative)
                     water[j] = thetasat[j] * depth[j]; //'reset to saturation
                  } //End if// //'deficit <>0 if
               } //for//j
               runoff = runoff + rain; //'whatever is left over will run off
               drainage = 0; //'no drainage if any layer is rising above field capacity
               rain = 0; //'reset rain to zero
            }
            //'sumdrain = sumdrain + drainage //'total drainage
         } //End if// //'rain if

           //'now check for exhausted layers
         if (water[0] <= 0)
            water[0] = 0.00001; //'set lower limit to surface water content
         for (z = 1; z <= layers; z++)//z = 1 To layers
         {
            if (water[z] < swclimit[z])
               layer[z] = 1; //'water exhausted
         } //for//z
           //'now get water content change over PREVIOUS time step
         if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'now get updated new water content
            waternew = 0;
            for (z = 0; z <= layers; z++)//z = 0 To layers //'check for exhausted layers
            {
               waternew = water[z] + waternew;
               //'Cells(16 + dd, 42 + z) = water[z]
            } //for//z



              //'waternew = waternew * 55555556# //'new water content at beginning of timestep
            waterchange = waternew - waterold; //'total water change in mm m-2 ground area
                                               //'waterchange = waterchange * 1 / baperga * 1 / laperba //'total water in mmol per m2 leaf area
            dSheet.Cells(rowD + dd, colD + dColF_End_watercontent) = waternew * 1000; //'root zone water content in mm

                                                                                      // always store the water content as the "final" -- not worth testing if it's really the last day of the year
            gs_ar_waterFinal[gs_yearIndex] = waternew * 1000;

            dSheet.Cells(rowD + dd, colD + dColF_End_waterchange) = waterchange * 1000; //'change in water content over PREVIOUS timestep
                                                                                        //'if raining(xx) = "y" { //if// Cells(16 + dd, 59) = Cells(16 + dd - 1, 4) //'rain input per previous timestep
            if (rainEnabled == true)// && dSheet.Cells(rowD + dd - 1, colD + dColRain) < 100.0)
               dSheet.Cells(rowD + dd, colD + dColF_End_rain) = dSheet.Cells(rowD + dd - 1, colD + dColRain); //'rain input per previous timestep
            dSheet.Cells(rowD + dd, colD + dColF_End_gwater) = gwflow; //'groundwater input in mm per timestep
                                                                       //'Cells(16 + dd, 61) = transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 //'transpiration per ground area in mm m-2 per timestop
            dSheet.Cells(rowD + dd, colD + dColF_End_drainage) = drainage; //'total drainage in mm per timestep
            dSheet.Cells(rowD + dd, colD + dColF_End_input) = dSheet.Cells(rowD + dd, colD + dColF_End_rain) + dSheet.Cells(rowD + dd, colD + dColF_End_gwater); //'total input per timestep in mm
                                                                                                                                                                 //'Cells(16 + dd, 64) = water(0) * 1000 //'water in top layer in mm
            gs_ar_input[gs_yearIndex] = gs_ar_input[gs_yearIndex] + dSheet.Cells(rowD + dd, colD + dColF_End_input);

            dSheet.Cells(rowD + dd, colD + dColF_End_runoff) = runoff * 1000; //'excess root zone water per timestep in mm
            waterold = waternew; //'reset to beginning of current timestep
         } //End if// //'dd>1 if
      } //End if// //'dd>1 if
        //'} //End if////'pet if
      if (dd > 1 && !isNewYear) { //if//
         tempDouble = transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000;
         dSheet.Cells(rowD + dd, o + colD + dColF_End_E) = tempDouble;//transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000; //'transpiration per ground area in mm m-2 per timestop
         if (gs_inGrowSeason) // only record growing season E (should not be any non-GS E, but just for safety)
            gs_ar_E[gs_yearIndex] = gs_ar_E[gs_yearIndex] + tempDouble;

         dSheet.Cells(rowD + dd, o + colD + dColF_End_soilEvap) = soilevap * 1 / 998.2 * timestep * 1000; //'evaporative water loss in mm per timestep

         tempDouble = transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 + soilevap * 1 / 998.2 * timestep * 1000;
         dSheet.Cells(rowD + dd, o + colD + dColF_End_ET) = tempDouble;//transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 + soilevap * 1 / 998.2 * timestep * 1000;
         gs_ar_ET[gs_yearIndex] = gs_ar_ET[gs_yearIndex] + tempDouble;

         tempDouble = atree * timestep * 3600 * 0.001;
         dSheet.Cells(rowD + dd, o + colD + dColF_End_ANet) = tempDouble;//atree * timestep * 3600 * 0.001; //'Anet per timestep in mmoles per leaf area
         if (!std::isnan(tempDouble) && gs_inGrowSeason) // only record growing season A
         {
            gs_ar_Anet[gs_yearIndex] = gs_ar_Anet[gs_yearIndex] + tempDouble;
            //anytime we record A, also record the ci-related outputs
            // TODO CRIT units???
            // Anet in mmoles per leaf area as in final output columns? (calc above)
            if (night == "n") // daytime only
            {
               gs_ar_cica[gs_yearIndex] += cinc / ca; // these are in mols/mol, but it's a ratio so not important
               gs_ar_cica_N[gs_yearIndex]++;

               gs_ar_Aci[gs_yearIndex] += tempDouble * cinc; // for Aci what units? pa would be (cinc * patm * 1000)
               gs_ar_AnetDay[gs_yearIndex] += tempDouble; // keep a daytime-only Anet tally
                                                          //this is for A-weighted Ci
                                                          //gs_ar_Acica[gs_yearIndex] += tempDouble * cinc / ca; // for this one, another ratio so ignore units
            }
         }
      } //End if// //'dd>1 if

      if (tod == 16 && !gs_doneFirstDay && gs_inGrowSeason && dSheet.Cells(rowD + dd - 3, colD + dColF_CP_kplant) > 0.000000001) { //if// //'get midday k//'s for day 1
                                                                                                                                   // VPD zero case -- if the stomata did not open on the first day of the GS, kplant won't have been set and will be zero... in which case, try again tomorrow
         gs_doneFirstDay = true;

         sum = 0;
         for (z = 1; z <= 3; z++)//z = 1 To 3
         {
            sum = sum + dSheet.Cells(rowD + dd - z, colD + dColF_CP_kplant);
         } //for//z
         kpday1 = sum / 3.0; //'average midday kplant on day one
         sum = 0;
         for (z = 1; z <= 3; z++)//z = 1 To 3
         {
            sum = sum + dSheet.Cells(rowD + dd - z, colD + dColF_CP_kxylem);
         } //for//z
         kxday1 = sum / 3.0; //'average midday kxylem on day one

                             // [HNT] first day of the new growing season! Record the starting water
         gs_ar_waterInitial_GS[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent); // no matter what happened above, this has been set by this point
                                                                                                       // and record this as the FINAL water for THIS YEAR's off season -- remember that this year's off season is the PRECEDING winter
         gs_ar_waterFinal_OFF[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent);
         // [/HNT]
      } //End if// //'dd=16 if
      if (gs_doneFirstDay) { //if// //'calculate plc relative to midday of day 1
         if (iter_refK < 0.000000001) // || iter_Counter == 0 // no longer appropriate to test for iter_Counter == 0 ... may be doing a seperate stress profile that refers to a saved refK, which will have been loaded at start of modelProgramMain
            tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1); //' done to avoid repeating this calculation
         else // if we haven't loaded a ref K it will be set to zero, and we need to fall back to the old method.. otherwise use refK to calculate PLC
         {
            tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / iter_refK); // PLCp calculated from refK if it exists
            if (tempDouble < 0.0) // if we're using the refK, it's possible for this to be negative briefly -- should be considered zero
               tempDouble = 0.0;
         }
         dSheet.Cells(rowD + dd, colD + dColF_End_PLCplant) = tempDouble; //'100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1) 'plc plant...prior timestep

                                                                          // no matter what, add it to the tally for calculating the mean over-season PLCp
         gs_ar_PLCSum[gs_yearIndex] = gs_ar_PLCSum[gs_yearIndex] + tempDouble; // yearly running total of PLC values
         gs_ar_PLCSum_N[gs_yearIndex] = gs_ar_PLCSum_N[gs_yearIndex] + 1; // total hours in GS
                                                                          // now test for highest PLC and hours > 85
         if (tempDouble > gs_ar_PLCp[gs_yearIndex])
            gs_ar_PLCp[gs_yearIndex] = tempDouble;
         if (tempDouble > 85.0)
            gs_ar_PLC85[gs_yearIndex] = gs_ar_PLC85[gs_yearIndex] + 1;

         tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1);
         dSheet.Cells(rowD + dd, colD + dColF_End_PLCxylem) = tempDouble; //'100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1) 'plc xylem...prior timestep
         if (tempDouble > gs_ar_PLCx[gs_yearIndex])
            gs_ar_PLCx[gs_yearIndex] = tempDouble;

         //dSheet.Cells(rowD + dd, colD + dColF_End_PLCplant) = (100.0 * (1.0 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1)); //'plc plant...prior timestep
         //dSheet.Cells(rowD + dd, colD + dColF_End_PLCxylem) = (100.0 * (1.0 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1)); //'plc xylem...prior timestep

         // [HNT] keep track of in-season input
         if (gs_inGrowSeason)
         {
            // if we've done the first GS day and we're in the growing season then it's "this year" during GS
            gs_ar_waterInput_GS[gs_yearIndex] += dSheet.Cells(rowD + dd, colD + dColF_End_input);
         }
         else
         {
            // if we've done the first GS day but we're NOT in the growing season, it's the winter following GS.
            // this is considered NEXT YEAR's off-season input! First check if we already have a value for the FINAL water for THIS YEAR, because if it's zero then
            // this is the first timestep of the winter and we need to store it
            if (gs_ar_waterFinal_GS[gs_yearIndex] <= 0.0 && gs_ar_waterInitial_OFF[gs_yearIndex + 1] <= 0.0) // ok to use == or <= with double here because this will be memset to 0 if it hasn't been set
            {
               gs_ar_waterFinal_GS[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent); // ok to overshoot by 1 hour, I think (instead of using dd - 1)
               gs_ar_waterInitial_OFF[gs_yearIndex + 1] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent); // also the off-season initial for next year
            }
            else // otherwise, we're in the middle of NEXT YEAR's off season ... note this +1 on an array index is super lazy and bad. Make sure to never run exactly this array size # of years
            {
               gs_ar_waterInput_OFF[gs_yearIndex + 1] += dSheet.Cells(rowD + dd, colD + dColF_End_input); // add the stored input to the input tally, for NEXT YEAR
            }
         }
         // [/HNT]
      } //End if// //'dd>16 if
      else if (!gs_inGrowSeason)// have NOT done first day and are NOT in growing season
      {
         // we must be in the pre-GS winter of what we called the NEXT year above... so gs_yearIndex is now that year, and we add to the off-season input for THIS year
         gs_ar_waterInput_OFF[gs_yearIndex] += dSheet.Cells(rowD + dd, colD + dColF_End_input);
      }
}

// Rhizosphere hydrodynamics
/* Rhizosphere E(P)
generate soil E(P) global curve--only need to do this once
*/
void MainProgram::get_rhizoPcrit(){
    for (z = 1; z <= layers; z++) {//z = 1 To layers
        p1 = 0;
        k = 1;
        e = 0; //flow integral
        erh[z][0] = 0; //first array position is layer number
        krh[z][0] = kmaxrh[z];
        do{
            p2 = p1 + pinc;
            qtrapvg(p1, p2, s);
            e = e + s;
            erh[z][k] = e;
            x = p2;
            krh[z][k] = get_vg(x); //instantaneous k at upper limit of integration = derivative of flow integral (fundamental theorem of calculus)
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000){
               break; //Then Exit Do //avoid crashing for extreme vc//s
            }
        } while (!(krh[z][k - 1] < kmin));//Loop Until krh(z, k - 1) < kmin
        pcritrh[z] = p2; //end of line for rhizo element z
    } //endfor// z
} //endsub//

/* Rhizosphere flow
Calculates rhizosphere flow using global E(P) curve
*/
void MainProgram::get_rhizoflow(){
   plow = int(p1 / pinc); //'pressure index below target
   k = int(p1 / pinc);
   elow = erh[z][k]; //'e below target
   klow = krh[z][k];
   ehigh = erh[z][k + 1]; //'e above target
   khigh = krh[z][k + 1];
   plow = plow * pinc; //'convert index to pressure
   estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //'linear interpolation of starting e
   klower = (p1 - plow) / pinc * (khigh - klow) + klow; //'linear interpolation of K(P)at lower limit of integration
   plow = int(p2 / pinc); //'pressure index below target
   k = int(p2 / pinc);
   elow = erh[z][k]; //'e below target
   klow = krh[z][k];
   ehigh = erh[z][k + 1]; //'e above target
   khigh = krh[z][k + 1];
   plow = plow * pinc; //'convert index to pressure
   efinish = (p2 - plow) / pinc * (ehigh - elow) + elow; //'linear interpolation of finishing e
   kupper = (p2 - plow) / pinc * (khigh - klow) + klow; //'linear interpolation of K(P) at upper limit of flow integration
   flow = efinish - estart; //'e upstream flow
}

// Hydraulics module 
// TO-DO: paste these functions to a separate file as its own class

/*Weibull parameter c
Calculates c parameter from p12 and p50
*/
double MainProgram::get_cweibull(double &p12, double &p50) {
   // c = (log(log(0.88)/log(0.50)))/(log(p12)-log(p50))
   double px1;
   double px2;
   px1 = p12*-1; // absolute value of p12
   px2 = p50*-1; // absolute value of p50
   return (log(log(0.88)/log(0.5)))/(log(px1)-log(px2));
}

/*Weibull parameter b
Calculates b parameter from parameter c and p12
*/
double MainProgram::get_bweibull(double &p12, double &c) {
    // b = p12/(-log(0.88))^(1/c)
    double px1;
    px1 = p12*-1; // absolute value of p12
    return px1/(pow(-log(0.88),(1/c)));
}

/*Weibull function for stem and leaves
Fits weibull function for a given water potential x.
It requires parameters b & c, and ksat
*/
double MainProgram::get_weibullfit(double &x, double &b, double &c, double &ksat) {
   //wb = ksat * Exp(-((x / b) ^ c))
   return ksat * exp(-(pow((x / b), c)));
}

// ROOT specific functions
/*Weibull function for root elements*/
double MainProgram::get_weibullfitroot(double &x){
   //wb = ksat * Exp(-((x / b) ^ c))
   return ksatr[z] * exp(-(pow((x / root_b), root_c)));
}

///*Functions to integrate root element z weibull*/
/*Root element z integration function 1:
*/
void MainProgram::trapzdwbr(double &p1, double &p2, double &s, long &t) {
    if (t == 1){
        s = 0.5 * (p2 - p1) * (get_weibullfitroot(p1) + get_weibullfitroot(p2));
        it = 1;
    } else {
        tnm = it;
        del = (p2 - p1) / tnm;
        x = p1 + 0.5 * del;
        sum = 0;
        for (j = 1; j <= it; j++){
            sum = sum + get_weibullfitroot(x);
            x = x + del;
        }
        s = 0.5 * (s + (p2 - p1) * sum / tnm);
        it = 2 * it;
    }
    //[HNT]
    //testcounterIntRoot = testcounterIntRoot + 1
    //[\HNT]
}

/*Root element integration function 2:
Evaluates accuracy of root element z integration
*/
void MainProgram::qtrapwbr(double &p1, double &p2, double &s){
    olds = -1; //'starting point unlikely to satisfy if statement below
    for (t = 1; t <= f; t++){
        trapzdwbr(p1, p2, s, t);
        if (std::abs(s - olds) <= epsx * std::abs(olds)){
            return;
        }
        olds = s;
    }
}

/* Root E(P) curves
Generates fresh global E(P) curve for the element(erases history)
*/
void MainProgram::get_rootPcrit() { 
    //clear arrays from earlier calls
    memset(er, 0, sizeof(er));//Erase er
    memset(kr, 0, sizeof(kr));//Erase kr
    //do root elements
    for (z = 1; z <= layers; z++){//z = 1 To layers
        kr[z][0] = ksatr[z]; //kr is instantaneous K from weibull function
        //kminr(z) = ksatr(z)
        p1 = 0;
        k = 1;
        e = 0; //value of integral
        er[z][0] = 0; //first array position is layer number
        do{
            p2 = p1 + pinc;
            qtrapwbr(p1, p2, s);
            e = e + s;
            er[z][k] = e;
            x = p2;
            kr[z][k] = get_weibullfitroot(x); //weibull k at upper limit of integration = derivative of flow integral
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000){
                break;
            }
            //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
        } while (!(kr[z][k - 1] < kmin));//Loop Until kr(z, k - 1) < kmin
        pcritr[z] = p2; //end of line for element z
    } //endfor// z
} //endsub//

/* Root E(P) virgin curves
Generates fresh global E(P) virgin curve for the element(erases history)
*/
void MainProgram::get_rootPcrit_v() { //generates fresh global E(P) curve for the element(erases history)
    //clear arrays from earlier calls
    memset(er_v, 0, sizeof(er_v));//Erase er_v
    memset(kr_v, 0, sizeof(kr_v));//Erase kr_v
    //do root elements
    for (z = 1; z <= layers; z++){//z = 1 To layers
        kr_v[z][0] = ksatr[z]; //kr_v is instantaneous K from weibull function
        //kminr(z) = ksatr(z)
        p1 = 0;
        k = 1;
        e = 0; //value of integral
        er_v[z][0] = 0; //first array position is layer number
        do{
            p2 = p1 + pinc;
            qtrapwbr(p1, p2, s);
            e = e + s;
            er_v[z][k] = e;
            x = p2;
            kr_v[z][k] = (x); //weibull k at upper limit of integration = derivative of flow integral
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000){
                break;
            }//If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
        } while (!(kr_v[z][k - 1] < kmin));//Loop Until kr_v(z, k - 1) < kmin
        pcritr[z] = p2; //end of line for element z
    } //endfor// z
} //endsub//

/* Root layer flow
Gets flow through root using global E(P) curve
*/
void MainProgram::get_rootflow(){
   plow = int(p1 / pinc); //'pressure index below target
   k = int(p1 / pinc);
   elow = er[z][k]; //'e below target
   klow = kr[z][k];
   ehigh = er[z][k + 1]; //'e above target
   khigh = kr[z][k + 1];
   plow = plow * pinc; //'convert index to pressure
   estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //'linear interpolation of starting e
   klower = (p1 - plow) / pinc * (khigh - klow) + klow; //'linear interpolation of starting K(P)(lower limit of integration)
   plow = int(p2 / pinc); //'pressure index below target
   k = int(p2 / pinc);
   elow = er[z][k]; //'e below target
   klow = kr[z][k];
   ehigh = er[z][k + 1]; //'e above target
   khigh = kr[z][k + 1];
   plow = plow * pinc; //'convert index to pressure
   efinish = (p2 - plow) / pinc * (ehigh - elow) + elow; //'linear interpolation of finishing e
   kupper = (p2 - plow) / pinc * (khigh - klow) + klow; //'linear interpolation of ending K(P)(upper limit of integration)
   flow = efinish - estart; //'e upstream flow
}

/* LU decomposition
does LU decomposition on the jacobian prior to solution by lubksb
*/
void MainProgram::ludcmp() {
    d = 1;
    for (i = 1; i <= unknowns; i++){//i = 1 To unknowns
        aamax = 0;
        for (j = 1; j <= unknowns; j++){//j = 1 To unknowns
            if (std::abs(jmatrix[i][j]) > aamax){
                aamax = std::abs(jmatrix[i][j]);
            }
        }
        if (aamax == 0){
            return;
        }
        vv[i] = 1 / aamax;
    }
    for (j = 1; j <= unknowns; j++){//j = 1 To unknowns
        for (i = 1; i < j; i++){//i = 1 To j - 1
            sum = jmatrix[i][j];
            for (k = 1; k < i; k++){//k = 1 To i - 1
                sum = sum - jmatrix[i][k] * jmatrix[k][j];
            }
            jmatrix[i][j] = sum;
        }
        aamax = 0;
        for (i = j; i <= unknowns; i++){//i = j To unknowns
            sum = jmatrix[i][j];
            for (k = 1; k < j; k++){//k = 1 To j - 1
                sum = sum - jmatrix[i][k] * jmatrix[k][j];
            }
            jmatrix[i][j] = sum;
            dum = vv[i] * std::abs(sum);
            if (dum > aamax){
                imax = i;
                aamax = dum;
            }
        }
        if (j != imax){// j <> imax
            for (k = 1; k <= unknowns; k++){//k = 1 To unknowns
                dum = jmatrix[imax][k];
                jmatrix[imax][k] = jmatrix[j][k];
                jmatrix[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (jmatrix[j][j] == 0){
            jmatrix[j][j] = 1E-25;
        }
        if (j != unknowns){
            dum = 1 / jmatrix[j][j];
            for (i = j + 1; i <= unknowns; i++){//i = j + 1 To unknowns
                jmatrix[i][j] = jmatrix[i][j] * dum;
            }
        }
    }
    memset(vv, 0, sizeof(vv));//Erase vv
}

/* Delta p's solver
solves the decomposed jacobian delta p's
*/
void MainProgram::lubksb(){
    ii = 0;
    for (i = 1; i <= unknowns; i++){//i = 1 To unknowns
        ll = int(indx[i]); //'indx array comes from ludcmp
        sum = func[ll]; //'the func array input is the right-hand vector
        func[ll] = func[i];
        if (ii != 0){
            for (j = ii; j < i; j++){//j = ii To i - 1
               sum = sum - jmatrix[i][j] * func[j];
            }
        } else {
            if (sum != 0){
                ii = i;
            }
        }
        func[i] = sum;
    }
    for (i = unknowns; i >= 1; i--){//i = unknowns To 1 Step -1
        sum = func[i];
        for (j = i + 1; j <= unknowns; j++) {//j = i + 1 To unknowns
            sum = sum - jmatrix[i][j] * func[j];
        }
        func[i] = sum / jmatrix[i][i];
    }
}

/* Newton Rhapson Function
returns rhizosphere pressures and root pressure, pr, as function of pd's and e
*/
void MainProgram::newtonrhapson(){
    //'prinitial = pr //record the original guess
    weird = 0; //tracks pr estimate
    check = 0; //restart loop counter
    do {//loop to reset guesses
        //'restore layer functioning
        for (z = 1; z <= layers; z++){
            if (kminroot[z] != 0){
                layer[z] = 0; //make sure to start with all layers functioning
                layerfailure[z] = "no failure"; //reset
            }
        }
        failspot = "no failure";

        if (weird == 1){//reset guesses
            double rFloat = (double)rand() / (double)RAND_MAX;
            k = int((layers - 1 + 1) * rFloat + 1);
            pr = pd[k]; //random choice of pd
            
            // can enable this is running into erroneous solutions -- but allowing these to vary instead of resetting results in more frequent solutions in my experience
            if (false){ // alternatively could randomize them properly
                for (k = 1; k <= layers; k++){//reset prhz(k)
                    prh[k] = pd[k];
                }
            }
            weird = 0; //reset cutoff
        } //end reset guesses loop
        
        check = check + 1; //number of restarts
        ticks = 0; //convergence counter
        
        do {//loop to seek convergence
            ticks = ticks + 1;
            if (ticks > 1000){
                weird = 1;
                std::cout << "NR ticks exceeded 1000 -- setting weird=1 for retry. Pinc too high? dd = " << dd << std::endl;
            }
            //'get top row of matrix and right-hand func vector
            //'zero out the jacobian first
            for (k = 1; k <= unknowns; k++){
                for (j = 1; j <= unknowns; j++){
                    jmatrix[k][j] = 0;
                }
            }
            //'fill up four arrays:
            //'func(i) is zero flow function...the right-hand flow vector
            //'dfrhdprh(i) is partial derivative of frh(i) for prh(i)(rhizo pressure)...the diagonal of the jacobian
            //'dfrhdpr(i) is partial derivative of frh(i)for pr (root pressure)...the last column of the jacobian
            //'dfrdprh(i) is partial derivative of fr for prh(i)...the last row of the jacobian
            frt = 0; //this is the last row of right-hand flow vector
            dfrdpr = 0; //this is lower right-hand partial for jacobian
            for (z = 1; z <= layers; z++){
               if (pd[z] >= pcritrh[z] || prh[z] >= pcritrh[z]){
                    layer[z] = 1; //layer's just gone out of function
                    layerfailure[z] = "rhizosphere";
                    //weird = 1; //could be result of non-convergence
                }
                if (layer[z] == 0){//it's functional
                    p1 = pd[z]; //p1 is not a guess
                    p2 = prh[z]; //prh[z] IS a guess and it is initially set prior to the RN routine
                    get_rhizoflow(); //gets flows through rhizosphere element from p2-p1 and E(p)curve; gets K's at p1 and p2 as well
                    func[z] = flow;
                    dfrhdprh[z] = kupper;
                }
                if (prh[z] >= pcritr[z] || pr >= pcritr[z]){
                    layer[z] = 1; //layer's just gone out of function
                    layerfailure[z] = "root";
                    //weird = 1; //could be result of non-convergence
                }
                if (layer[z] == 0) {//it's functional
                    p1 = prh[z]; //now re-set p1 to prh[z]...the guess
                    p2 = pr; //guess must come from before NR routine?
                    get_rootflow(); //gets flows through root element, and K's at upper and lower bounds
                    func[z] = func[z] - flow;
                    dfrhdprh[z] = dfrhdprh[z] + klower;
                    dfrhdpr[z] = -kupper;
                    dfrdprh[z] = -klower;
                }
                if (layer[z] == 1){ //layer's out of function...zero out values
                    func[z] = 0;
                    dfrhdprh[z] = 1E-25;
                    dfrdprh[z] = 1E-25;
                    dfrhdpr[z] = 1E-25;
                    kupper = 1E-25;
                    flow = 0;
                }
                dfrdpr = dfrdpr + kupper;
                frt = frt + flow;
            }
            frt = frt - e;
            //'now load jacobian
            for (k = 1; k <= layers; k++){
                jmatrix[k][unknowns] = dfrhdpr[k]; //last column with dFrh/dPr partials
            }
            for (k = 1; k <= layers; k++){
                jmatrix[unknowns][k] = dfrdprh[k]; //last row with dFr/dPrh partials
            }
            for (k = 1; k <= layers; k++){
                jmatrix[k][k] = dfrhdprh[k]; //diagonal of dFrh/dPrh partials
            }
            jmatrix[unknowns][unknowns] = dfrdpr; //lower right corner with dFr/dPr partial
            func[unknowns] = frt; //last position in right-hand flow vector

            //'ok, jacobian and righthand vector are loaded
            //'test for total failure
            sum = 0;
            for (k = 1; k <= layers; k++){
                sum = sum + layer[k];
            }
            if (sum == layers){//total failure
                failspot = "belowground";
                weird = 1; //trigger a restart
            }
            //'test for flow conservation (steady-state)
            threshold = 0;
            for (k = 1; k <= unknowns; k++){//k = 1 To unknowns
               threshold = threshold + std::abs(func[k]);
            }
            if (ticks == 1){
                initialthreshold = threshold;
            }
            //'remember to replace "n" with "unknowns" in ludcmp and lubksb
            ludcmp(); //numerical recipe for doing LU decomposition of jacobian prior to solving
            lubksb(); //solves the decomposed jacobian for delta p's
            //'print out solution vector of pressures
            //'revise unknown pressures
            for (k = 1; k <= layers; k++){//k = 1 To layers
                prh[k] = prh[k] - func[k]; //NOTE lubksb replaces original right-side func()vector with the solution vector
            }
            pr = pr - func[unknowns];
            //'check for jumping lower bound
            for (k = 1; k <= layers; k++){//k = 1 To layers
                if (prh[k] < 0){
                    prh[k] = 0;
                }
            }
            if (pr < 0){
                pr = 0;
            }
            //'if pr > pcritr Then
            //'pr = prinitial
            //'weird = 1 //trigger a re-start
            //'}
            if (ticks > 1){//check for non convergence
                if (threshold > initialthreshold){
                    weird = 1;//pr is spiraling, restart NR with new guesses
                }
                if (pr >= pcrits){
                    weird = 1;
                }
            }
        } while (!(threshold < 0.02 || weird == 1)); //Loop Until threshold < 0.01 Or weird = 1 //weird = 1 restarts NR
    } while (!((threshold < 0.02 && weird == 0) || check > 500));

    if (check > 500){
        // disable this output if it's causing too much spam
        //std::cout << "NR Failure " << threshold << " check = " << check << " dd = " << dd << " watercontent = " << waterold * 1000.0 << " weird = " << weird << std::endl;
        // keep track of the frequency of NR failures
        gs_ar_nrFailConverge[gs_yearIndex]++; // non-convergent failure
        gs_ar_nrFailConverge_Water[gs_yearIndex] += waterold * 1000.0;
        if (waterold * 1000.0 > gs_ar_nrFailConverge_WaterMax[gs_yearIndex]){
            gs_ar_nrFailConverge_WaterMax[gs_yearIndex] = waterold * 1000.0;
        }
    }

    //Loop Until threshold < 0.01 And weird = 0 Or check > 500 //give up after 2000 restarts
    //'if check >= 500 Then Stop
    //final step -- recheck the layers
    for (z = 1; z <= layers; z++){
        if (kminroot[z] != 0){
            layer[z] = 0; //make sure to start with all layers functioning
            layerfailure[z] = "no failure"; //reset
        }
        if (pd[z] >= pcritrh[z] || prh[z] >= pcritrh[z]){
            layer[z] = 1; //layer's just gone out of function
            layerfailure[z] = "rhizosphere";
        }
        if (prh[z] >= pcritr[z] || pr >= pcritr[z]){
            layer[z] = 1; //layer's just gone out of function
            layerfailure[z] = "root";
        }
    }
}

// STEM specific funtions
/* Functions to integrate stem element*/
/* Stem element integration
*/
void MainProgram::trapzdwbs(double &p1, double &p2, double &s, long &t) {
    if (t == 1){
        s = 0.5 * (p2 - p1) * (get_weibullfit(p1,stem_b,stem_c,ksats) + get_weibullfit(p2,stem_b,stem_c,ksats));
        it = 1;
    } else {
        tnm = it;
        del = (p2 - p1) / tnm;
        x = p1 + 0.5 * del;
        sum = 0;
        for (j = 1; j <= it; j++) {
            sum = sum + get_weibullfit(x,stem_b,stem_c,ksats);
            x = x + del;
        }
        s = 0.5 * (s + (p2 - p1) * sum / tnm);
        it = 2 * it;
    }
}

/* Integration accuracy
Evaluates accuracy of stem element integration
*/
void MainProgram::qtrapwbs(double &p1, double &p2, double &s) {
    olds = -1; //starting point unlikely to satisfy if statement below
    for (t = 1; t <= f; t++){
        trapzdwbs(p1, p2, s, t);
        if (std::abs(s - olds) <= epsx * std::abs(olds)){
            return;
        }
        olds = s;
    }
}

/* Stem E(P) curve
Generates fresh global E(P) curve for the stem element(erases history)
*/
void MainProgram::get_stemPcrit(bool vCurve = false){
    double *es_ptr = es;
    if (vCurve){
        memset(es_v, 0, sizeof(es_v));
        es_ptr = es_v;
    } else {
        memset(es, 0, sizeof(es));
    }
    //memset(es, 0, sizeof(es));//Erase es //eliminate values from previous calls
    p1 = 0;
    k = 1;
    e = 0; //value of integral
    es_ptr[0] = 0;
    do{
        p2 = p1 + pinc;
        qtrapwbs(p1, p2, s);
        e = e + s;
        es_ptr[k] = e;
        x = p2;
        ksh = get_weibullfit(x,stem_b,stem_c,ksats); //weibull k
        p1 = p2; //reset p1 for next increment
        k = k + 1;
        if (k == 100000){
            break;
        }//If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
    } while (!(ksh < kmin));//Loop Until ksh < kmin
    pcrits = p2;//end of line for element z
} //endsub//

/* Stem hydraulic tension
gets stem pressure and conductance from pr and e.
*/
void MainProgram::get_stem() {
    //'start with stem
    p1 = pr + pgrav; //add gravity drop before integration
    plow = int(p1 / pinc); //pressure index below target
    k = int(p1 / pinc);
    elow = es[k]; //e below target
    ehigh = es[k + 1]; //e above target
    plow = plow * pinc; //convert index to pressure
    estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //linear interpolation of starting e
    efinish = estart + e;
    j = k;
    do {//find efinish
        j = j + 1;
        if (es[j] == 0){
            test = 1;
            failspot = "stem";
            return;
        }
    } while (!(es[j] > efinish)); //Loop Until es(j) > efinish
    ehigh = es[j];
    elow = es[j - 1];
    p2 = ((efinish - elow) / (ehigh - elow)) * pinc + pinc * (j - 1); //pstem
    ps = p2; //ps is downstream stem pressure
    if (ps >= pcrits) {
        test = 1;
        failspot = "stem";
        return;
    }
}

// LEAF specific functions
/* Function to integrate leaf element
*/
void MainProgram::trapzdwbl(double &p1, double &p2, double &s, long &t){
    if (t == 1){
        s = 0.5 * (p2 - p1) * (get_weibullfit(p1,leaf_b,leaf_c,ksatl) + get_weibullfit(p2,leaf_b,leaf_c,ksatl));
        it = 1;
    } else {
        tnm = it;
        del = (p2 - p1) / tnm;
        x = p1 + 0.5 * del;
        sum = 0;
        for (j = 1; j <= it; j++){
            sum = sum + get_weibullfit(x,leaf_b,leaf_c,ksatl);
            x = x + del;
        }
        s = 0.5 * (s + (p2 - p1) * sum / tnm);
        it = 2 * it;
    }
}

/* Integration accuracy
Evaluates accuracy of leaf element integration
*/
void MainProgram::qtrapwbl(double &p1, double &p2, double &s){
    olds = -1; //starting point unlikely to satisfy if statement below
    for (t = 1; t <= f; t++){
        trapzdwbl(p1, p2, s, t);
        if (std::abs(s - olds) <= epsx * std::abs(olds)){
            return;
        }
        olds = s;
    }
}

/* Leaf E(P) curve
Generates fresh global E(P) curve for the leaf element(erases history)
*/
void MainProgram::get_leafPcrit(bool vCurve = false) {
    double *el_ptr = el;
    if (vCurve) {
        memset(el_v, 0, sizeof(el_v));
        el_ptr = el_v;
    } else {
        memset(el, 0, sizeof(el));
    }
    //memset(el_ptr, 0, sizeof(el_ptr));//Erase el_ptr //eliminate values from previous calls
    p1 = 0;
    k = 1;
    e = 0; //value of integral
    el_ptr[0] = 0;
    do {
        p2 = p1 + pinc;
        qtrapwbl(p1, p2, s);
        e = e + s;
        el_ptr[k] = e;
        x = p2;
        ksh = get_weibullfit(x,leaf_b,leaf_c,ksatl); //weibull k
        p1 = p2; //reset p1 for next increment
        k = k + 1;
        if (k == 100000){
            break;
        }//If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
    } while (!(ksh < kmin)); //Loop Until ksh < kmin
    pcritl = p2; //end of line for element z
} //endsub//

/* Leaf hydraulic tension
//gets leaf pressure from stem pressure and e
*/
void MainProgram::get_leaf(){
    p1 = ps;
    plow = int(p1 / pinc); //pressure index below target
    k = int(p1 / pinc);
    elow = el[k]; //e below target
    ehigh = el[k + 1]; //e above target
    plow = plow * pinc; //convert index to pressure
    estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //linear interpolation of starting e
    efinish = estart + e;
    j = k;
    do { //find efinish
        j = j + 1;
        if (el[j] == 0) { //= Empty
            test = 1;
            failspot = "leaf";
            return;
        }
    } while (!(el[j] > efinish));//Loop Until el(j) > efinish
    ehigh = el[j];
    elow = el[j - 1];
    p2 = ((efinish - elow) / (ehigh - elow)) * pinc + pinc * (j - 1); //pleaf
    pl = p2; //pl is leaf pressure
    if (pl >= pcritl){
        test = 1;
        failspot = "leaf";
    }
}

// WHOLE PLANT
/*Get pcrits*/
void MainProgram::componentpcrits(){ //'gets pcrits
    get_rhizoPcrit();
    get_rootPcrit(); //'gets root element curves
    get_stemPcrit(); //'gets stem element curve
    pcrits;
    get_leafPcrit(); //'gets leaf element curve
    pcritl;

    get_rootPcrit_v(); //erases history for md solution
    get_stemPcrit(true);
    get_leafPcrit(true);

    memset(ter, 0, sizeof(ter));
    memset(tkr, 0, sizeof(tkr));
    memset(tes, 0, sizeof(tes));
    memset(tel, 0, sizeof(tel));
}

/* Composite curve
stores composite E(P)curve and the element conductances
*/
void MainProgram::compositecurve() {
    elayer[0][p] = 0; //'no root mediated flow in topmost layer
    for (z = 1; z <= layers; z++){//z = 1 To layers
         if (layer[z] == 0) { //if//
            prhizo[z][p] = prh[z];
            p1 = prh[z];
            p2 = pr;
            get_rootflow();
            elayer[z][p] = flow; //'flow through layer
            if (flow != 0){
               kroot[z][p] = std::abs(elayer[z][p] / (pr - prhizo[z][p]));
            }
            if (flow == 0) { //if//
                if (refilling == true) { //if// //'for refilling, starting point is always weibull
                    x = pd[z];
                    kroot[z][p] = get_weibullfitroot(x);
                } //End if//
                if (refilling == false){
                    kroot[z][p] = kminroot[z];
                }
            } //End if//
         } //End if//
         if (layer[z] == 1) { //if//
            elayer[z][p] = 0; //'no flow
            if (layerfailure[z] == "root") { //if// //'root element has failed
                kroot[z][p] = 0; //'total cavitation in root
                prhizo[z][p] = pd[z]; //'rhizosphere pressure returns to the predawn value
            } //End if//
            if (layerfailure[z] == "rhizosphere") { //if// //'rhizosphere element has failed
                x = pr;
                kroot[z][p] = get_weibullfitroot(x); //'root element conductance = instantaneous conductance from weibull curve at pr
                prhizo[z][p] = pcritrh[z];
            } //End if//
        } //End if//
    } //for//z
    proot[p] = pr;
    pstem[p] = ps;
    pleaf[p] = pl;
    if (e > 0) { //if//
        kleaf[p] = e / (pl - ps); //'leaf element conductance
        kstem[p] = e / (ps - pr - pgrav); //'stem element conductance subtracting extra gravity drop
        kplant[p] = e / (pl - pleaf[0] - pgrav); //'whole plant k, subtracting extra gravity drop
    } else {
        if (refilling == false) { //if//
            kleaf[p] = kminleaf;
            kstem[p] = kminstem;
            kplant[p] = kminplant;
        } //End if//
        if (refilling == true) { //if//
            x = pl;
            kleaf[p] = get_weibullfit(x,leaf_b,leaf_c,ksatl);
            x = ps;
            kstem[p] = get_weibullfit(x,stem_b,stem_c,ksats);
            //'kplant[p]=??? ignore kplant setting for e=0...probably not important
        } //End if//
        //'dedp[p] = kminplant
    } //End if//
    eplant[p] = e; //'total flow
    if (p > 0) { //if//
        if (pleaf[p] - pleaf[p - 1] == 0) { //if//
            test = 1;
        } else {
            dedp[p] = einc / (pleaf[p] - pleaf[p - 1]); //'dedp=instantaneous K of system
            dedpf[p] = dedp[p] / dedp[1]; //'fractional canopy conductance
        } //End if//
    } //End if//
    pcritsystem = pl;
    ecritsystem = e;
    total = p;
}

/* Store history
stores historical element e(P) curves
*/
void MainProgram::storehistory(){
    /*memset(ter, 0, sizeof(ter));
    memset(tkr, 0, sizeof(tkr));
    memset(tes, 0, sizeof(tes));
    memset(tel, 0, sizeof(tel));*/

    /*memcpy(ter, er, sizeof(ter));
    memcpy(tkr, kr, sizeof(tkr));
    memcpy(tes, es, sizeof(tes));
    memcpy(tel, el, sizeof(tel));

    memcpy(er, er_v, sizeof(er));
    memcpy(kr, kr_v, sizeof(kr));
    memcpy(es, es_v, sizeof(es));
    memcpy(el, el_v, sizeof(el));*/

    for (z = 1; z <= layers; z++){
        k = 0;
        do{
            ter[z][k] = er[z][k];
            tkr[z][k] = kr[z][k];
            er[z][k] = er_v[z][k];
            kr[z][k] = kr_v[z][k];
            k = k + 1;
            if (k == 100000){
               break;
            }
        } while (er[z][k] != 0 || er_v[z][k] != 0 || ter[z][k] != 0);//Loop Until er(z, k) = Empty
    }
    k = 0;
    do{
        tes[k] = es[k];
        es[k] = es_v[k];
        k = k + 1;
        if (k == 100000){
            break;
        }
    } while (es[k] != 0 || es_v[k] != 0 || tes[k] != 0); //Loop Until es(k) = Empty
    k = 0;
    do{
        tel[k] = el[k];
        el[k] = el_v[k];
        k = k + 1;
        if (k == 100000){
            break;
        }
    } while (el[k] != 0 || el_v[k] != 0 || tel[k] != 0); //Loop Until el(k) = Empty
    for (z = 1; z <= layers; z++){
        tlayer[z] = layer[z];
        tlayerfailure[z] = layerfailure[z];
    }
    for (z = 1; z <= layers; z++){// 'get all layers functioning to start with
        if (layerfailure[z] == "ok"){
            layer[z] = 0;
        } else if (layerfailure[z] == "rhizosphere"){
            layer[z] = 0;
        } else if (layerfailure[z] == "root" && kminroot[z] == 0){
            layer[z] = 1; // 'take out root layer if failed at start of composite curve
        } else{
            layer[z] = 0; //'still functional
        }
    }
}

/* Gets history
restores historical element e(P) curves
*/
void MainProgram::gethistory(){
    for (z = 1; z <= layers; z++){
        k = 0;
        do{
            er[z][k] = ter[z][k];
            kr[z][k] = tkr[z][k];
            k = k + 1;
            if (k == 100000){
                break;
            }
        } while (ter[z][k] != 0 || er[z][k] != 0);//Loop Until ter[z][k] = Empty
    }
    k = 0;
    do{
        es[k] = tes[k];
        k = k + 1;
        if (k == 100000){
            break;
        }
    } while (tes[k] != 0 || es[k] != 0);//Loop Until tes[k] = Empty
    k = 0;
    do{
        el[k] = tel[k];
        k = k + 1;
        if (k == 100000){
            break;
        }
    } while (tel[k] != 0 || el[k] != 0);//Loop Until tel[k] = Empty
    //'re-set failure status (this at the critical point)
    for (z = 1; z <= layers; z++){
        layer[z] = tlayer[z];
        layerfailure[z] = tlayerfailure[z];
    }
}

/* Update E(P) curves
resets element E(P) curves
*/
void MainProgram::updatecurves(){
    //'if k<kmin, re-assign e//'s on element curve by back-calculating
    for (z = 1; z <= layers; z++){ //z = 1 To layers
        if (true){
            if (kroot[z][halt] < kminroot[z]){
                kminroot[z] = kroot[z][halt];
                phigh = int(proot[halt] / pinc) + 1; //'pressure datum just above the target
                for (k = phigh; k >= 0; k--) {//k = phigh To 0 Step -1 //'back-calculate e//'s
                    er[z][k] = er[z][phigh] - kminroot[z] * pinc * (phigh - k);
                    kr[z][k] = kminroot[z]; //'back-calculate KR(Z,K) too for roots (not stem or leaves)
                } //EndFor  k
            } //EndIf//
        }
    } //EndFor  z
    if (kstem[halt] < kminstem){
        kminstem = kstem[halt];
        phigh = int(pstem[halt] / pinc) + 1;
        for (k = phigh; k >= 0; k--){ //k = phigh To 0 Step -1 //'back-calculate e//'s
            es[k] = es[phigh] - kminstem * pinc * (phigh - k);
        } //EndFor  k
    } //EndIf//
    if (kleaf[halt] < kminleaf){
        kminleaf = kleaf[halt];
        phigh = int(pleaf[halt] / pinc) + 1;
        for (k = phigh; k >= 0; k--){//k = phigh To 0 Step -1 //'back-calculate e//'s
            el[k] = el[phigh] - kminleaf * pinc * (phigh - k);
        } //EndFor  k
    } //EndIf//
    //'if kplant[halt] < kminplant Then kminplant = kplant[halt]NOTE: kplant CAN go up because of rhizosphere recovery!
}

/* Canopy pressure
computes carbon-based middays; pressures and cost curve from virgin supply function, gas exchange from historical values
store history, get MD and cost function from virgin curve
*/
void MainProgram::canopypressure(){  
    if (ecritsystem == 0){
        k = 0;
        goto tenMarker; //don't bother to find midday
    }
    storehistory(); //stores xylem element curves and failure status, resets layers to full functioning
    //rootcurves(); //erases history for md solution
    //stemcurve();
    //leafcurve();
    sum = 0;
    t = 0;
    for (k = 1; k <= layers; k++) {//assign source pressures, set layer participation
        if (layer[k] == 0){
            prh[k] = pd[k]; //initial guess of unknown rhizosphere pressures
            sum = sum + pd[k];
        } else {
            prh[k] = pcritr[k];
            t = t + 1;
        }
    }
    if (t < layers) {
        pr = sum / (layers - t); //set unknown proot to average pd
    } else {
        failure = 1; //system is critical
        return;
    }
    test = 0; //=1 if stem or leaf fails
    //virgin gain function params
    psynmaxmd = 0;
    psynmaxshmd = 0;
    //now loop through virgin risk curve
    e = -einc;
    p = -1;
    dedplmin = ksatp; //insures the kloss function is monotonic
    do {
        e = e + einc;
        p = p + 1;
        newtonrhapson(); //pd's already assigned...this solves for p's and e's in fingers of chicken foot

        if (check >= 500){
            //p = p - 1
            //std::cout << "NR failed on virgin curve at e = " << e << " on timestep dd = " << dd << std::endl;
            break; //gone as far as can go
        }
        get_stem(); //gets stem and leaf pressures
        get_leaf();
        pleafv[p] = pl; //pleaf from virgin curve
        if (p == 0){
            predawn = pl; //set the predawn...pl returned by "leaf" routine
            plold = pl;
        }
        if (p > 0){
            if ((pl - plold) == 0){
                break; //gone to failure
            }
            if (p == 1){
                dedplzero = einc / (pl - plold); //note: pl is returned by "leaf" routine
            }
            dedpl = einc / (pl - plold);
            if (dedpl < dedplmin){
               dedplmin = dedpl; //insure that kloss only goes up
            }
            klossv[p] = dedplzero - dedpl; //non-normalized kloss from virgin curve
            //dedpf(p) = dedpl / dedplzero //fractional k/kmax canopy from virgin curve (units don't matter!)
        }
        if (pl >= pcritsystem){
            break; //gone to failure
        }          //now get virgin A curve

        leaftempsmd(); //gets virgin sun layer leaf temperature from energy balance
        leaftempsshademd(); //gets virgin shade layer leaf temperature
        assimilationmd(); //gets virgin sun layer photosynthesis
        assimilationshademd(); //gets virgin shade layer photosynthesis
        //by now we have assigned psynmd(p) and psynshmd(p) and reset psynmaxes

        plold = pl;
        //e = e + einc
        //p = p + 1
    } while (!(test == 1 || p > 99900)); //loop to failure or out of "p"s
    //Loop Until test = 1 Or p > 99900 'loop to failure or out of "p"s
    //If check >= 2000 Then Exit Sub

    if (p <= 2){
        k = 0; //too close to critical, shut down system
        goto tenMarker;
    }

    klossv[0] = 0;
    maxkloss = dedplzero - dedpl; //maximum kloss...may not be kmax
    totalv = p - 1;
    klossv[totalv] = maxkloss;
    //now normalize klossv for virgin pleafv
    for (p = 0; p <= totalv; p++){
        klossv[p] = klossv[p] / maxkloss;
    }

    //now, find the middday from virgin risk and gain
    //First do for sun layer
    p = -1; //p is still index for virgin curve
    dpmax = 0;
    dpamax = -100;
    amaxmax = 0; //insures the gain function is monotonic

    dpamin = 0.0; // keep track of low values to avoid extreme negative profit curves producing a result
    do{
        p = p + 1;
        amaxfrac[p] = psynmd[p] / psynmaxmd; //this is the normalized revenue function from the virgin curve
        if (amaxfrac[p] > amaxmax){
            amaxmax = amaxfrac[p];
        } else {
            amaxfrac[p] = amaxmax;
        } //this insures that amaxfrac monotonically increases
        if (amaxfrac[p] < 0){
            amaxfrac[p] = 0; //no negative gains
        }
        if (klossv[p] < 0){
            klossv[p] = 0; //no negative risks
        }
        dpa[p] = amaxfrac[p] - klossv[p]; //profit, with revenue from virgin curve and cost from virgen one
        if (p < runmean - 1){
            rmean = 0;
        }
        if (p >= runmean - 1){ //get running mean
            sum = 0;
            for (i = p - runmean + 1; i <= p; i++){
                sum = sum + dpa[i];
            }
            rmean = sum / runmean; //the running mean
        }
        if (rmean < dpamin){
            dpamin = rmean;
        }
        if (rmean < 0){
            rmean = 0; //avoid negative rmean
        }
        if (rmean > dpamax){
            dpamax = rmean;
            md = pleafv[p]; //midday pressure for sun layer from virgin curves
        }//print out gain and cost
    } while (!(einc * p >= gmax * lavpd[p] || total == 0 || (rmean < dpamax / cutoff && p > runmean && p > 15) || klossv[p] > 0.9 || p >= totalv));

    //std::cout << "DPA MIN = " << dpamin << " DPA MAX = " << dpamax << std::endl;
    if (dpamin < 0.0 && dpamax > 0.0 && std::abs(dpamin) > dpamax){
        // the profit went more negative than positive, so reset mid-day to predawn
        md = pleafv[0];
    }
    // [HNT] debug
    if (!(rmean < dpamax / cutoff)){
        //std::cout << "Terminated sun layer opt without finding peak! At timestep dd = " << dd << std::endl;
        if (einc * p >= gmax * lavpd[p]){
            ;// std::cout << "Terminated sun layer opt without finding peak: end case 1 einc * p >= gmax * lavpd[p]" << std::endl;
        }
        if (total == 0){
            std::cout << "Terminated sun layer opt without finding peak: end case 2 total == 0" << std::endl;
        }
        if (p >= totalv){
            std::cout << "Terminated sun layer opt without finding peak: end case 3 exceeded end of virgin curve" << std::endl;
        }
        if (klossv[p] > 0.9){
            std::cout << "Terminated sun layer opt without finding peak: end case 4 klossv[p] > 0.9" << std::endl;
        }
    }
    // [/HNT]
    //while (!(einc * p >= gmax * lavpd[p] || total == 0 || rmean < dpamax / cutoff && p > runmean || klossv[p] > 0.9 || p >= totalv));
    //Loop Until einc * p >= gmax * lavpd(p) Or total = 0 Or rmean < dpamax / cutoff And p > runmean Or klossv(p) > 0.9 Or p >= totalv //loop until g maxed out or to failure...note e and g in kg hr-1

    dpasun = dpamax; // unused?

    //now do for shade layer
    if (psynmaxshmd == 0) { //shade layer's below light compensation point, don't open
        mdsh = pleafv[0]; //set midday = predawn
    } else {
        p = -1; //p is still index for virgin curve
        dpmax = 0;
        dpamax = -100;
        amaxmax = 0;
        dpamin = 0.0;
        do {
            p = p + 1;
            amaxfracsh[p] = psynshmd[p] / psynmaxshmd; //this is the normalized shade revenue function from the virgin curve
            if (amaxfracsh[p] > amaxmax){
                amaxmax = amaxfracsh[p];
            } else {
                amaxfracsh[p] = amaxmax;
            }//this insures that amaxfrac monotonically increases
            //now loop to find kloss from virgin curve that matches historical pleaf
            if (amaxfracsh[p] < 0){
                amaxfracsh[p] = 0; //no negative gains
            }
            if (klossv[p] < 0){
                klossv[p] = 0; //no negative risks
            }
            dpa[p] = amaxfracsh[p] - klossv[p]; //profit, with revenue from historical curve and cost from virgen one
            if (p < runmean - 1){
                rmean = 0;
            }
            if (p >= runmean - 1){ //get running mean
                sum = 0;
                for (i = p - runmean + 1; i <= p; i++){
                    sum = sum + dpa[i];
                }
                rmean = sum / runmean; //the running mean
            }
            if (rmean < dpamin){
                dpamin = rmean;
            }
            if (rmean < 0){
                rmean = 0; //avoid negative dpa
            }
            if (rmean > dpamax){
                dpamax = rmean;
                mdsh = pleafv[p]; //midday pressure for shade layer from virgin curve
            }
        } while (!(einc * p >= gmax * lavpdsh[p] || total == 0 || (rmean < dpamax / cutoff && p > runmean && p > 15) || klossv[p] > 0.9 || p >= totalv));
        //Loop Until einc * p >= gmax * lavpdsh[p] Or total = 0 Or rmean < dpamax / cutoff And p > runmean Or klossv[p] > 0.9 Or p >= totalv //loop until g maxed out or to failure...note e and g in kg hr-1
        if (dpamin < 0.0 && dpamax > 0.0 && std::abs(dpamin) > dpamax){
            // the profit went more negative than positive, so reset mid-day to predawn
            mdsh = pleafv[0];
        }
    } //psynmaxsh if

    k = -1;
    //Range("c17:f10000").ClearContents
    do {
        k = k + 1;
    } while (!(pleaf[k] >= md || pleaf[k] == 0));
    //Loop Until pleaf(k) >= md Or pleaf(k) = Empty //pleaf from historical curve must match pleaf from virgin curve
    transpiration = eplantl[k]; //all gas exchange values are from most recent historical values
    psynact = psyn[k];
    gcmd = gcanw[k]; //g for water in mmol
    lavpdmd = lavpd[k] * patm;
    cinc = cin[k];
    //If k > 1 Then deda = (eplantl(k) - eplantl(k - 1)) / (psyn(k) - psyn(k - 1))
    halt = k; //halt is index of midday datum
    //now do shade layer
    k = -1;
    do{
        k = k + 1;
    } while (!(pleaf[k] >= mdsh || pleaf[k] == 0));
    //Loop Until pleaf(k) >= mdsh Or pleaf(k) = Empty //pleaf for historical curve must match pleaf from virgin curve
    transpirationsh = eplantl[k]; //all gas exchange values are from most recent historical values
    psynactsh = psynsh[k];
    gcmdsh = gcanwsh[k]; //g for water in mmol
    lavpdmdsh = lavpdsh[k] * patm;
    cincsh = cinsh[k];
    //If k > 1 Then deda = (eplantl(k) - eplantl(k - 1)) / (psyn(k) - psyn(k - 1))
    haltsh = k; //halt is index of midday datum
    if (ecritsystem == 0){
        tenMarker:     //no midday
        k = 0;
        transpiration = eplantl[k]; //all gas exchange values are from most recent historical values
        psynact = psyn[k];
        gcmd = gcanw[k]; //g for water in mmol
        lavpdmd = lavpd[k] * patm;
        cinc = cin[k];
        halt = k;
        transpirationsh = eplantl[k]; //all gas exchange values are from most recent historical values
        psynactsh = psynsh[k];
        gcmdsh = gcanwsh[k]; //g for water in mmol
        lavpdmdsh = lavpdsh[k] * patm;
        cincsh = cinsh[k];
        haltsh = k; //halt is index of midday datum
    }
    gethistory(); //reinstates historical element curves and failure status prior to updating
    if (refilling == true){ //need to record midday kmins, uses sunlit pressures
        for (z = 1; z <= layers; z++){
            kminroot[z] = kroot[z][halt];
        }
        kminstem = kstem[halt];
        kminleaf = kleaf[halt];
    }
}

// Time-step simulation functions
/* Time step iterator function */
long MainProgram::modelTimestepIter(long& VBA_dd) {
    dd = VBA_dd;
        
    if (dd == 1 || isNewYear) {
        failure = 0;//'=1 for system failure at midday...terminates run
        failspot = "no failure";
        componentpcrits();//'gets pcrits for each component
        failspot = "no failure";

        for (k = 1; k <= layers; k++) {// k = 1 To layers ;//'exclude the top layer
            kminroot[k] = ksatr[k];
        }

        kminstem = ksats;
        kminleaf = ksatl;
        kminplant = ksatp;

        gwflow = 0; //'inflow to bottom of root zone
        drainage = 0; //'drainage from bottom of root zone

        gs_prevDay = 0;
        gs_inGrowSeason = false;
        gs_doneFirstDay = false; // prevents PLC misses on first year
        //[/HNT]
    }

    //dd++;
    long testCount = 0;
    // yearVal = std::lround(dSheet.Cells(rowD + dd, colD + dColYear));
    yearVal = std::lround(dataCells[rowD + dd][colD + dColYear]);

    if (yearVal != year_cur) 
    {// all of these year values get zeroed out during initModelVars, so we can assume they will be zero at the start of a new iteration and test against 0 meaning "not set"    
        if (yearVal > year_cur && yearVal > year_start) {
            // if the start year is zero, this is the first year we've processed
            if (year_start == 0)// model runs starting in the year zero are not supported, I guess?
                year_start = yearVal;
            year_cur = yearVal;

            gs_yearIndex = year_cur - year_start;
            if (gs_yearIndex >= 100)
                gs_yearIndex = 99; //safety check

            gs_ar_years[gs_yearIndex] = year_cur; // correct the year listing in the "growing seasons" array 
            // TODO make that input data able to handle the new system?? Otherwise just eliminate
            // on VBA side it might be sufficient to pull the start year from the first line of data
            // though this still assumes that our data is in linear time order ...

            if (gs_yearIndex > 0) {
                // if we've set the year and it was anythign other than the first timestep of this iteration (in which case we set the year index from zero TO zero)
                // then return to VBA and let it handle the model reset
                // it will re-run this timestep afterwards
                gs_prevDay = 0;
                gs_inGrowSeason = false;
                gs_doneFirstDay = false;
                isNewYear = true; // this can be used to override the dd==/>1 cases -- it will be set to false upon successful completion of a timestep
                return gs_yearIndex;
            } // if we detected a year change, but the yearIndex comes out to zero, then this must be the first timestep
            // of the first year in the data -- so we can continue without returning and re-setting the model
         } else { // going back in time means something went wrong -- may be able to handle out of order years later though 
            return 0; // failure
        }
    }

    // if ( dd == 1 || isNewYear){
    //     // Get CO2 for current year
    //     ca = getCarbonByYear(yearVal,GSCells,maxYears); // get current year atmospheric CO2
    //     std::cout << " Atmospheric CO2 concentration for " << yearVal << ": " << ca << std::endl;
    //     ca = ca * 0.000001;
    // }

    // if (dd == 0){
    //     double caDefault = std::atof(GSCells[2][4].c_str());
    //     std::cout << "WARNING: Initial CO2 for day 0: " << caDefault << std::endl;
    //     ca = caDefault * 0.000001;
    // }

    //jd = dSheet.Cells(rowD + dd, colD + dColDay); //'julian day
    jd = dataCells[rowD + dd][colD + dColDay]; // julian day

    //Debug.Print "DOING A LOOP-1 " & dd
    if (dd > 1 && !isNewYear) { //if// //'set timestep
        //if (tod < dSheet.Cells(rowD + dd, colD + dColTime)) {
        if (tod < dataCells[rowD + dd][colD + dColTime]) { //if// //'same old day getting older
            //timestep = dSheet.Cells(rowD + dd, colD + dColTime) - tod;
            timestep = dataCells[rowD + dd][colD + dColTime] - tod;
        } else {
            timestep = (24 - tod) + dataCells[rowD + dd][colD + dColTime];
            //timestep = (24 - tod) + dSheet.Cells(rowD + dd, colD + dColTime); //'a new day has started
            //[HNT] multi-year support
            // following method of new year detection has been replaced with the method above
            // now that we are including the year as a data input
            gs_prevDay = jd;
            gs_inGrowSeason = true; // isInGrowSeasonSimple(); //it's a new day, so let's see if this is in the growing season or not
            //[/HNT]
        } //End if// //'tod if
    } else {//End if// //'dd>1 if [HNT] multi-year support
        gs_inGrowSeason = true; // isInGrowSeasonSimple(); //it's the first data point, test if we're starting in growing season
    }

    gs_inGrowSeason = isInGrowSeasonSimple(); // just always call this!
    
    tod = dataCells[rowD + dd][colD + dColTime]; //'time of day, standard local time in hour fraction
    //tod = dSheet.Cells(rowD + dd, colD + dColTime); //'time of day, standard local time in hour fraction
    obssolar = dataCells[rowD + dd][colD + dColSolar]; //'observed total solar, horizontal, wm-2
    //obssolar = dSheet.Cells(rowD + dd, colD + dColSolar); //'observed total solar, horizontal, wm-2
    vpd = dataCells[rowD + dd][colD + dColD]; //'midday vpd in kPa
    //vpd = dSheet.Cells(rowD + dd, colD + dColD); //'midday vpd in kPa
    vpd = vpd / patm; //'vpd in mole fraction
    airtemp = dataCells[rowD + dd][colD + dColTAir]; //'in C
    //airtemp = dSheet.Cells(rowD + dd, colD + dColTAir); //'in C
    maxvpd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    wind = dataCells[rowD + dd][colD + dColWind]; //'wind speed
    //wind = dSheet.Cells(rowD + dd, colD + dColWind); //'wind speed
    if (wind < minwind) { //if//
        dataCells[rowD+dd][colD+dColWind] = minwind;
        //dSheet.Cells(rowD + dd, colD + dColWind) = minwind; //'set to minimum wind
        wind = minwind;
    } //End if//
    us = wind * 0.1; //'understory windspeed in m s-1
    soiltemp = dataCells[rowD + dd][colD + dColTSoil]; //'surface temp of soil
    //soiltemp = dSheet.Cells(rowD + dd, colD + dColTSoil); //'surface temp of soil
    if (vpd > maxvpd) { //if//
        vpd = maxvpd;
        dataCells[rowD + dd][colD + dColD] = maxvpd * patm; 
        // dSheet.Cells(rowD + dd, colD + dColD) = maxvpd * patm; //'print out maximum vpd
    } //End if//

    //'after initializing, start updating water contents of soil layers
    soilcalculator.get_soilwetness(drainage,runoff,waterold,x,thetafracres,a,n,thetafracfc,thetafc,depth,fieldcapfrac,thetasat,
        water,fc,ffc,dataCells,gs_ar_input,gs_ar_waterInitial_OFF,gs_ar_waterInitial,layerflow,elayer,laisl,lai,laish,baperga,
        timestep,soilredist,deficit,swclimit,tod,rain,waternew,waterchange,gs_ar_waterFinal,gwflow,transpirationtree,laperba,gs_ar_E,
        soilevap,gs_ar_ET,atree,gs_ar_Anet,cinc,ca,gs_ar_cica,gs_ar_cica_N,gs_ar_Aci,gs_ar_AnetDay,sum,kpday1,kxday1,gs_ar_waterInitial_GS,
        gs_ar_waterFinal_OFF,iter_refK,gs_ar_PLCSum,gs_ar_PLCSum_N,gs_ar_PLCp,gs_ar_PLC85,gs_ar_PLCx,gs_ar_waterInput_GS,gs_ar_waterFinal_GS,
        gs_ar_waterInput_OFF,dd,gs_yearIndex,z,layers,rowD,colD,dColF_End_watercontent,halt,haltsh,iter_Counter,stage_ID,dColRain,j,layer,
        dColF_End_waterchange,dColF_End_rain,dColF_End_gwater,dColF_End_drainage,dColF_End_input,dColF_End_runoff,o,dColF_End_E, dColF_End_soilEvap,
        dColF_End_ET,dColF_End_ANet,dColF_CP_kplant,dColF_CP_kxylem,dColF_End_PLCplant,dColF_End_PLCxylem,night,isNewYear,useGSData,
        mode_predawns,rainEnabled,ground,gs_inGrowSeason,gs_doneFirstDay); 
    //'get radiation for timestep
    photosyncalculator.get_solarcalc(fet,et,pi,sm,longitude,tsn,tsncorr,sindec,dec,cosdec,tim,tod,coszen,lat,zen,cosaz,az,m,patm,
        sp,solar,tau,sb,sd,st,cloud,obssolar,fcd,xang,kbe,kbezero,mleafang,rad,sum,k1,t1,lai,told,t2,kd,qd,qds,abspar,qdt,qb,qbt,
        qsc,qsh,qsl,laisl,laish,parsh,parsl,parbottom,absnir,absolar,nirsh,nirsl,sshade,ssun,sbottom,ssunb,ssund,sref,par,ppfd,ea,
        maxvpd,vpd,eac,airtemp,la,sbc,lg,jd);
    // night yet?
    if (qsl > lightcomp /*[HNT] multiyear*/ && gs_inGrowSeason) {
        night = "n"; //'it//'s light enough for sun layer to do business
    } else {
        night = "y"; //'too dark
    } //End if// //'end night if

    gwflow = 0; //'re-set inflow to bottom of root zone
    drainage = 0; //'re-set drainage from bottom of root zone
    //'Call getpredawns //'update soil pressure of each layer
    //'if failure = 1 { //if// Exit do
    chalk = 0; //'einc counter

    twentyMarker:

    //'passed initializing...update soil pressure of each layer
    soilcalculator.get_predawns(kminroot,theta,water,depth,thetasat,x,pd,dataCells,pgrav,a,n,prh,pcritrh,pcritr,sum,pr,prinitial,k,layers,
        layer,z,rowD,dd,colD,dColRain,t,failure,dColF_p1,o,layerfailure,failspot,mode_predawns);

    if (failure == 1){
        std::cout << "WARNING: Unrecoverable model failure!" << std::endl;
        std::cout << "SOURCE: Errors when calculating soil pressures" << std::endl;
        std::cout << "ACTION: Model stops " << std::endl;
        std::cout << std::endl;
        return -1;//break;//Exit do
    }

    test = 0; //'=1 if stem or leaf fails
    p = -1; //'E(Pleaf) point counter
    //'set initial guesses of the three unknowns
    e = -einc; //'total e: einc and e are still in kg hr-1 m-2 basal area, as are conductances
    //'Range("c18:i10000").ClearContents
    psynmax = -100;
    psynmaxsh = -100;
    skip = 0; //'this turns off psynthesis

    do{//'this loop obtains and stores the entire composite curve
        p = p + 1;
        e = e + einc;
        //'solves for p//'s and e//'s in fingers of chicken foot
        newtonrhapson();
        memset(vv, 0, sizeof(vv));//Erase vv
        //'if check >= 400 { //if// GoTo 60: //'NR can//'t find a solution
        if (check > 500){ //if//
            break;//Exit do
        } //End if// //'test for total failure
        
        sum = 0;
        for (k = 1; k <= layers; k++) {//k = 1 To layers
            sum = sum + layer[k];
        } //end for k
        
        if (sum == layers) { //if//
            failspot = "below ground"; //'total failure below ground
            break;//Exit do
        } //End if//
        
        get_stem(); //'gets stem and leaf pressures
        get_leaf();
        
        if (test == 1) {
            break;//Exit do
        } 
        
        compositecurve(); //'stores the entire composite curve //if skip = 0 { //if//
        
        // gets sun layer leaf temperature from energy balance
        photosyncalculator.get_leaftemps(rabs,ssun,sref,emiss,la,lg,lambda,airtemp,grad,gha,wind,
            leafwidth,laperba,eplantl,eplant,numerator,denominator,sbc,sha,leaftemp,lavpd,
            patm,vpd,p);
        // gets shade layer leaf temperature
        photosyncalculator.get_leaftempsshade(rabs,sshade,sref,emiss,lg,numerator,sbc,
            airtemp,lambda,eplantl,denominator,sha,grad,gha,leaftempsh,lavpdsh,patm,
            vpd,p);
        // gets sun layer photosynthesis
        photosyncalculator.get_assimilation(lavpd,gcanw,gmax,eplantl,gcanc,comp25,comp,gas,leaftemp,
            numerator,svvmax,hdvmax,havmax,denominator,vmax25,vmax,svjmax,hdjmax,hajmax,jmax,jmax25,
            kc25,ko25,kc,ko,rday25,rday,ci,jact,qmax,qsl,lightcurv,je,oa,jc,var,thetac,ca,psyn,cin,
            marker,psynmax,p,night);
        // gets shade layer photosynthesis
        photosyncalculator.get_assimilationshade(lavpdsh,gcanwsh,gmax,eplantl,gcancsh,comp,comp25,gas,
            leaftempsh,numerator,svvmax,hdvmax,havmax,denominator,vmax25,vmax,svjmax,hdjmax,hajmax,jmax,
            jmax25,kc25,ko25,kc,ko,rday25,rdaysh,ci,jact,qmax,qsh,lightcurv,je,oa,jc,var,thetac,ca,psynsh,
            cinsh,marker,psynmaxsh,p,night);
    } while (!(sum == layers || test == 1 || night == "y" && (dd > 1 && !isNewYear) || check >= 500)); //'loop to complete failure unless it//'s night

    if (chalk > 0) { //if//
        weird = 0; //'done our best
        failspot = "convergence";
        for (z = 1; z <= layers; z++) {//z = 1 To layers //'restore layers to functioning if they//'ve been turned off by convergence failure
            if (kminroot[z] != 0) {
               layer[z] = 0;
            }
        } //end for z
        goto fortyMarker; //'got as much of the composite curve as is going to happen
    } //End if//
    
    if (dd == 1 || isNewYear || night == "n") { //if//
        if (check >= 500) { //if// //'try once more //'Stop
            chalk = chalk + 1;
            if (ecritsystem == 0) {
                einc = ksatp / 500.0;
                std::cout << "ecritsystem is zero... try resetting to ksatp/500, dd = " << dd << std::endl;
            }
            goto twentyMarker;
        } //End if//

        if (total > 500 || total < 400) { //if//
            einc = ecritsystem / 450.0; //'re-set Einc
            if (ecritsystem == 0) {
                einc = ksatp / 500.0;
                std::cout << "ecritsystem is zero... try resetting to ksatp/500, dd = " << dd << std::endl;
            }
            testCount++; // [DEBUG]
            if (testCount > 10) {
                testCount = 0;
                goto fortyMarker;
            }
            goto twentyMarker; //'recalculate the composite curve
        } //End if// //'total ok
    } //End if// //'night <>"n" or it//'s not the first round

    fortyMarker:

    bool isNight = true;
    if (night == "n"){
        isNight = false;
    }
    if (night == "n" && psynmax > 0 && psynmaxsh > 0 && weird == 0) { //if//
    //DoEvents //'[HNT] this was required to prevent a hard lock -- this portion of the loop is the most intensive, so let Excel take a "breath" by processing system events to prevent lockup
        hydraulicscalculator.get_canopypressure(ecritsystem,ter,er,kr,tkr,er_v,kr_v,tes,es,es_v,tel,el,el_v,
            kminroot,sum,prh,pd,pcritr,pr,psynmaxmd,psynmaxshmd,e,einc,dedplmin,ksatp,jmatrix,frt,dfrdpr,
            pcritrh,p1,p2,plow,pinc,elow,erh,klow,krh,ehigh,khigh,estart,klower,efinish,kupper,flow,func,
            dfrhdprh,dfrhdpr,dfrdprh,threshold,initialthreshold,aamax,vv,dum,indx,pcrits,waterold,
            gs_ar_nrFailConverge_Water,gs_ar_nrFailConverge_WaterMax,pgrav,ps,pl,pcritl,pleafv,predawn,
            plold,dedplzero,dedpl,klossv,pcritsystem,rabs,ssun,sref,emiss,la,lg,lambda,airtemp,grad,gha,wind, 
            leafwidth,emd,laperba,sbc,numerator,denominator,sha,leaftmd,lavpdmd,patm,vpd,sshade,leaftshmd,
            lavpdshmd,gcanwmd,gmax,gcancmd,comp,comp25,gas,svvmax,havmax,hdvmax,vmax25,vmax,svjmax,hdjmax,hajmax,
            jmax,jmax25,kc25,ko25,kc,ko,rday25,rdaymd,ci,jact,qmax,qsl,lightcurv,je,oa,jc,var,thetac,ca,psynmd,
            cinmd,marker,gcanwshmd,gcancshmd,rdayshmd,qsh,psynshmd,cinshmd,maxkloss,dpmax,dpamax,amaxmax,dpamin,
            amaxfrac,dpa,rmean,md,lavpd,dpasun,mdsh,amaxfracsh,lavpdsh,pleaf,eplantl,transpiration,psynact,psyn,
            gcmd,gcanw,cinc,cin,transpirationsh,psynactsh,psynsh,gcmdsh,gcanwsh,lavpdmdsh,cincsh,cinsh,kminstem,
            kminleaf,kroot,kstem,kleaf,k,layer,tlayer,t,failure,test,p,gs_ar_nrFailConverge,weird,check,ticks,
            unknowns,d,imax,ii,ll,gs_yearIndex,dd,totalv,runmean,total,cutoff,halt,haltsh,layers,layerfailure,
            tlayerfailure,failspot,night,refilling); //'returns canopy P and associated output
                           //'if check >= 2000 { //if// GoTo 60:
        if (refilling == true) {
            hydraulicscalculator.update_curves(kroot,kminroot,pinc,proot,er,kr,kstem,kminstem,pstem,es,kleaf,
                kminleaf,pleaf,el,halt,phigh,layers); //'updates element E(P) curves as required for midday exposure for no refilling
        }
    } //End if// //'night <> "n", psynmax IF

    if (soilred == true) { //if//
       soilcalculator.get_soilflow(store,kmaxrh,kkmax,depth,pd,soilredist,p1,pend,e,p2,pinc,olds,a,n,x,s,del,sum,eps,
        soilf,z,layers,tmax,t,it,tnm,j); //'gets vertical soil flow between layers in m3/m2
    } else {
        for (z = 0; z <= layers; z++) {//z = 0 To layers
            soilredist[z] = 0;
        } //end for z
    } //End if// //'soil red <> y
    
    if (ground == true) {
        soilcalculator.get_deepflow(store,kmaxrh,kkmax,grounddistance,pground,pd,p1,pend,e,p2,pinc,olds,a,n,x,s,
            del,sum,eps,groundflow,gwflow,timestep,drainage,soilredist,z,layers,tmax,t,it,tnm,j);//'gets groundwater flow into bottom layer in m3/m2
    } //End if// //'pet <> y or n
    
    if (gs_inGrowSeason && sevap == true) { //if//
            soilcalculator.get_soilevaporation(emission,emiss,sbc,soiltemp,lc,leaftempsh,rabssoil,soilabssol,sbottom,mdensair,
                patm,airtemp,gha,us,xh,zdispl,rough,zh,soilep,sha,grad,lambda,rha,vpd,maxvpd,rhs,pd,gas,soilevap,baperga,soilredist,haltsh); //'gets soil evaporation rate
    } else {
        soilevap = 0;
    } //End if//

    if (failure == 0 || weird == 1) { //if// //Debug.Print "DOING A LOOP-8 " & dd
        if (night == "y" || psynmax == 0 || psynmaxsh == 0) { //if// //'set everything to starting point
            k = 0;
            transpiration = eplantl[k]; //'all gas exchange values are for closed stomata
            md = pleaf[k];
            psynact = psyn[k];
            gcmd = gcanw[k]; //'g for water in mmol
            lavpdmd = lavpd[k] * patm;
            cinc = cin[k];
            halt = k;
            transpirationsh = eplantl[k]; //'all gas exchange values are from most recent historical values
            mdsh = pleaf[k];
            psynactsh = psynsh[k];
            gcmdsh = gcanwsh[k]; //'g for water in mmol
            lavpdmdsh = lavpdsh[k] * patm;
            cincsh = cinsh[k];
            haltsh = k; //'halt is index of midday datum
        } //End if// //'night<>y

        // Writing time-step outputs
        dataCells[rowD + dd][colD + o + dColF_predawn] = pleaf[0]; // predawn water potential (MPa)
    
        // SUN LAYER OUTPUT
        dataCells[rowD + dd][colD + o + dColF_P] = md; // middady water potential (MPa)
        dataCells[rowD + dd][colD + o + dColF_E] = transpiration; // midday transpiration (mmol s-1 m-2 leaf area)
        dataCells[rowD + dd][colD + o + dColF_Gw] = gcmd; // midday canopy diffusive conductance to water (mmol s-1m-2)
        dataCells[rowD + dd][colD + o + dColF_laVPD] = lavpdmd; // leaf to air vpd (kPa)
        dataCells[rowD + dd][colD + o + dColF_leaftemp] = leaftemp[halt]; // leaf temperature (Celsius)
        dataCells[rowD + dd][colD + o + dColF_ANet] = psynact; // net carbon assimilation by leaf area (umol s-1m-2)
        //'anet = psynact //'keeping this old setting just in case downstream Anet is needed
        //'Cells(16 + dd, o + 21) = chalk //'number of NR restarts
        dataCells[rowD + dd][colD + o + dColF_ci] = cinc * patm * 1000; // partial pressure of CO2 (Pa)
        dataCells[rowD + dd][colD + o + dColF_PPFD] = qsl; // photon flux density (umol s-1m-2)
    
        // SHADE LAYER OUTPUT
        dataCells[rowD + dd][colD + o + dColF_S_P] = mdsh; // midday water potential (MPa)
        dataCells[rowD + dd][colD + o + dColF_S_E] = transpirationsh; // midday transpiration (mmol s-1 m-2 leaf area)
        dataCells[rowD + dd][colD + o + dColF_S_Gw] = gcmdsh; // midday canopy diffusive conductance to water (mmol s-1m-2)
        dataCells[rowD + dd][colD + o + dColF_S_laVPD] = lavpdmdsh; // leaf to air vpd (kPa)
        dataCells[rowD + dd][colD + o + dColF_S_leaftemp] = leaftempsh[haltsh]; // leaf temperature (Celsius)
        dataCells[rowD + dd][colD + o + dColF_S_Anet] = psynactsh; // net carbon assimilation by leaf area (umol s-1m-2)
        //'anetsh = psynactsh //'see above
        //'Cells(16 + dd, o + 30) = dpasun //'for debugging purposes
        dataCells[rowD + dd][colD + o + dColF_S_ci] = cincsh * patm * 1000; // partial pressure of CO2 in Pa
        dataCells[rowD + dd][colD + o + dColF_S_PPFD] = qsh; // photon flux density (umol s-1m-2)
    
        //'WHOLE TREE OUTPUT
        if (night == "n"){
            transpirationtree = laisl / lai * transpiration + laish / lai * transpirationsh; //'weighted mean
        }    
        if (night == "y"){
            transpirationtree = transpiration;
        }
        if (night == "n"){
            atree = laisl / lai * psynact + laish / lai * psynactsh; //'weighted mean
        }
        if (night == "y"){
            atree = (psynact + psynactsh) / 2.0; //'simple average at night when there//'s no sun or shade leaves
        }
        dataCells[rowD + dd][colD + o + dColF_T_E] = transpirationtree; // Whole tree transpiration
        dataCells[rowD + dd][colD + o + dColF_T_ANet] = atree; // Whole tree Anet
        //'Cells(16 + dd, o + 35) = dpamax //'shade leaf dpa
    
        //'HYDRAULIC OUTPUT (BASED ON SUN MD)
        dataCells[rowD + dd][colD + o + dColF_T_pcrit] = pcritsystem;
        dataCells[rowD + dd][colD + o + dColF_T_Ecrit] = ecritsystem * (1 / laperba) * (1.0 / 3600.0) * 55.4 * 1000; // ecrit (mmol s-1m-2)
        dataCells[rowD + dd][colD + o + dColF_CP_Pstem] = pstem[halt];
        dataCells[rowD + dd][colD + o + dColF_CP_Proot] = proot[halt];
        dataCells[rowD + dd][colD + o + dColF_CP_kstem] = kstem[halt]; // k stem at midday (kg hr-1m-2MPa-1)
        dataCells[rowD + dd][colD + o + dColF_CP_kleaf] = kleaf[halt]; // k in leaf at midday

        if (transpiration > 0) { //if//
            kplantold = eplant[halt] / (pleaf[halt] - pleaf[0]);  // whole plant k at midday (kg hr-1 m-2 basal area...sun value)
            dataCells[rowD + dd][colD + o + dColF_CP_kplant] = kplantold;
        } //End if//
        if (transpiration == 0){
            dataCells[rowD + dd][colD + o + dColF_CP_kplant] = kplantold; // use most recent kplant
        }
    
        if (kplantold < gs_ar_kPlant[gs_yearIndex] || gs_ar_kPlant[gs_yearIndex] == 0){
            gs_ar_kPlant[gs_yearIndex] = kplantold;
        }
        k = o + dColF_CP_kroot1 - 1;//43;
        sum = 0;
        for (z = 1; z <= layers; z++) {//z = 1 To layers
            dataCells[rowD + dd][colD + k + z] = kroot[z][halt]; // root k at midday, sun fluxes
            sum = sum + kroot[z][halt];
        } //end for z
        dataCells[rowD + dd][colD + k + 1 + layers] = sum; // total root k at midday
        if (failure == 0) { //if//
            double tempDouble = 0.0;
            tempDouble = 1 / (1 / kminleaf + 1 / kminstem + 1 / dataCells[rowD + dd][colD + k + 1 + layers]); // total xylem k
            dataCells[rowD + dd][colD + k] = tempDouble;
            if (tempDouble < gs_ar_kXylem[gs_yearIndex] || gs_ar_kXylem[gs_yearIndex] == 0)
                gs_ar_kXylem[gs_yearIndex] = tempDouble;
        } //End if//
        
        for (z = 1; z <= layers; z++){//z = 1 To layers
            dataCells[rowD + dd][colD + k + 1 + layers + z] = elayer[z][halt] * (1 / laperba) * (1.0 / 3600.0) * 55.56 * 1000; //'uptake in mmol s-1m-2 leaf area...sun rate
        } //end for z
        //'Cells(16 + dd, k + 2 + 2 * layers) = failspot //'position of failure at critical point

        // TODO since all IO is being handled as double currently, cannot add these failure notes
        // temporarily putting in an obvious number to flag failure
        for (z = 1; z <= 1; z++){//z = 1 To 1
            if (layer[z] == 1){
                dataCells[rowD + dd][colD + k + 2 + 2 * layers + z] = -1137;// layerfailure[z]; //'layers failed at critical point
            }
        } //end for z
        //Debug.Print "DOING A LOOP-9 " & dd
    } //End if// //'failure IF (basically...failure can//'t happen!)


    if (dd == 1 || isNewYear) { //if// //'NOTE: must be sure that pcritsystem is computed for dd=1!!! (i.e., it//'s not computed at night)
        x = pcritsystem; //'estimate of "permanent wilting point"
        for (z = 1; z <= layers; z++) { //z = 1 To layers
            swclimit[z] = swc(x); //'theta/thetasat at critical point
            swclimit[z] = swclimit[z] * thetasat[z]; //'convert to water content
            swclimit[z] = swclimit[z] * depth[z]; //'water content left over in m3/m2 ground area
            //'sumsoil = sumsoil + (fc[z] - swclimit[z]) //'sum is total m3 water per m2 ground withdrawn
        } //end for z
        //Debug.Print "DOING A LOOP-11 " & dd
    } //End if//

    //'now...need to reset layer failure status to midday values (not critical point)
    for (z = 1; z <= layers; z++) {//z = 1 To layers
        if (layer[z] == 1) { //if// //'check to see if kminroot[z]=0; otherwise, restore it to function
            if (layerfailure[z] == "root" && kminroot[z] > 0) { //if//
                layer[z] = 0;
                layerfailure[z] = "ok";
            } //End if//
            if (layerfailure[z] == "rhizosphere" && kminroot[z] > 0) { //if//
                layer[z] = 0;
                layerfailure[z] = "ok";
            } //End if//
        } //End if//
    } //end for z

    if (isNewYear){
        isNewYear = false; // always set this
    }
    return -1;
}

/* New year parameterization*/
void MainProgram::modelProgramNewYear(){
    // save the iterate-able water system states
    bool oldGround;
    bool oldRaining;
    bool oldRainEnabled;
    oldGround = ground;
    oldRaining = raining;
    oldRainEnabled = rainEnabled;

    // clear all the model variables
    cleanModelVars(); //initialize all used variables to avoid any memory effect from previous runs
    InitialConditions(); // get all the parameters again

    gs_prevDay = 0;
    gs_inGrowSeason = false;

    // set up the iteration details based on what we saved
    if (iter_runSupplyCurve) {
        ground = false; //disable ground water for the first two iterations
        raining = false; //disable rain for the first iteration
        rainEnabled = false; //disable rain for the first iteration
        ground = oldGround;
        raining = oldRaining;
        rainEnabled = oldRainEnabled;

        if (iter_bagaEnable) {
            baperga = iter_baga * 0.0001;
            //if we set the BA:GA { also set the LAI
            lai = (laperba * baperga) / treeToPhotoLAI;
            std::cout << "DEBUG: set BA:GA to " << iter_baga << " m2/ha -> " << baperga << "m2/m2. treeLAI/fotoLAI = " << treeToPhotoLAI << " and LAI set to " << lai << std::endl;
        } //End if
    } //End if

    for (k = 0; k <= layers; k++){// To layers //assign source pressures, set layer participation  
        layerfailure[k] = "ok";
        layer[k] = 0; //1 if out of function
    } //Next k

    //[HNT] all this is done in C now
    if (true) { //Not useDLL {
        failure = 0; //=1 for system failure at midday...terminates run
        failspot = "no failure";
        componentpcrits(); //gets pcrits for each component
        failspot = "no failure";

        for (k = 1; k <= layers; k++) {//k = 1 To layers //exclude the top layer
            kminroot[k] = ksatr[k];
        } //Next k

        kminstem = ksats;
        kminleaf = ksatl;
        kminplant = ksatp;

        gwflow = 0; //inflow to bottom of root zone
        drainage = 0; //drainage from bottom of root zone
    } //End if
    //if kmaxset = False Or leafpercentset = False { dd = 0 //start with initializing row
    //if kmaxset = True And leafpercentset = True { dd = 1 //skip initializing row
    //dd = 0
    //liveGraphNextUpdate = liveGraphUpdate
    //cleanModelVars(); //initialize all used variables to avoid any memory effect from previous runs
    //readin();
    //Call loadParamsToC
    //Call CPP_setIterationCount(iter_Counter)
}

/* Save output */
void MainProgram::saveOutputSheet(std::string filename, std::string sheetType){
    std::cout << "Called saveOutputSheet! Filename " << filename << " type " << sheetType << std::endl;
    // output!
    std::ofstream dataOut;
    dataOut.open(filename + ".csv", std::ios::out);
    if (!dataOut.is_open()){
        std::cout << "FAILED to open file " << filename + ".csv" << std::endl;
    }

    dataOut.precision(FIO_PRECISION);

    //int rowCount = -1;
    if (sheetType == "timesteps"){
        long maxRow = (gs_yearIndex + 1) * 366 * 24 + 10;
        //if (useGSData)
        //   maxRow = 310000; // 306,600 hours in 35 years (307,440 in 35 leap years) ... leave a buffer

        for (int rowCount = 1; rowCount < maxRow; rowCount++){// excel can only load ~ 1 mil rows
            if (rowCount < rowD){
                for (int hRC = 0; hRC < dataHeaderRow.size(); hRC++){
                    dataOut << dataHeaderRow[hRC] << ",";
                }
            } else if (rowCount >= rowD){
                for (int rC = 1; rC < 80; rC++){
                    dataOut << dataCells[rowCount][rC] << ",";
                }
            }
            dataOut << std::endl; // no matter what, write a new line at the end of the row output
        }
        std::cout << "Write complete!" << std::endl;
    } else if (sheetType == "summary"){
        for (int rowCount = 1; rowCount < MAX_SUMMARY_ROWS; rowCount++){// excel can only load ~ 1 mil rows
            if (rowCount == rowD){
                for (int hRC = 0; hRC < summaryHeaderRow.size(); hRC++){
                  dataOut << summaryHeaderRow[hRC] << ",";
                }
            } else if (rowCount > rowD){
                if (finalOutCells[rowCount][dColF_GS_year] > 0.001){
                    for (int rC = 1; rC < MAX_SUMMARY_COLS; rC++){
                        dataOut << finalOutCells[rowCount][rC] << ",";
                    }
                } else{
                    break;
                }
            }
            dataOut << std::endl; // no matter what, write a new line at the end of the row output
        }
        std::cout << "Write complete!" << std::endl;
    } else if (sheetType == "bagahist"){
        // write this one to the root!
        // species/runtype/site/scen/model
        dataOut << stage_OptHistBAGA << "," << std::endl;
        std::cout << "Write complete!" << std::endl;
    } else if (sheetType == "bagafut"){
        // write this one to the root!
        // species/runtype/site/scen/model
        dataOut << stage_OptFutBAGA << "," << std::endl;
        std::cout << "Write complete!" << std::endl;
    }
}

/* Running the model */ 
long MainProgram::Rungarisom(){
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << "|  CARBON GAIN VS HYDRAULIC RISK MODEL V 2.0  |" << std::endl;
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << std::endl;

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
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << "|        READING ClIMATE FORCING FILES        |" << std::endl;
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << "- Growing Season Data: " << std::endl;
    MainProgram::readGSSheet();
    MainProgram::readGrowSeasonData();
    std::cout << std::endl;  
    std::cout << "- Climate Forcing Data: " << std::endl;
    MainProgram::readDataSheet();
    std::cout << std::endl;
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << "|             INITIALIZING MODEL              |" << std::endl;
    std::cout << " ---------------------------------------------" << std::endl;
    std::cout << std::endl;

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
    componentpcrits(); //gets pcrits for each component
    failspot = "no failure";
    
    for (k = 1; k <= layers; k++){ // k = 1 To layers //exclude the top layer
        kminroot[k] = ksatr[k];
    } // Next k

    kminstem = ksats;
    kminleaf = ksatl;
    kminplant = ksatp;

    gwflow = 0; //inflow to bottom of root zone
    drainage = 0; //drainage from bottom of root zone

    dd = 0;
    long ddMod = 0;
    long successCode = 0;
    std::cout << "------------------ TESTING AREA ------------------"<< std::endl; 
    std::cout << std::endl;
    std::cout << "Patm: " << patm << std::endl;
    std::cout << "No Soil layers: " << layers << std::endl;
    std::cout << "Root beta: " << beta << std::endl;
    std::cout << "layer: " << 1 <<std::endl;
    std::cout << "alpha: " << a[1] <<std::endl;
    std::cout << " VG n: " << n[1] <<std::endl;
    std::cout << " kkmax: " << kkmax[1] <<std::endl;
    std::cout << "thetas: " << thetasat[1] <<std::endl;
    std::cout << std::endl;
    std::cout << "layer: " << 3 <<std::endl;
    std::cout << "alpha: " << a[3] <<std::endl;
    std::cout << " VG n: " << n[3] <<std::endl;
    std::cout << " kkmax: " << kkmax[3] <<std::endl;
    std::cout << "thetas: " << thetasat[3] <<std::endl;
    std::cout << std::endl;
    std::cout << "------------------ TESTING table ------------------"<< std::endl; 
    std::cout << std::endl;
    std::cout << "layer: " << soillayersTable[rowLR+1][colLR+1] <<std::endl;
    std::cout << "alpha: " << soillayersTable[rowLR+1][colLR+3]  <<std::endl;
    std::cout << " VG n: " << soillayersTable[rowLR+1][colLR+4]  <<std::endl;
    std::cout << " kkmax: " << soillayersTable[rowLR+1][colLR+5]  <<std::endl;
    std::cout << "thetas: " << soillayersTable[rowLR+1][colLR+6]  <<std::endl;
    std::cout << "depth: " << soillayersTable[rowLR+1][colLR+11]  <<std::endl;
    std::cout << std::endl;

    do { // loop through time steps
        dd = dd + 1;
        successCode = modelTimestepIter(dd);

        if(successCode == 0){
            std::cout << "Unrecoverable model failure!" << std::endl;
            std::cout << "Model stops " << std::endl;
            std::cout << std::endl;
            return 0; // failure, unrecoverable
        } else if (successCode > 0) { // this returns the year if we've incremented it -- not necessary in the full C version (also only supports 1 year right now)
            dd = dd - 1; //we need to repeat this timestep because we bailed early when finding a new year
            gs_yearIndex = successCode;

            /* if we're running without growing season limits, we need to record the "end of GS" water content now
            because we did not complete the previous timestep, back up 1 more to grab a value */
            if(!useGSData && gs_yearIndex > 0){
                gs_ar_waterFinal_GS[gs_yearIndex - 1] = dataCells[rowD + dd - 1][colD + dColF_End_watercontent]; // make sure this goes with the previous year
            }
            std::cout << "------------------ TESTING AREA ------------------"<< std::endl; 
            std::cout << std::endl;
            modelProgramNewYear();
        } else { // -1 = success, VBA bool convention
            int breakpoint = 1137; // success, in the C version we just continue instead of outputting
            // do all CSV writing at the end
            if(!useGSData && gs_yearIndex > 0){
                gs_ar_waterFinal_GS[gs_yearIndex] = dataCells[rowD + dd - 1][colD + dColF_End_watercontent];
            } // if this was the end of the set of years, gs_yearIndex will not have been changed so use as-is
        }
        
        if(dd % 1000 == 0){
            std::cout << "Timestep " << dd << " completed" << std::endl;
        }
    } while (!(dataCells[rowD + 1 + dd][colD + dColDay] < 0.01));
    //Dim gsCount As Long
    long gsCount = 0;
    for (gsCount = 0; gsCount <= gs_yearIndex; gsCount++){// gsCount = 0 To gs_yearIndex
        if (gs_ar_years[gsCount] > 0) { //don//t bother with years that don//t exist
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_year] = gs_ar_years[gsCount]; //gs_ar_Anet(gs_yearIndex) //[HNT] todo improve I don//t like using this hard-coded constant for the size of the C interface array here, maybe check how many data columns there really are in the sheet
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_input] = gs_ar_input[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_Anet] = gs_ar_Anet[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_E] = gs_ar_E[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_PLCp] = gs_ar_PLCp[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_PLCx] = gs_ar_PLCx[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_kPlant] = gs_ar_kPlant[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_kXylem] = gs_ar_kXylem[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET] = gs_ar_ET[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 1] = rainEnabled;
            if(ground == true){
                finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 2] = 1.0;
            } else {
                finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 2] = 0.0;
            }
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 3] = ffc;
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 4] = grounddistance;
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 5] = gs_ar_PLC85[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 6] = baperga / 0.0001;
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 7] = laperba; // need this for the final calcs
            //dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 7) = gs_ar_kPlantMean[gsCount] / gs_ar_kPlantMean_N[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 8] = gs_ar_PLCSum[gsCount] / gs_ar_PLCSum_N[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 9] = gs_ar_waterFinal[gsCount] - gs_ar_waterInitial[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 10] = gs_ar_waterInitial[gsCount]; // we want to know what the initial was too, in case it's not FC
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 11] = gs_ar_waterInitial_OFF[gsCount]; // initial content for preceding off-season
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 12] = gs_ar_waterInput_OFF[gsCount]; // input for preceding off-season
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 13] = gs_ar_waterFinal_OFF[gsCount]; // final content for preceding off-season
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 14] = gs_ar_waterInitial_GS[gsCount]; // initial content for growing season
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 15] = gs_ar_waterInput_GS[gsCount]; // input for growing season
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 16] = gs_ar_waterFinal_GS[gsCount]; // final content for growing season
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 17] = gs_ar_nrFailConverge[gsCount]; // number of convergence failures
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 18] = gs_ar_nrFailConverge_Water[gsCount] / gs_ar_nrFailConverge[gsCount]; // avg water content during convergence failure
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 19] = gs_ar_nrFailConverge_WaterMax[gsCount]; // MAX water content during convergence failure
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 20] = gs_ar_nrFailThreshold[gsCount]; // convergence failure threshold
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 21] = gs_ar_cica[gsCount] / gs_ar_cica_N[gsCount];
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 23] = (gs_ar_Aci[gsCount] / gs_ar_AnetDay[gsCount]) * patm * 1000.0;
            finalOutCells[rowD + iter_Counter + 1 + gsCount /* *40 */][colD /*+ gsCount * 16*/ + dColF_GS_ET + 24] = (gs_ar_Aci[gsCount] / gs_ar_AnetDay[gsCount]) / ca;
            /*memcpy(&testProg, this, sizeof(ModelProgram));
            std::cout << "TEST! My REAL ci/ca " << gs_ar_cica[gsCount] / gs_ar_cica_N[gsCount] << std::endl;
            std::cout << "TEST! My COPIED ci/ca " << testProg.gs_ar_cica[gsCount] / testProg.gs_ar_cica_N[gsCount] << std::endl;*/
        } else{
            break; //quit the loop if we reach the end of the growing seasons list early somehow
        } // End If
    }// Next gsCount
    saveOutputSheet("./" + stageNames[stage_ID] + "_OUTPUT_timesteps", "timesteps");
    saveOutputSheet("./" + stageNames[stage_ID] + "_OUTPUT_summary", "summary");
    return 1;
}