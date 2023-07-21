// Functions to calculate soil parameters at start
#include "02Soils.h"


/* Soil structure */
// Layer depths and % root ksat
double soils::get_depths(long &k,double &layers, double &beta){
    //Cells(8 + k, 2) = 100 * (0.995 / layers) //equal % roots per layer
    //Cells(8 + k, 11) = 0.01 * Log(1 - k * 0.995 / layers) / Log(beta) //lower depth of each layer converted to m
    return 0.01 * log(1.0 - k * 0.995 / layers) / log(beta);//lower depth of each layer converted to m
    //paramCells[rowLR + k][colLR + 11] = std::to_string(layerDepths[k]);
}

// Maximum depth
//depthmax = Cells(8 + layers, 11) //max depth in meters
//depthmax = lSheet.Cells(rowLR + layers, colLR + 11) //max depth in meters

// Get half depths
double soils::get_halfdepths(long &k, double &layers, double &beta) {
    return 0.01 * log(1 - k * 0.995 / (layers * 2.0)) / log(beta);
}

// Layer volume
double soils::get_layvolume(double &depth, double &pi, double &depthmax, double &aspect){
    return depth * pi * pow((depthmax * aspect), 2.0);
}

// Layer radial width
double soils::get_radius(double &vol, double &depth,double &pi) {
    return pow((vol / (depth * pi)), 0.5);
}

/* Rhizosphere resistance
    Solve for what rhizor has to be at the target
*/

double soils::get_rhizor(double &rplant, double &rhizotarg)
{
    return rplant * (rhizotarg / (1.0 - rhizotarg));
}

// Max resistance in the rhizosphere
double soils::get_kmaxrh(double &rhizor, double &vgterm){
   return (1.0 / rhizor) / vgterm;
}

/*Get Van Genuchten alpha // override if provided*/
// This function obtains the Van Genuchten parameters for a given texture
void soils::get_vgparams(std::string &texture, long &layers, std::string (&soillayersTable)[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS], long &rowLR, long &colLR) {
    std::string a, n,soilkmax, thetasat;
    if (texture == "sand"){
        a = "1479.5945"; 
        n = "2.68";
        soilkmax = "30305.88"; 
        thetasat = "0.43";
    } else if (texture == "loamy sand"){
        a = "1265.3084";
        n = "2.28";
        soilkmax = "14897.84";
        thetasat = "0.41";
    } else if (texture == "sandy loam") {
        a = "765.3075";
        n = "1.89";
        soilkmax = "4510.168";
        thetasat = "0.41";
    } else if (texture == "loam") {
        a = "367.3476";
        n = "1.56";
        soilkmax = "1061.216";
        thetasat = "0.43";
    } else if (texture == "silt") {
        a = "163.2656";
        n = "1.37";
        soilkmax = "255.1";
        thetasat = "0.46";
    } else if (texture == "silt loam") {
        a = "204.082";
        n = "1.41";
        soilkmax = "459.18";
        thetasat = "0.45";
    } else if (texture == "sandy clay loam") {
        a = "602.0419";
        n = "1.48";
        soilkmax = "1336.724";
        thetasat = "0.39";
    } else if (texture == "clay loam") {
        a = "193.8779";
        n = "1.31";
        soilkmax = "265.304";
        thetasat = "0.41";
    } else if (texture == "silty clay loam") {
        a = "102.041";
        n = "1.23";
        soilkmax = "71.428";
        thetasat = "0.43";
    } else if (texture == "sandy clay") {
        a = "275.5107";
        n = "1.23";
        soilkmax = "122.448";
        thetasat = "0.38";
    } else if (texture == "silty clay") {
        a = "51.0205";
        n = "1.09";
        soilkmax = "20.408";
        thetasat = "0.36";
    } else if (texture == "clay") {
        a = "81.6328";
        n = "1.09";
        soilkmax = "204.08";
        thetasat = "0.38";
    } else {
        std::cout << "WARNING: Unrecoverable model failure!" << std::endl;
        std::cout << "SOURCE: Incorrect soil texture category" << std::endl;
        std::cout << "ACTION: Model stops " << std::endl;
        std::cout << std::endl;
        abort();
    }
    for(long k=1;k<=layers;k++){
        soillayersTable[rowLR + k][colLR + 3] = a; 
        soillayersTable[rowLR + k][colLR + 4] = n;
        soillayersTable[rowLR + k][colLR + 5] = soilkmax;
        soillayersTable[rowLR + k][colLR + 6] = thetasat;
    }
}

// EQ1
// vp = 1 / ((a(z) * x) ^ n(z) + 1) 
double soils::get_vgp(double (&a)[6], double (&n)[6], double &x, long &z) {
    return 1.0 / (pow((a[z] * x), n[z]) + 1);
}

// EQ2
// vp ^ ((n[z] - 1) / (2 * n[z])) * ((1 - vp) ^ ((n[z] - 1) / n[z]) - 1) ^ 2
double soils::get_vgterm(double (&n)[6], double &vp, long &z)
{
    return pow(vp, ((n[z] - 1) / (2.0 * n[z]))) * pow((pow((1 - vp), ((n[z] - 1) / n[z])) - 1), 2.0);
}

// EQ3 // does function does everything than the other 2, might keep just 1
double soils::get_vg(double (&a)[6], double (&n)[6], double &x, double (&kmaxrh)[6], long &z) //'the van genuchten function for soil k
{
    //vp = 1 / ((a(z) * x) ^ n(z) + 1)
    //vg = kmaxrh(z) * vp ^ ((n(z) - 1) / (2 * n(z))) * ((1 - vp) ^ ((n(z) - 1) / n(z)) - 1) ^ 2
    double vp;
    vp = 1 / (pow((a[z] * x), n[z]) + 1);
    return kmaxrh[z] * pow(vp, ((n[z] - 1) / (2 * n[z]))) * pow((pow((1 - vp), ((n[z] - 1) / n[z])) - 1), 2.0);
}

/* Van Genuchten Integration */
// integrates van genuchten
void soils::trapzdvg(double (&a)[6], double (&n)[6], double &x, double (&kmaxrh)[6], double &p1, double &p2, double &s,double& del, double& sum,
    long &t, long &z, long &it, long &tnm, long &j) {
    if (t == 1) {
        double vg_p1 = get_vg(a,n,p1,kmaxrh,z);
        double vg_p2 = get_vg(a,n,p2,kmaxrh,z);
        s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
        it = 1;
    } else {
        tnm = it;
        del = (p2 - p1) / tnm;
        x = p1 + 0.5 * del;
        sum = 0;
        
        for (j = 1; j <= it; j++) {
            sum = sum + get_vg(a, n, x, kmaxrh,z);
            x = x + del;
        }
        s = 0.5 * (s + (p2 - p1) * sum / tnm);
        it = 2 * it;
    }
}

// evaluates accuracy of van genuchten integration
void soils::qtrapvg(double &olds, double (&a)[6], double (&n)[6], double &x, double (&kmaxrh)[6], double &p1, double &p2,
    double &s,double &del, double &sum, double &eps, long &tmax, long &t, long &z, long &it, long &tnm, long &j) {
    
    eps = 0.001; //fractional refinement threshold for integral
    tmax = 70; //limit of calls to trapzdvg
    olds = -1; //starting point unlikely to satisfy if statement below
    
    for (t = 1; t <= tmax; t++){
        trapzdvg(a,n,x,kmaxrh,p1,p2,s,del,sum,t,z,it,tnm,j);
        if (std::abs(s - olds) < eps * std::abs(olds)){
            return;
        }
        olds = s;
    }
}

//Generate soil E(P) global curve--only need to do this once
void soils::get_rhizoPcrit(long &z, long &layers, double &p1, double &p2,long &k, double &e,double (&erh)[6][100001], double (&krh)[6][100001],
    double &pinc, double &kmin, double &olds, long &t,double (&a)[6], double (&n)[6], double (&kmaxrh)[6], double &s,
    long &it, long &tnm, double &del, double &x, double &sum, double (&pcritrh)[6], long &j, double &eps, long &tmax){
    for (z = 1; z <= layers; z++) { //z = 1 To layers
        p1 = 0;
        k = 1;
        e = 0; //flow integral
        erh[z][0] = 0; //first array position is layer number
        krh[z][0] = kmaxrh[z];
            
        do {
            p2 = p1 + pinc;
            qtrapvg(olds,a,n,x,kmaxrh,p1,p2,s,del,sum,eps,tmax,t,z,it,tnm,j);
            e = e + s;
            erh[z][k] = e;
            x = p2;
            krh[z][k] = get_vg(a,n,x,kmaxrh,z); //instantaneous k at upper limit of integration = derivative of flow integral (fundamental theorem of calculus)
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000){
                break; //Then Exit Do //avoid crashing for extreme vc//s
            }
        } while (!(krh[z][k - 1] < kmin)); //Loop Until krh(z, k - 1) < kmin
        pcritrh[z] = p2; //end of line for rhizo element z
    } //endfor// z
} //endsub//

/* Soil water status*/
//'gives soil Y in MPa from soil theta/thetasat=x
double soils::get_rvg(double (&a)[6], double (&n)[6], double &x, long &z) {
    //rvg = (x ^ (1 / (1 - 1 / n(z))) + 1) ^ (1 / n(z)) / (x ^ (1 / (n(z) - 1)) * a(z))
    double aa = pow((pow(x, (1 / (1 - 1 / n[z]))) + 1), (1 / n[z]));
    double bb = (pow(x, (1 / (n[z] - 1))) * a[z]);
    return aa / bb;
}

/* Soil water potential*/
// Soil water content from soil water potential in MPa=x   
double soils::get_swc(double (&a)[6], double (&n)[6], double &x, long &z) {
    //swc = (1 / (1 + (a(z) * x) ^ n(z))) ^ (1 - 1 / n(z))
    return pow((1 / (1 + pow((a[z] * x), n[z]))), (1 - 1 / n[z]));
}

/* Soil wetness*/
//'gets predawns accounting for rain, groundwater flow, transpiration, and redistribution via soil and roots
void soils::get_soilwetness(double &drainage, double &runoff, double &waterold, double &x, double (&thetafracres)[6], double (&a)[6],
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
    bool &ground, bool &gs_inGrowSeason, bool &gs_doneFirstDay) {

    double tempDouble = 0.0;
    drainage = 0;
    runoff = 0;
    
    if (dd == 1 || isNewYear) { //if// //'every layer starts at initial % of field capacity
        if (!(useGSData && gs_yearIndex > 0)){
            waterold = 0; //'total root zone water content (mmol m-2 ground area)
            for (z = 0; z <= layers; z++){//z = 0 To layers //'
                x = 10; //'MPa water potential for getting residual thetafrac
                thetafracres[z] = get_swc(a,n,x,z); //'residual thetafrac
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
            dataCells[rowD + dd][colD + dColF_End_watercontent] = waterold * 1000;
            //dSheet.Cells(rowD + dd, colD + dColF_End_watercontent) = waterold * 1000; //'root zone water content in mm m-2 ground area
                                                                                      //'waterold = waterold * 55555556# //'convert m3 water per ground area to mmol water per ground area
                                                                                      // [HNT] starting water now counts as an input
            gs_ar_input[gs_yearIndex] = gs_ar_input[gs_yearIndex] + waterold * 1000;
        }

        if (gs_yearIndex == 0) // if it's the first year, store this as the off-season starting water because we don't have a real value
            gs_ar_waterInitial_OFF[gs_yearIndex] = waterold * 1000;

        // store the initial water content to check how much we consume at the end
        gs_ar_waterInitial[gs_yearIndex] = waterold * 1000;
        // [/HNT]
    } //End if// //'dd=1 if
    
    //'if pet = "y" Or pet = "n" { //if// //'do the original routine...doesn//'t run for PET scenario
    if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'get flows that happened during previous timestep
        for (z = 0; z < layers; z++){ //z = 0 To layers - 1 //'transpiration, root and soil redistribution
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
            for (z = layers; z >= 0; z--){ //z = layers To 0 Step -1 //'start at bottom, go up
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
        if (water[0] <= 0)
            water[0] = 0.00001; //'set lower limit to surface water content
        for (z = 1; z <= layers; z++){ //z = 1 To layers
            if (water[z] < swclimit[z])
                water[z] = swclimit[z]; //'water at limit
        } //for//z

        bool rainOverride = false;
        if (iter_Counter == 0 && tod == 23 && (stage_ID == STAGE_ID_HIST_OPT || stage_ID == STAGE_ID_FUT_OPT))
            rainOverride = true;

        //predawns mode
        if (mode_predawns)
            rainOverride = true;
        // end predawns mode
        
        //'now check for rain during PREVIOUS TIME STEP
        if (dataCells[rowD + dd - 1][colD + dColRain] >= 0.0 || rainOverride == true) { //if// //'there//'s been rain!
            rain = dataCells[rowD + dd - 1][colD + dColRain] * 0.001; //'convert mm of rain to m depth...equivalent to m3/m2 water volume per ground area
                                                                         //'if raining(xx) = "n" { //if// rain = 0 //'rain turned off
            if (rainOverride){ // just force everything to FC on every timestep
                for (z = 0; z <= layers; z++){//z = 0 To layers //'add in rain
                    water[z] = fc[z];
                }
            } // do this before the rain addition, so that we see all rain as runoff

            if (rainEnabled == false || mode_predawns) // predawns mode also turns off rain
                rain = 0;
            // [HNT] temp
            //if (rain > 0.1) // more than 1/10th meter per hour is a data anomoly, however this is a poor assumption and these should be checked in the weather curator
            //   rain = 0.0;
            // [/HNT]
            //'sumrain = rain * 1000 + sumrain //'total rain input in mm m-2
            for (z = 0; z <= layers; z++){ //z = 0 To layers //'add in rain
                if (rain <= 0)
                    break; //'rain//'s used up
                deficit = fc[z] - water[z]; //'m3 of water per m2 to wet up layer to fc
                                           //if (deficit >= -0.001) { //if// //'layer is at or below field capacity
                if (rain - deficit >= 0) { //if// //'there//'s enough to wet the layer
                    water[z] = fc[z];
                    rain = rain - deficit;
                    //'sumrain = sumrain + deficit //'absorbed rain
                    drainage = rain * 1000.0; //'any left over will drain out the bottom unless layer is rising above field capacity
                } else {
                    water[z] = water[z] + rain; //'it all goes into the first layer
                                                //'sumrain = sumrain + rain
                    rain = 0; //'rain used up
                    drainage = 0;
                } //End if// //'wetting up to field capacity "if"
                 //}
            } //for//z

            // If there's drainage, and the ground water is on, that should be used to fill up to saturation
            if (rain > 0.0 && ground == true) { // the remaining "drainage" is also still stored in the rain variable
                // this is kind of inefficient, but the rain routine actually drained all the layers to FC even if GW was on and we should have been filling to sat
                // now we start at the bottom and fill the layers to saturation using the drainage
                for (j = layers; j >= 0; j--){//j = z - 1 To 0 Step -1 //'go back up to fill profile to saturation
                    if (rain <= 0) //if// Exit for
                        break;
                    deficit = thetasat[j] * depth[j] - water[j];
                    if (deficit >= 0) { //if// //'got capacity
                        if (rain - deficit >= 0) { //if// //'enough rain to saturate the layer
                            water[j] = thetasat[j] * depth[j]; //'saturate the layer
                            rain = rain - deficit; //'reduce rain
                        } else { //'rain absorbed by layer
                            water[j] = water[j] + rain;
                            rain = 0; //'use up rain
                        } //End if// //'deficit=>0 "if"
                    } else { //'deficit<0...layer//'s saturated
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
        
        for (z = 1; z <= layers; z++){//z = 1 To layers
            if (water[z] < swclimit[z])
                layer[z] = 1; //'water exhausted
        } //for//z

        //'now get water content change over PREVIOUS time step        
        if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'now get updated new water content
            waternew = 0;
            for (z = 0; z <= layers; z++){ //z = 0 To layers //'check for exhausted layers
               waternew = water[z] + waternew;
               //'Cells(16 + dd, 42 + z) = water[z]
            } //for//z

            //'waternew = waternew * 55555556# //'new water content at beginning of timestep
            waterchange = waternew - waterold; //'total water change in mm m-2 ground area
                                               //'waterchange = waterchange * 1 / baperga * 1 / laperba //'total water in mmol per m2 leaf area
            dataCells[rowD + dd][colD + dColF_End_watercontent] = waternew * 1000;
            //dSheet.Cells(rowD + dd, colD + dColF_End_watercontent) = waternew * 1000; //'root zone water content in mm

            // always store the water content as the "final" -- not worth testing if it's really the last day of the year
            gs_ar_waterFinal[gs_yearIndex] = waternew * 1000;
            dataCells[rowD + dd][colD + dColF_End_waterchange] = waterchange * 1000;
            //dSheet.Cells(rowD + dd, colD + dColF_End_waterchange) = waterchange * 1000; //'change in water content over PREVIOUS timestep
            //'if raining(xx) = "y" { //if// Cells(16 + dd, 59) = Cells(16 + dd - 1, 4) //'rain input per previous timestep
            if (rainEnabled == true)// && dSheet.Cells(rowD + dd - 1, colD + dColRain) < 100.0)
                dataCells[rowD + dd][colD + dColF_End_rain] = dataCells[rowD + dd - 1][colD + dColRain];
                //dSheet.Cells(rowD + dd, colD + dColF_End_rain) = dSheet.Cells(rowD + dd - 1, colD + dColRain); //'rain input per previous timestep
            dataCells[rowD + dd][colD + dColF_End_gwater] = gwflow;
            //dSheet.Cells(rowD + dd, colD + dColF_End_gwater) = gwflow; //'groundwater input in mm per timestep
                                                                       //'Cells(16 + dd, 61) = transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 //'transpiration per ground area in mm m-2 per timestop
            dataCells[rowD + dd][colD + dColF_End_drainage] = drainage;
            //dSheet.Cells(rowD + dd, colD + dColF_End_drainage) = drainage; //'total drainage in mm per timestep
            dataCells[rowD + dd][colD + dColF_End_input] = dataCells[rowD + dd][colD + dColF_End_rain] + dataCells[rowD + dd][colD + dColF_End_gwater];
            //dSheet.Cells(rowD + dd, colD + dColF_End_input) = dSheet.Cells(rowD + dd, colD + dColF_End_rain) + dSheet.Cells(rowD + dd, colD + dColF_End_gwater); //'total input per timestep in mm
                                                                                                                                                                 //'Cells(16 + dd, 64) = water(0) * 1000 //'water in top layer in mm
            gs_ar_input[gs_yearIndex] = gs_ar_input[gs_yearIndex] + dataCells[rowD + dd][colD + dColF_End_input];// dSheet.Cells(rowD + dd, colD + dColF_End_input);

            dataCells[rowD + dd][colD + dColF_End_runoff] = runoff * 1000;
            //dSheet.Cells(rowD + dd, colD + dColF_End_runoff) = runoff * 1000; //'excess root zone water per timestep in mm
            waterold = waternew; //'reset to beginning of current timestep
        } //End if// //'dd>1 if
    } //End if// //'dd>1 if
    //'} //End if////'pet if
    
    if (dd > 1 && !isNewYear) { //if//
        tempDouble = transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000;
        dataCells[rowD + dd][o + colD + dColF_End_E] = tempDouble;
        //dSheet.Cells(rowD + dd, o + colD + dColF_End_E) = tempDouble;//transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000; //'transpiration per ground area in mm m-2 per timestop
        if (gs_inGrowSeason) // only record growing season E (should not be any non-GS E, but just for safety)
            gs_ar_E[gs_yearIndex] = gs_ar_E[gs_yearIndex] + tempDouble;

        dataCells[rowD + dd][o + colD + dColF_End_soilEvap] = (soilevap * 1) / (998.2 * timestep * 1000);
        //dSheet.Cells(rowD + dd, o + colD + dColF_End_soilEvap) = soilevap * 1 / 998.2 * timestep * 1000; //'evaporative water loss in mm per timestep
        tempDouble = (transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 + soilevap * 1) / (998.2 * timestep * 1000);
        dataCells[rowD + dd][o + colD + dColF_End_ET] = tempDouble;
        //dSheet.Cells(rowD + dd, o + colD + dColF_End_ET) = tempDouble;//transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 + soilevap * 1 / 998.2 * timestep * 1000;
        gs_ar_ET[gs_yearIndex] = gs_ar_ET[gs_yearIndex] + tempDouble;
        tempDouble = atree * timestep * 3600 * 0.001;
        dataCells[rowD + dd][o + colD + dColF_End_ANet] = tempDouble;
        //dSheet.Cells(rowD + dd, o + colD + dColF_End_ANet) = tempDouble;//atree * timestep * 3600 * 0.001; //'Anet per timestep in mmoles per leaf area
        
        if (!std::isnan(tempDouble) && gs_inGrowSeason) { // only record growing season A
            gs_ar_Anet[gs_yearIndex] = gs_ar_Anet[gs_yearIndex] + tempDouble;
            //anytime we record A, also record the ci-related outputs
            // TODO CRIT units???
            // Anet in mmoles per leaf area as in final output columns? (calc above)
            if (night == "n") {// daytime only
                gs_ar_cica[gs_yearIndex] += cinc / ca;          // these are in mols/mol, but it's a ratio so not important
                gs_ar_cica_N[gs_yearIndex]++;
                gs_ar_Aci[gs_yearIndex] += tempDouble * cinc;   // for Aci what units? pa would be (cinc * patm * 1000)
                gs_ar_AnetDay[gs_yearIndex] += tempDouble;      // keep a daytime-only Anet tally
                                                                //this is for A-weighted Ci
                                                                //gs_ar_Acica[gs_yearIndex] += tempDouble * cinc / ca; // for this one, another ratio so ignore units
            }
        }
    } //End if// //'dd>1 if

    // dSheet.Cells(rowD + dd - 3, colD + dColF_CP_kplant)
    if (tod == 16 && !gs_doneFirstDay && gs_inGrowSeason && dataCells[rowD + dd - 3][colD + dColF_CP_kplant] > 0.000000001) { //if// //'get midday k//'s for day 1
        // VPD zero case -- if the stomata did not open on the first day of the GS, kplant won't have been set and will be zero... in which case, try again tomorrow
        gs_doneFirstDay = true;
        sum = 0;
        for (z = 1; z <= 3; z++){//z = 1 To 3
            //sum = sum + dSheet.Cells(rowD + dd - z, colD + dColF_CP_kplant);
            sum = sum + dataCells[rowD + dd - z][colD + dColF_CP_kplant];
        } //for//z
        kpday1 = sum / 3.0; //'average midday kplant on day one
        sum = 0;
        for (z = 1; z <= 3; z++){ //z = 1 To 3
            //sum = sum + dSheet.Cells(rowD + dd - z, colD + dColF_CP_kxylem);
            sum = sum + dataCells[rowD + dd - z][colD + dColF_CP_kxylem];
        } //for//z
        kxday1 = sum / 3.0; //'average midday kxylem on day one

        // [HNT] first day of the new growing season! Record the starting water
        //gs_ar_waterInitial_GS[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent);
        gs_ar_waterInitial_GS[gs_yearIndex] = dataCells[rowD + dd][colD + dColF_End_watercontent];   // no matter what happened above, this has been set by this point
                                                                                                        // and record this as the FINAL water for THIS YEAR's off season -- remember that this year's off season is the PRECEDING winter
        //gs_ar_waterFinal_OFF[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent);
        gs_ar_waterFinal_OFF[gs_yearIndex] = dataCells[rowD + dd][colD + dColF_End_watercontent];
        // [/HNT]
    } //End if// //'dd=16 if
    
    if (gs_doneFirstDay) { //if// //'calculate plc relative to midday of day 1
        if (iter_refK < 0.000000001){ // || iter_Counter == 0 // no longer appropriate to test for iter_Counter == 0 ... may be doing a seperate stress profile that refers to a saved refK, which will have been loaded at start of modelProgramMain
            //tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1); 
            tempDouble = 100 * (1 - dataCells[rowD - 1 + dd][colD + dColF_CP_kplant] / kpday1); //' done to avoid repeating this calculation
        } else { // if we haven't loaded a ref K it will be set to zero, and we need to fall back to the old method.. otherwise use refK to calculate PLC
            //tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / iter_refK);
            tempDouble = 100 * (1 - dataCells[rowD - 1 + dd][colD + dColF_CP_kplant] / iter_refK); // PLCp calculated from refK if it exists
            if (tempDouble < 0.0) // if we're using the refK, it's possible for this to be negative briefly -- should be considered zero
                tempDouble = 0.0;
        }

        //dSheet.Cells(rowD + dd, colD + dColF_End_PLCplant) = tempDouble;
        dataCells[rowD + dd][colD + dColF_End_PLCplant] = tempDouble; //'100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1) 'plc plant...prior timestep
                                                                                // no matter what, add it to the tally for calculating the mean over-season PLCp
        gs_ar_PLCSum[gs_yearIndex] = gs_ar_PLCSum[gs_yearIndex] + tempDouble;   // yearly running total of PLC values
        gs_ar_PLCSum_N[gs_yearIndex] = gs_ar_PLCSum_N[gs_yearIndex] + 1;        // total hours in GS
                                                                                // now test for highest PLC and hours > 85
        if (tempDouble > gs_ar_PLCp[gs_yearIndex])
            gs_ar_PLCp[gs_yearIndex] = tempDouble;
        
        if (tempDouble > 85.0)
            gs_ar_PLC85[gs_yearIndex] = gs_ar_PLC85[gs_yearIndex] + 1;

        //tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1);
        tempDouble = 100 * (1 - dataCells[rowD - 1 + dd][colD + dColF_CP_kxylem] / kxday1);
        //dSheet.Cells(rowD + dd, colD + dColF_End_PLCxylem) = tempDouble; //'100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1) 'plc xylem...prior timestep
        dataCells[rowD + dd][colD + dColF_End_PLCxylem] = tempDouble;

        if (tempDouble > gs_ar_PLCx[gs_yearIndex])
            gs_ar_PLCx[gs_yearIndex] = tempDouble;

        //dSheet.Cells(rowD + dd, colD + dColF_End_PLCplant) = (100.0 * (1.0 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1)); //'plc plant...prior timestep
        //dSheet.Cells(rowD + dd, colD + dColF_End_PLCxylem) = (100.0 * (1.0 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1)); //'plc xylem...prior timestep

        // [HNT] keep track of in-season input
        if (gs_inGrowSeason) {
            // if we've done the first GS day and we're in the growing season then it's "this year" during GS
            //gs_ar_waterInput_GS[gs_yearIndex] += dSheet.Cells(rowD + dd, colD + dColF_End_input);
            gs_ar_waterInput_GS[gs_yearIndex] += dataCells[rowD + dd][colD + dColF_End_input];
        } else {
            // if we've done the first GS day but we're NOT in the growing season, it's the winter following GS.
            // this is considered NEXT YEAR's off-season input! First check if we already have a value for the FINAL water for THIS YEAR, because if it's zero then
            // this is the first timestep of the winter and we need to store it
            if (gs_ar_waterFinal_GS[gs_yearIndex] <= 0.0 && gs_ar_waterInitial_OFF[gs_yearIndex + 1] <= 0.0) { // ok to use == or <= with double here because this will be memset to 0 if it hasn't been set
                //gs_ar_waterFinal_GS[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent); // ok to overshoot by 1 hour, I think (instead of using dd - 1)
                gs_ar_waterFinal_GS[gs_yearIndex] = dataCells[rowD + dd][colD + dColF_End_watercontent]; // ok to overshoot by 1 hour, I think (instead of using dd - 1)
                //gs_ar_waterInitial_OFF[gs_yearIndex + 1] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent); // also the off-season initial for next year
                gs_ar_waterInitial_OFF[gs_yearIndex + 1] = dataCells[rowD + dd][colD + dColF_End_watercontent]; // also the off-season initial for next year           
            } else {// otherwise, we're in the middle of NEXT YEAR's off season ... note this +1 on an array index is super lazy and bad. Make sure to never run exactly this array size # of years
                //gs_ar_waterInput_OFF[gs_yearIndex + 1] += dSheet.Cells(rowD + dd, colD + dColF_End_input); // add the stored input to the input tally, for NEXT YEAR
                gs_ar_waterInput_OFF[gs_yearIndex + 1] += dataCells[rowD + dd][colD + dColF_End_input]; // add the stored input to the input tally, for NEXT YEAR
            }
        }
        // [/HNT] //End if// //'dd>16 if
    // have NOT done first day and are NOT in growing season
    } else if (!gs_inGrowSeason){
        // we must be in the pre-GS winter of what we called the NEXT year above... so gs_yearIndex is now that year, and we add to the off-season input for THIS year
        //gs_ar_waterInput_OFF[gs_yearIndex] += dSheet.Cells(rowD + dd, colD + dColF_End_input);
        gs_ar_waterInput_OFF[gs_yearIndex] += dataCells[rowD + dd][colD + dColF_End_input];
    }
}

/* Pre-Dawn Water Potential*/
void soils::get_predawns(double (&kminroot)[6], double &theta, double (&water)[6], double (&depth)[6], double (&thetasat)[6],
    double &x, double (&pd)[6], double (&dataCells)[DATAFILE_MAXROWS][DATAFILE_MAXCOLS], double &pgrav, double (&a)[6], 
    double (&n)[6], double (&prh)[6], double (&pcritrh)[6], double (&pcritr)[6], double &sum, double &pr, double &prinitial,
    long &k, long &layers, long (&layer)[6], long &z, long &rowD, long &dd, long &colD, long &dColRain, long &t, long &failure,
    long &dColF_p1, long &o, std::string (&layerfailure)[6], std::string &failspot, bool &mode_predawns){
   
    //'first check for layer participation...only rooted layers, not layer 0
    for (k = 1; k <= layers; k++){//k = 1 To layers //'assign source pressures, set layer participation
        if (layerfailure[k] == "root") { //if//
            if (kminroot[k] == 0)
                layer[k] = 1; //'gone from scene if roots cavitated at midday
            if (kminroot[k] != 0){//if// //'roots still around
                layerfailure[k] = "ok";
                layer[k] = 0;
            }//End if//
        } //End if//
        if (layerfailure[k] == "rhizosphere")
            layer[k] = 0; //'layer can come back to life
    }//for//k
   
    //'after getting water[z] and layer participation, get predawns
    for (z = 0; z <= layers; z++){//z = 0 To layers
        if (layer[z] == 0) { //if//
            theta = water[z] / depth[z]; //'convert m3 water per m2 ground back to m3 water / m3 soil
            x = theta / thetasat[z]; //'remember, VG function takes theta/thetasat as input
            //predawns mode
            if (mode_predawns)
                pd[z] = dataCells[rowD + dd][colD + dColRain] - pgrav;
                //pd[z] = dSheet.Cells(rowD + dd, colD + dColRain) - pgrav; // read the predawns+pgrav from the "rain" column
            else 
                pd[z] = get_rvg(a,n,x,z); //'soil pressure of layer
            //end predawns mode changes
            prh[z] = pd[z]; //'guess for NR solution
            if (pd[z] >= pcritrh[z] && z > 0) { //if// //'only rooted layers // [HNT] >= instead of > for consistency w/ Newton Rhapson update
                layer[z] = 1;
                layerfailure[z] = "rhizosphere";
            } //End if//
            if (pd[z] >= pcritr[z] && z > 0) { //if// //'only rooted layers // [HNT] >= instead of > for consistency w/ Newton Rhapson update
                layer[z] = 1;
                layerfailure[z] = "root";
                kminroot[z] = 0;
            } //End if//
        } else { //'layer//'s disconnected
            pd[z] = pcritr[z];
            prh[z] = pcritr[z];
        } //End if//
    } //for//z

    //'now get guess of proot
    sum = 0;
    t = 0;
   
    for (k = 1; k <= layers; k++){ //k = 1 To layers
        if (layer[k] == 0) { //if//
            sum = sum + pd[k];
        } else { //'predawn is not seen by the roots
            t = t + 1;
        } //End if//
    } //for//k

    failspot = "no failure";
   
    if (t < layers) { //if//
        pr = sum / (layers - t); //'set unknown proot to average pd
        prinitial = pr; //'store initial value if NR gets off the rails
    } else {
        failure = 1;
    } //End if//
   
    for (z = 1; z <= layers; z++){//z = 1 To layers
        //dSheet.Cells(rowD + dd, colD + dColF_p1 - 1 + o + z) = pd[z]; //'soil pressures by layer (only for rooted layers)
        dataCells[rowD + dd][colD + dColF_p1 - 1 + o + z] = pd[z]; //'soil pressures by layer (only for rooted layers)
    } //for//z
}

/* Soil water vertical flow */
//'gets flow out of each layer via soil, not including the groundwater basement
void soils::get_soilflow(double &store, double (&kmaxrh)[6], double (&kkmax)[6], double (&depth)[6], double (&pd)[6],
    double (&soilredist)[6], double &p1, double &pend, double &e, double &p2, double &pinc, double &olds, double (&a)[6],
    double (&n)[6], double &x, double &s, double &del, double &sum, double &eps, double (&soilf)[6], long &z, long &layers,
    long &tmax, long &t, long &it, long &tnm, long &j){
    
    for (z = 0; z < layers; z++)//z = 0 To layers - 1
      {
         store = kmaxrh[z]; //'store the rhizosphere kmax for now
         kmaxrh[z] = kkmax[z + 1] * 1 / (depth[z] / 2.0 + depth[z + 1] / 2.0); //'reset it to the vertical soil kmax per m2 ground using distance between midpoints of each layer. Use properties of fatter layer.
         if (pd[z] == pd[z + 1])
            soilredist[z] = 0; //'no redistribution (gravity ignored!)
         if (pd[z] != pd[z + 1]) { //if//
            if (pd[z] < pd[z + 1]) { //if// //'flow is out of layer z (positive)
               p1 = pd[z];
               pend = pd[z + 1];
            }
            else { //'flow is into layer z(negative)
               p1 = pd[z + 1];
               pend = pd[z];
            } //End if//
            e = 0; //'flow integral
            do
            {
               p2 = p1 + pinc;
               qtrapvg(olds,a,n,x,kmaxrh,p1,p2,s,del,sum,eps,tmax,t,z,it,tnm,j);
               e = e + s;
               p1 = p2; //'reset p1 for next increment
            } while (!(p2 >= pend || p2 > 5));
            //Loop Until p2 >= pend Or p2 > 5 //'cutoff for really dry soil
            if (pd[z] < pd[z + 1]) { //if// //'flow is out of layer z (positive)
               soilf[z] = e; //'flow in kg hr-1 m-2 ground area
            }
            else { //'flow is into layer z(negative)
               soilf[z] = -e;
            } //End if//
         } //End if//
         kmaxrh[z] = store; //'reset rhizosphere kmaxrh (though i don//'t THINK it//'s used again?)
      } //for//z
      soilredist[0] = soilf[0]; //'set upper layer which has only one flux
      soilredist[layers] = -1 * soilf[layers - 1]; //'set water flowing into/out of top of the bottom layer...soilredist is net soil flow of layer
                                                   //'now calculate net flows of internal layers
      for (z = 1; z < layers; z++)//z = 1 To layers - 1
      {
         soilredist[z] = -1 * soilf[z - 1] + soilf[z]; //'add up to get net flow...remember negative flow is flow into the layer
      } //for//z
}

/* Deepsoil water flow*/
//'gets flow into/out of the lowermost layer
void soils::get_deepflow(double &store, double (&kmaxrh)[6], double (&kkmax)[6], double &grounddistance, double &pground,
    double (&pd)[6], double &p1, double &pend, double &e, double &p2, double &pinc, double &olds, double (&a)[6],
    double (&n)[6], double &x, double &s,double &del, double &sum, double &eps, double &groundflow, double &gwflow, double &timestep,
    double &drainage, double (&soilredist)[6],long &z, long &layers, long &tmax, long &t, long &it, long &tnm, long &j){
    z = layers; //'just the bottom layer
    store = kmaxrh[z]; //'store the rhizosphere kmax for now
    kmaxrh[z] = kkmax[z] / grounddistance; //'reset it to the vertical soil kmax using distance to groundwater
    if (pd[z] >= pground) { //if// //'groundwater will flow in (negative flow)
        p1 = pground;
        pend = pd[z];
        e = 0; //'flow integral
        do{
            p2 = p1 + pinc;
            qtrapvg(olds,a,n,x,kmaxrh,p1,p2,s,del,sum,eps,tmax,t,z,it,tnm,j);
            e = e + s;
            p1 = p2; //'reset p1 for next increment
        } while (!(p2 >= pend));
        //Loop Until p2 >= pend
        groundflow = -1 * e; //'negative flow (into bottom layer)
        gwflow = -1 * groundflow * 1 / 998.2 * timestep * 1000; //'groundwater flow totals in mm m-2 ground area per time step
    } else { //'flow is out of bottom layer(positive)
        p1 = pd[z];
        pend = pground;
        e = 0;//'flow integral
        do {
            p2 = p1 + pinc;
            qtrapvg(olds,a,n,x,kmaxrh,p1,p2,s,del,sum,eps,tmax,t,z,it,tnm,j);
            e = e + s;
            p1 = p2; //'reset p1 for next increment
        } while (!(p2 >= pend));
        //Loop Until p2 >= pend
        groundflow = e; //'positive flow (out of bottom layer)
        drainage = drainage + groundflow * 1 / 998.2 * timestep * 1000; //'drainage totals in mm m-2 ground area per time step
    } //End if//
    kmaxrh[z] = store; //'reset rhizosphere kmaxrh (though i don//'t THINK it//'s used again?)
    //'now calculate net flow into lowermost layer
    soilredist[z] = soilredist[z] + groundflow; //'add up to get net flow...remember negative flow is flow into the layer
}

/* Soil evaporation*/
// gets soil evaporation based on understory vegetation height and shade
void soils::get_soilevaporation(double &emission, double &emiss, const double &sbc, double& soiltemp, double& lc,
    double (&leaftempsh)[100001], double &rabssoil, double &soilabssol, double &sbottom, double &mdensair,
    double &patm, double &airtemp, double &gha, double &us, double &xh, double &zdispl, double &rough, double &zh,
    double &soilep, const double &sha, double &grad, double &lambda, double &rha, double &vpd, double &maxvpd, 
    double &rhs, double (&pd)[6], const double &gas, double &soilevap, double &baperga, double (&soilredist)[6],
    long &haltsh){
    emission = emiss * sbc * pow((soiltemp + 273.15), 4); //'long wave emissivity of soil
    lc = 0.97 * sbc * pow((leaftempsh[haltsh] + 273.15), 4); //'thermal emissivity from bottom of canopy
    rabssoil = soilabssol * sbottom + 0.97 * lc; //'rabs of soil surface; does not account for slope (beam only)
    mdensair = 44.6 * (patm / 101.3) * 273.15 / (airtemp + 273.15); //'molar density of air, CN 3.3, moles/m3
    gha = 0.16 * mdensair * us / (log((xh - zdispl) / rough) * log((xh - zdispl) / zh)); //'CN 7.28/14.9,heat and vapor aerodynamic conductance for soil, mol m-2 s-1, stability factors = 0
    soilep = (rabssoil - emission - sha * (soiltemp - airtemp) * (grad + gha)) / lambda; //'potential evaporation in moles m-2 s-1; 12.8
    if (soilep < 0)
        soilep = 0; //'set to zero if negative
    rha = 1 - vpd / maxvpd; //'relative humidity of air
    rhs = exp((-pd[0] * 0.00001805) / (gas * 0.000001 * (soiltemp + 273.15))); //'relative humidity in equilibrium with soil surface water potential
    if (rha == 1) { //if//
        soilevap = 0;
    } else {
        soilevap = soilep * (rhs - rha) / (1 - rha); //'campbell 1985, eqn 9.14; actual soil evaporation rate in moles m-2s-1
    } //End if//
    if (soilevap < 0)
        soilevap = 0; //'don//'t let it go negative
    soilevap = soilevap * (1 - baperga); //'reduce for basal area
    soilevap = soilevap * 0.0180153 * 3600; //'converts from moles m-2s-1 to kg m-2 hr-1
    soilredist[0] = soilredist[0] + soilevap; //'add evaporative loss to redistribution for top layer
}
