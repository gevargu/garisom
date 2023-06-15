// Functions to calculate soil parameters at start
#include "02Soils.h"
#include "00MainProgram.h"

/* Soil structure */
// Layer depths and % root ksat
double soils::get_depths(long k,double layers, double beta){
    //Cells(8 + k, 2) = 100 * (0.995 / layers) //equal % roots per layer
    //Cells(8 + k, 11) = 0.01 * Log(1 - k * 0.995 / layers) / Log(beta) //lower depth of each layer converted to m
    return 0.01 * log(1.0 - k * 0.995 / layers) / log(beta);//lower depth of each layer converted to m
    //paramCells[rowLR + k][colLR + 11] = std::to_string(layerDepths[k]);
}

// Maximum depth
//depthmax = Cells(8 + layers, 11) //max depth in meters
//depthmax = lSheet.Cells(rowLR + layers, colLR + 11) //max depth in meters

// Get half depths
double soils::get_halfdepths(long k, double layers, double beta) {
    return 0.01 * log(1 - k * 0.995 / (layers * 2.0)) / log(beta);
}

// Layer volume
double soils::get_layvolume(double depth, double pi, double depthmax, double aspect){
    return depth * pi * pow((depthmax * aspect), 2.0);
}

// Layer radial width
double soils::get_radius(double vol, double depth, double pi) {
    return pow((vol / (depth * pi)), 0.5);
}

/* Rhizosphere resistance
    Solve for what rhizor has to be at the target
*/

double soils::get_rhizor(double rplant, double rhizotarg)
{
    return rplant * (rhizotarg / (1.0 - rhizotarg));
}

/* 
Get Van Genuchten alpha // override if provided
 - Get soil parameters from USDA texture
 - from Leij, Alves van Genuchten, 1996
 */

std::vector <double> soils::get_vgparams(string texture) {
    /* 
    This function calculates the Van Genuchten parameters for a given texture
    Inputs:
    texture = USDA texture classificaion. Object type string. 
    Outputs:
    a = van Genuchten alpha. Object type double.
    n = van Genuchten n. Object type double.
    soilkmax = saturate conductivity of soil kmax in kg h-1 MPa-1 m-1. Object type double.
    thetasat = theta sat in volume/volume. Object type double.
    */
    double a, n,soilkmax, thetasat;
    vector <double> params;
    if (texture == "sand"){
        a = 1479.5945; 
        n = 2.68;
        soilkmax = 30305.88; 
        thetasat = 0.43;
    } else if (texture == "loamy sand"){
        a=1265.3084;
        n=2.28;
        soilkmax=14897.84;
        thetasat=0.41;
    } else if (texture == "sandy loam") {
        a=765.3075;
        n=1.89;
        soilkmax=4510.168;
        thetasat=0.41;
    } else if (texture == "loam") {
        a=367.3476;
        n=1.56;
        soilkmax=1061.216;
        thetasat=0.43;
    } else if (texture == "silt") {
        a=163.2656;
        n=1.37;
        soilkmax=255.1;
        thetasat=0.46;
    } else if (texture == "silt loam") {
        a=204.082;
        n=1.41;
        soilkmax=459.18;
        thetasat=0.45;
    } else if (texture == "sandy clay loam") {
        a=602.0419;
        n=1.48;
        soilkmax=1336.724;
        thetasat=0.39;
    } else if (texture == "clay loam") {
        a=193.8779;
        n=1.31;
        soilkmax=265.304;
        thetasat=0.41;
    } else if (texture == "silty clay loam") {
        a=102.041;
        n=1.23;
        soilkmax=71.428;
        thetasat=0.43;
    } else if (texture == "sandy clay") {
        a=275.5107;
        n=1.23;
        soilkmax=122.448;
        thetasat=0.38;
    } else if (texture == "silty clay") {
        a=51.0205;
        n=1.09;
        soilkmax=20.408;
        thetasat=0.36;
    } else if (texture == "clay") {
        a=81.6328;
        n=1.09;
        soilkmax=204.08;
        thetasat=0.38;
    } else {
        std::cout << "Incorrect soil texture category" << endl;
    }
    params.push_back(a);
    params.push_back(n);
    params.push_back(soilkmax);
    params.push_back(thetasat);

  return params;
}

// EQ1
// vp = 1 / ((a(z) * x) ^ n(z) + 1) 
double soils::get_vgp(double a[6], double n[6], double x, long z)
{
    return 1.0 / (pow((a[z] * x), n[z]) + 1);
}

// EQ2
// vp ^ ((n[z] - 1) / (2 * n[z])) * ((1 - vp) ^ ((n[z] - 1) / n[z]) - 1) ^ 2
double soils::get_vgterm(double n[6], double vp, long z)
{
    return pow(vp, ((n[z] - 1) / (2.0 * n[z]))) * pow((pow((1 - vp), ((n[z] - 1) / n[z])) - 1), 2.0);
}

// EQ3 // does function does everything than the other 2, might keep just 1
double soils::get_vg(double a[6], double n[6], double &x,double kmaxrh[6], long z) //'the van genuchten function for soil k
{
    //vp = 1 / ((a(z) * x) ^ n(z) + 1)
    //vg = kmaxrh(z) * vp ^ ((n(z) - 1) / (2 * n(z))) * ((1 - vp) ^ ((n(z) - 1) / n(z)) - 1) ^ 2
    double vp;
    vp = 1 / (pow((a[z] * x), n[z]) + 1);
    return kmaxrh[z] * pow(vp, ((n[z] - 1) / (2 * n[z]))) * pow((pow((1 - vp), ((n[z] - 1) / n[z])) - 1), 2.0);
}

/* Soil water status*/
//'gives soil Y in MPa from soil theta/thetasat=x
double soils::get_rvg(double a[6], double n[6], double &x, long z) {
    //rvg = (x ^ (1 / (1 - 1 / n(z))) + 1) ^ (1 / n(z)) / (x ^ (1 / (n(z) - 1)) * a(z))
    double aa = pow((pow(x, (1 / (1 - 1 / n[z]))) + 1), (1 / n[z]));
    double bb = (pow(x, (1 / (n[z] - 1))) * a[z]);
    return aa / bb;
}

// Soil water content from soil water potential in MPa=x   
double soils::get_swc(double a[6], double n[6], double &x, long z) {
    //swc = (1 / (1 + (a(z) * x) ^ n(z))) ^ (1 - 1 / n(z))
    return pow((1 / (1 + pow((a[z] * x), n[z]))), (1 - 1 / n[z]));
}
    
/* Van Genuchten Integration */
// integrates van genuchten
void soils::trapzdvg(double &p1, double &p2, double &s, long &t) {
    MainProgram garisom;
    if (t == 1) {
        double vg_p1 = get_vg(garisom.a, garisom.n, p1,garisom.kmaxrh,garisom.z);
        double vg_p2 = get_vg(garisom.a, garisom.n, p2,garisom.kmaxrh,garisom.z);
        s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
        garisom.it = 1;
    } else {
        garisom.tnm = garisom.it;
        garisom.del = (p2 - p1) / garisom.tnm;
        garisom.x = p1 + 0.5 * garisom.del;
        garisom.sum = 0;
        
        for (garisom.j = 1; garisom.j <= garisom.it; garisom.j++) {
            garisom.sum = garisom.sum + get_vg(garisom.a, garisom.n, garisom.x, garisom.kmaxrh,garisom.z);
            garisom.x = garisom.x + garisom.del;
        }
        s = 0.5 * (s + (p2 - p1) * garisom.sum / garisom.tnm);
        garisom.it = 2 * garisom.it;
    }
}

// evaluates accuracy of van genuchten integration
void soils::qtrapvg(double &p1, double &p2, double &s) {
    MainProgram garisom;
    garisom.eps = 0.001; //fractional refinement threshold for integral
    garisom.tmax = 70; //limit of calls to trapzdvg
    garisom.olds = -1; //starting point unlikely to satisfy if statement below
    
    for (garisom.t = 1; garisom.t <= garisom.tmax; garisom.t++){
        trapzdvg(p1, p2, s, garisom.t);
        if (std::abs(s - garisom.olds) < garisom.eps * std::abs(garisom.olds))
        return;
        garisom.olds = s;
    }
}

