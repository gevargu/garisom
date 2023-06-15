// Functions to extract weibull parameters and fit weibull function
#include "03Hydraulics.h"
#include "00MainProgram.h"

// Weibull parameter c
double hydraulics::get_cweibull(double p12, double p50) {
   // c = (log(log(0.88)/log(0.50)))/(log(p12)-log(p50))
   double px1;
   double px2;
   px1 = p12*-1; // absolute value of p12
   px2 = p50*-1; // absolute value of p50
   return (log(log(0.88)/log(0.5)))/(log(px1)-log(px2));
}

// Weibull parameter b
double hydraulics::get_bweibull(double p12, double c) {
   // b = p12/(-log(0.88))^(1/c)
   double px1;
   px1 = p12*-1; // absolute value of p12
   return px1/(pow(-log(0.88),(1/c)));
}

// Weibull function for stem and leaves
double hydraulics::get_weibullfit(double &x, double b, double c, double ksat) {
   //wb = ksat * Exp(-((x / b) ^ c))
   return ksat * exp(-(pow((x / b), c)));
}

// ROOT block
// % resistance in roots
double hydraulics::get_rootpercent(double leafpercent){
   return 2.0 / 3.0 * (100.0 - leafpercent);
}

// kmax of root system; assumes zero % rhizosphere resistance in WET soil
double hydraulics::get_rootkmax(double rootpercent, double rsatp) {
   return 1.0 / ((rootpercent / 100.0) * rsatp);
}
// Weibull function for root elements
double hydraulics::get_weibullfitroot(double &x, double b, double c, double ksatr[6], int z){
   //wb = ksat * Exp(-((x / b) ^ c))
   return ksatr[z] * exp(-(pow((x / b), c)));
}

//integrates root element z weibull
void hydraulics::trapzdwbr(double &p1, double &p2, double &s, long &t) {
   MainProgram garisom;
   if (t == 1){
      double vg_p1 = get_weibullfitroot(p1, garisom.root_b, garisom.root_c, garisom.ksatr, garisom.z);
      double vg_p2 = get_weibullfitroot(p2, garisom.root_b, garisom.root_c, garisom.ksatr, garisom.z);
      s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
      garisom.it = 1;
   } else {
      garisom.tnm = garisom.it;
      garisom.del = (p2 - p1) / garisom.tnm;
      garisom.x = p1 + 0.5 * garisom.del;
      garisom.sum = 0;
      
      for (garisom.j = 1; garisom.j <= garisom.it; garisom.j++) {
         garisom.sum = garisom.sum + get_weibullfitroot(garisom.x, garisom.root_b, garisom.root_c, garisom.ksatr, garisom.z);
         garisom.x = garisom.x + garisom.del;
      }
      s = 0.5 * (s + (p2 - p1) * garisom.sum / garisom.tnm);
      garisom.it = 2 * garisom.it;
   }
   //[HNT]
   //testcounterIntRoot = testcounterIntRoot + 1
   //[\HNT]
}

void hydraulics::qtrapwbr(double &p1, double &p2, double &s){ //'evaluates accuracy of root element z integration
   MainProgram garisom;
   garisom.olds = -1; //'starting point unlikely to satisfy if statement below
   
   for (garisom.t = 1; garisom.t <= garisom.f; garisom.t++) {
      trapzdwbr(p1, p2, s, garisom.t);
      if (std::abs(s - garisom.olds) <= garisom.epsx * std::abs(garisom.olds))
         return;
      garisom.olds = s;
   }
}

// Max resistance in the rhizosphere
double hydraulics::get_kmaxrh(double rhizor, double vgterm){
   return (1.0 / rhizor) / vgterm;
}

// STEM block
// integrates stem element z weibull
void hydraulics::trapzdwbs(double &p1, double &p2, double &s, long &t){ //integrates root element z weibull
   MainProgram garisom;
   if (t == 1) {
      double vg_p1 = get_weibullfit(p1,garisom.stem_b,garisom.stem_c,garisom.ksats);
      double vg_p2 = get_weibullfit(p2,garisom.stem_b,garisom.stem_c,garisom.ksats);
      s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
      garisom.it = 1;
   } else {
      garisom.tnm = garisom.it;
      garisom.del = (p2 - p1) / garisom.tnm;
      garisom.x = p1 + 0.5 * garisom.del;
      garisom.sum = 0;
      for (garisom.j = 1; garisom.j <= garisom.it; garisom.j++){
         garisom.sum = garisom.sum + get_weibullfit(garisom.x,garisom.stem_b,garisom.stem_c,garisom.ksats);
         garisom.x = garisom.x + garisom.del;
      }
      s = 0.5 * (s + (p2 - p1) * garisom.sum / garisom.tnm);
      garisom.it = 2 * garisom.it;
   }
}

void hydraulics::qtrapwbs(double &p1, double &p2, double &s){ //evaluates accuracy of root element z integration
   MainProgram garisom;
   garisom.olds = -1; //starting point unlikely to satisfy if statement below
   
   for (garisom.t = 1; garisom.t <= garisom.f; garisom.t++){
      
      trapzdwbs(p1, p2, s, garisom.t);
      
      if (std::abs(s - garisom.olds) <= garisom.epsx * std::abs(garisom.olds))
         return;
      garisom.olds = s;
   }
}

// % resistance in stems
double hydraulics::get_stempercent(double leafpercent){
   return 1.0 / 3.0 * (100.0 - leafpercent);
}

// kmax of stem system
// Call setValueFromName("o_ksatRoot", ksatroot) //Cells(2, 2) = ksatroot //whole root system kmax
// Call setValueFromName("o_ksatStem", ksats) //Cells(2, 3) = ksats //stem network kmax
// Call setValueFromName("o_ksatLeaf", ksatl) //Cells(2, 4) = ksatl //parallel kmax for leaves
// dbh = ((ksatp * 0.785) / Cells(6, 9)) ^ (1 / (Cells(7, 9) - 2)) //basal diameter in m from k by D allometry
// height = (99 / Cells(5, 9)) * dbh ^ (2 / 3) //height from H by D allometry and safety factor
// Cells(3, 9) = height //height in m
// Cells(3, 7) = dbh
double hydraulics::get_stemkmax(double stempercent, double rsatp){
   return 1.0 / ((stempercent / 100.0) * rsatp);
}


// LEAVES block
// gmax per leaf area in mmol m-2s-1
double hydraulics::get_gmaxl(double gmax, double laperba){
   return gmax * (1.0 / laperba) * (1 / 3600.0) * 55.56 * 1000.0;
}
// Leaf conductance per basal area
double hydraulics::get_leafkmax(double ksatp, double leafpercent){
   return ksatp * (100.0 / leafpercent);
}
// Leaf specific conductance in kg hr-1m-2MPa-1
double hydraulics::get_LSC(double ksatl, double laperba){
   //lsc in kg hr-1m-2MPa-1
   //lsc per leaf area in kg hr - 1 
   //Call setValueFromName("o_leafLSC", lsc / (3600 * 0.00001805)) 
   //Cells(2, 9) = lsc / (3600 * 0.00001805) //lsc converted to mmol s - 1m - 2MPa - 1
   return ksatl * 1.0 / laperba; 
}