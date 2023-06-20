// Functions to extract weibull parameters and fit weibull function
#include "03Hydraulics.h"

// Weibull parameter c
double hydraulics::get_cweibull(const double& p12, const double& p50) {
   // c = (log(log(0.88)/log(0.50)))/(log(p12)-log(p50))
   double px1;
   double px2;
   px1 = p12*-1; // absolute value of p12
   px2 = p50*-1; // absolute value of p50
   return (log(log(0.88)/log(0.5)))/(log(px1)-log(px2));
}

// Weibull parameter b
double hydraulics::get_bweibull( const double& p12, const double& c) {
   // b = p12/(-log(0.88))^(1/c)
   double px1;
   px1 = p12*-1; // absolute value of p12
   return px1/(pow(-log(0.88),(1/c)));
}

// Weibull function for stem and leaves
double hydraulics::get_weibullfit(const double& x,const double& b, const double& c, const double& ksat) {
   //wb = ksat * Exp(-((x / b) ^ c))
   return ksat * exp(-(pow((x / b), c)));
}

// ROOT block
// % resistance in roots
double hydraulics::get_rootpercent(const double& leafpercent){
   return 2.0 / 3.0 * (100.0 - leafpercent);
}

// kmax of root system; assumes zero % rhizosphere resistance in WET soil
double hydraulics::get_rootkmax(const double& rootpercent,const double& rsatp) {
   return 1.0 / ((rootpercent / 100.0) * rsatp);
}
// Weibull function for root elements
double hydraulics::get_weibullfitroot(const double& x, const double& b, const double& c, double* ksatr, const int& z){
   //wb = ksat * Exp(-((x / b) ^ c))
   return ksatr[z] * exp(-(pow((x / b), c)));
}

//integrates root element z weibull
void hydraulics::trapzdwbr(double &p1, double &p2, double &s, long &t) {
   if (t == 1){
      double vg_p1 = get_weibullfitroot(p1, root_b, root_c, ksatr, z);
      double vg_p2 = get_weibullfitroot(p2, root_b, root_c, ksatr, z);
      s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
      it = 1;
   } else {
      tnm = it;
      del = (p2 - p1) / tnm;
      x = p1 + 0.5 * del;
      sum = 0;
      
      for (j = 1; j <= it; j++) {
         sum = sum + get_weibullfitroot(x, root_b, root_c, ksatr, z);
         x = x + del;
      }
      s = 0.5 * (s + (p2 - p1) * sum / tnm);
      it = 2 * it;
   }
   //[HNT]
   //testcounterIntRoot = testcounterIntRoot + 1
   //[\HNT]
}

void hydraulics::qtrapwbr(double &p1, double &p2, double &s){ //'evaluates accuracy of root element z integration
   
   olds = -1; //'starting point unlikely to satisfy if statement below
   
   for (t = 1; t <= f; t++) {
      trapzdwbr(p1, p2, s, t);
      if (std::abs(s - olds) <= epsx * std::abs(olds))
         return;
      olds = s;
   }
}

   void rootcurves() //generates fresh global E(P) curve for the element(erases history)
   {
      //clear arrays from earlier calls
      memset(er, 0, sizeof(er));//Erase er
      memset(kr, 0, sizeof(kr));//Erase kr
                                //do root elements
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         kr[z][0] = ksatr[z]; //kr is instantaneous K from weibull function
                              //kminr(z) = ksatr(z)
         p1 = 0;
         k = 1;
         e = 0; //value of integral
         er[z][0] = 0; //first array position is layer number
         do
         {
            p2 = p1 + pinc;
            qtrapwbr(p1, p2, s);
            e = e + s;
            er[z][k] = e;
            x = p2;
            kr[z][k] = wbr(x); //weibull k at upper limit of integration = derivative of flow integral
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000)
               break;
            //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
         } while (!(kr[z][k - 1] < kmin));
         //Loop Until kr(z, k - 1) < kmin
         pcritr[z] = p2; //end of line for element z
      } //endfor// z
   } //endsub//

   void rootcurves_v() //generates fresh global E(P) curve for the element(erases history)
   {
      //clear arrays from earlier calls
      memset(er_v, 0, sizeof(er_v));//Erase er_v
      memset(kr_v, 0, sizeof(kr_v));//Erase kr_v
                                //do root elements
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         kr_v[z][0] = ksatr[z]; //kr_v is instantaneous K from weibull function
                              //kminr(z) = ksatr(z)
         p1 = 0;
         k = 1;
         e = 0; //value of integral
         er_v[z][0] = 0; //first array position is layer number
         do
         {
            p2 = p1 + pinc;
            qtrapwbr(p1, p2, s);
            e = e + s;
            er_v[z][k] = e;
            x = p2;
            kr_v[z][k] = wbr(x); //weibull k at upper limit of integration = derivative of flow integral
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000)
               break;
            //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
         } while (!(kr_v[z][k - 1] < kmin));
         //Loop Until kr_v(z, k - 1) < kmin
         pcritr[z] = p2; //end of line for element z
      } //endfor// z
   } //endsub//




// STEM block
// integrates stem element z weibull
void hydraulics::trapzdwbs(double &p1, double &p2, double &s, long &t){ //integrates root element z weibull
   
   if (t == 1) {
      double vg_p1 = get_weibullfit(p1,stem_b,stem_c,ksats);
      double vg_p2 = get_weibullfit(p2,stem_b,stem_c,ksats);
      s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
      it = 1;
   } else {
      tnm = it;
      del = (p2 - p1) / tnm;
      x = p1 + 0.5 * del;
      sum = 0;
      for (j = 1; j <= it; j++){
         sum = sum + get_weibullfit(x,stem_b,stem_c,ksats);
         x = x + del;
      }
      s = 0.5 * (s + (p2 - p1) * sum / tnm);
      it = 2 * it;
   }
}

void hydraulics::qtrapwbs(double &p1, double &p2, double &s){ //evaluates accuracy of root element z integration
   
   olds = -1; //starting point unlikely to satisfy if statement below
   
   for (t = 1; t <= f; t++){
      
      trapzdwbs(p1, p2, s, t);
      
      if (std::abs(s - olds) <= epsx * std::abs(olds))
         return;
      olds = s;
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