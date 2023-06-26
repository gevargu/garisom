// Functions to extract weibull parameters and fit weibull function
// TO-DO: merge stem and leaves integration functions as these could be just generic functions
// TO-DO: after merging stem and leaves integration functions check possibility of merging root integration functions

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
void hydraulics::trapzdwbr(const long& t, double &p1, double &p2, const double& b, const double& c,double* ksatr, const int& z, double& s, long& it, long& tnm, double& del, double& x, double& sum) {
   if (t == 1){
      double vg_p1 = get_weibullfitroot(p1, b, c, ksatr, z);
      double vg_p2 = get_weibullfitroot(p2, b, c, ksatr, z);
      s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
      it = 1;
   } else {
      tnm = it;
      del = (p2 - p1) / tnm;
      x = p1 + 0.5 * del;
      sum = 0;
      
      for (long j = 1; j <= it; j++) {
         sum = sum + get_weibullfitroot(x, b, c, ksatr, z);
         x = x + del;
      }
      s = 0.5 * (s + (p2 - p1) * sum / tnm);
      it = 2 * it;
   }
   //[HNT]
   //testcounterIntRoot = testcounterIntRoot + 1
   //[\HNT]
}

void hydraulics::qtrapwbr(double& olds, long& t, double& p1, double& p2, const double& b, const double& c,double* ksatr, const int& z, double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx){ //'evaluates accuracy of root element z integration
   
   olds = -1; //'starting point unlikely to satisfy if statement below
   
   for (t = 1; t <= f; t++) {
      trapzdwbr(t, p1, p2, b, c, ksatr, z, s, it, tnm, del, x, sum);
      if (std::abs(s - olds) <= epsx * std::abs(olds)){
         return;
      }
      olds = s;
   }
}

// root E(P) curves
void hydraulics::get_rootPcrit(double er[][100001],double kr[][100001],long& z, const int& layers,double* ksatr, double& p1, double& p2, const double& pinc,long& k,double& e,double& olds, long& t, const double& b, const double& c, double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx, const double& kmin,double* pcritr) { //generates fresh global E(P) curve for the element(erases history)
   //clear arrays from earlier calls
   memset(er, 0, sizeof(er));//Erase er
   memset(kr, 0, sizeof(kr));//Erase kr
   
   //do root elements
   for (z = 1; z <= layers; z++) {//z = 1 To layers
      kr[z][0] = ksatr[z]; //kr is instantaneous K from weibull function
      //kminr(z) = ksatr(z)
      p1 = 0;
      k = 1;
      e = 0; //value of integral
      er[z][0] = 0; //first array position is layer number
      
      do {
         p2 = p1 + pinc;
         qtrapwbr(olds, t, p1, p2, b, c, ksatr, z, s, it, tnm, del, x, sum, f, epsx);
         e = e + s;
         er[z][k] = e;
         x = p2;
         kr[z][k] = get_weibullfitroot(x,b,c,ksatr,z); //weibull k at upper limit of integration = derivative of flow integral
         p1 = p2; //reset p1 for next increment
         k = k + 1;
         if (k == 100000){
            break;//If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
         }
      } while (!(kr[z][k - 1] < kmin)); //Loop Until kr(z, k - 1) < kmin
      
      pcritr[z] = p2; //end of line for element z
   } //endfor// z
} //endsub//

void hydraulics::get_rootPcrit_v(double er_v[][100001],double kr_v[][100001],long& z, const int& layers,double* ksatr, double& p1, double& p2, const double& pinc,long& k,double& e,double& olds, long& t, const double& b, const double& c, double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx, const double& kmin,double* pcritr) { //generates fresh global E(P) curve for the element(erases history)
   //clear arrays from earlier calls
   memset(er_v, 0, sizeof(er_v));//Erase er_v
   memset(kr_v, 0, sizeof(kr_v));//Erase kr_v
   
   //do root elements
   for (z = 1; z <= layers; z++) {//z = 1 To layers 
      kr_v[z][0] = ksatr[z]; //kr_v is instantaneous K from weibull function
      //kminr(z) = ksatr(z)
      p1 = 0;
      k = 1;
      e = 0; //value of integral
      er_v[z][0] = 0; //first array position is layer number
      
      do {
         p2 = p1 + pinc;
         qtrapwbr(olds, t, p1, p2, b, c, ksatr, z, s, it, tnm, del, x, sum, f, epsx);
         e = e + s;
         er_v[z][k] = e;
         x = p2;
         kr_v[z][k] = get_weibullfitroot(x,b,c,ksatr,z); //weibull k at upper limit of integration = derivative of flow integral
         p1 = p2; //reset p1 for next increment
         k = k + 1;
         if (k == 100000){
            break; //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
         }
      } while (!(kr_v[z][k - 1] < kmin)); //Loop Until kr_v(z, k - 1) < kmin
      
      pcritr[z] = p2; //end of line for element z
   } //endfor// z
} //endsub//

// rhizosphere flow using global E(P) curve
void hydraulics::get_rhizoflow(double& plow, double& p1, double& pinc, double& elow, double erh[6][100001], double& klow,
               double krh[6][100001],double& ehigh,  double& khigh, double& estart, double& klower, double& p2,
               double& efinish, double& kupper, double& flow,
               long& k, long& z){
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

// flow through root using global E(P) curve
void hydraulics::get_rootflow(double& plow, double& p1, double& pinc, double& elow, double er[6][100001], double& klow,
              double kr[6][100001],double& ehigh,  double& khigh, double& estart, double& klower, double& p2,
              double& efinish, double& kupper, double& flow, long& k, long& z){
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

//'does LU decomposition on the jacobian prior to solution by lubksb
void hydraulics::ludcmp(double& sum, double& aamax, double jmatrix[7][7], double* vv, double& dum,
            double* indx, long& d, long& unknowns, long& imax){
   d = 1;
   for (int i = 1; i <= unknowns; i++){//i = 1 To unknowns
      aamax = 0;
      for (int j = 1; j <= unknowns; j++){//j = 1 To unknowns
         if (std::abs(jmatrix[i][j]) > aamax){
            aamax = std::abs(jmatrix[i][j]);
         }
      }
      if (aamax == 0){
         return;
      }
      vv[i] = 1 / aamax;
   }
   
   for (int j = 1; j <= unknowns; j++){//j = 1 To unknowns
      
      for (int i = 1; i < j; i++){//i = 1 To j - 1
         sum = jmatrix[i][j];
         for (int k = 1; k < i; k++){//k = 1 To i - 1
               sum = sum - jmatrix[i][k] * jmatrix[k][j];
         }
         jmatrix[i][j] = sum;
      }
      aamax = 0;
      for (int i = j; i <= unknowns; i++){//i = j To unknowns
         sum = jmatrix[i][j];
            for (int k = 1; k < j; k++) //k = 1 To j - 1
            {
               sum = sum - jmatrix[i][k] * jmatrix[k][j];
            }
            jmatrix[i][j] = sum;
            dum = vv[i] * std::abs(sum);
            if (dum > aamax)
            {
               imax = i;
               aamax = dum;
            }
      }
      if (j != imax){ // j <> imax
         for (int k = 1; k <= unknowns; k++){ //k = 1 To unknowns
            
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
         for (int i = j + 1; i <= unknowns; i++){//i = j + 1 To unknowns
            jmatrix[i][j] = jmatrix[i][j] * dum;
         }
      }
   }
   memset(vv, 0, sizeof(vv));
   //Erase vv
}

// solves the decomposed jacobian delta p's
void hydraulics::lubksb(double& sum, double jmatrix[7][7], double* indx, double* func, long& ii, long& unknowns, long& ll){
   ii = 0;
   for (int i = 1; i <= unknowns; i++){//i = 1 To unknowns
      ll = int(indx[i]); //'indx array comes from ludcmp
      sum = func[ll]; //'the func array input is the right-hand vector
      func[ll] = func[i];
      if (ii != 0){
         for (int j = ii; j < i; j++){//j = ii To i - 1
            sum = sum - jmatrix[i][j] * func[j];
         }
      } else {
         if (sum != 0){
            ii = i;
         }
      }
      func[i] = sum;
   }
   for (int i = unknowns; i >= 1; i--){//i = unknowns To 1 Step -1
      sum = func[i];
      for (int j = i + 1; j <= unknowns; j++){//j = i + 1 To unknowns   
         sum = sum - jmatrix[i][j] * func[j];
      }
      func[i] = sum / jmatrix[i][i];
   }
}

// returns rhizosphere pressures and root pressure as a function of pd's and e
void hydraulics::newtonrhapson(double* kminroot, double& pr, double* pd, double* prh, double jmatrix[7][7],
                                 double& frt, double& dfrdpr, double* pcritrh, double& p1, double& p2,
                                 double& plow, double& pinc, double& elow, double erh[6][100001], double& klow,
                                 double krh[6][100001],double& ehigh, double& khigh, double& estart, double& klower,
                                 double& efinish, double& kupper, double& flow, double* func, double* dfrhdprh,
                                 double* pcritr, double er[6][100001], double kr[6][100001], double* dfrhdpr,
                                 double* dfrdprh, double& e, double& sum, double& threshold, double& initialthreshold,
                                 double& aamax, double* vv, double& dum, double* indx, double& pcrits, double& waterold,
                                 long* gs_ar_nrFailConverge, double* gs_ar_nrFailConverge_Water, double* gs_ar_nrFailConverge_WaterMax,
                                 long& weird, long& check, long* layer, long& k, long& ticks, long& unknowns,
                                 long& d, long& imax, long& ii,long& ll, long& gs_yearIndex, long& dd,
                                 int& layers, std::string* layerfailure, std::string& failspot){
   //'prinitial = pr //record the original guess
   weird = 0; //tracks pr estimate
   check = 0; //restart loop counter
   do{//loop to reset guesses
      //'restore layer functioning
      for (long z = 1; z <= layers; z++){
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
         if (false){// can enable this is running into erroneous solutions -- but allowing these to vary instead of resetting results in more frequent solutions in my experience
            // alternatively could randomize them properly
            for (k = 1; k <= layers; k++){//reset prhz(k)
               prh[k] = pd[k];
            }
         }
         weird = 0; //reset cutoff
      }//end reset guesses loop
      check = check + 1; //number of restarts
      ticks = 0; //convergence counter
      do{//loop to seek convergence
         ticks = ticks + 1;
         if (ticks > 1000) {
            weird = 1;
            std::cout << "NR ticks exceeded 1000 -- setting weird=1 for retry. Pinc too high? dd = " << dd << std::endl;
         }
         //'get top row of matrix and right-hand func vector
         //'zero out the jacobian first
         for (long k = 1; k <= unknowns; k++){
            for (long j = 1; j <= unknowns; j++){
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
         for (long z = 1; z <= layers; z++){
            if (pd[z] >= pcritrh[z] || prh[z] >= pcritrh[z]){
               layer[z] = 1; //layer's just gone out of function
               layerfailure[z] = "rhizosphere";
               //weird = 1; //could be result of non-convergence
            }
            if (layer[z] == 0){//it's functional
               p1 = pd[z]; //p1 is not a guess
               p2 = prh[z]; //prh[z] IS a guess and it is initially set prior to the RN routine
               get_rhizoflow(plow,p1,pinc,elow,erh,klow,krh,ehigh,khigh,estart,klower,p2,efinish,kupper,flow,k,z); //gets flows through rhizosphere element from p2-p1 and E(p)curve; gets K's at p1 and p2 as well
               func[z] = flow;
               dfrhdprh[z] = kupper;
            }
            if (prh[z] >= pcritr[z] || pr >= pcritr[z]){
               layer[z] = 1; //layer's just gone out of function
               layerfailure[z] = "root";
               //weird = 1; //could be result of non-convergence
            }
            if (layer[z] == 0){//it's functional
               p1 = prh[z]; //now re-set p1 to prh[z]...the guess
               p2 = pr; //guess must come from before NR routine?
               get_rootflow(plow,p1,pinc,elow,er,klow,kr,ehigh,khigh,estart,klower,p2,efinish,kupper,flow,k,z); //gets flows through root element, and K's at upper and lower bounds
               func[z] = func[z] - flow;
               dfrhdprh[z] = dfrhdprh[z] + klower;
               dfrhdpr[z] = -kupper;
               dfrdprh[z] = -klower;
            }
            if (layer[z] == 1){//layer's out of function...zero out values
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
         for (long k = 1; k <= layers; k++){
            jmatrix[k][unknowns] = dfrhdpr[k]; //last column with dFrh/dPr partials
         }
         for (long k = 1; k <= layers; k++){
               jmatrix[unknowns][k] = dfrdprh[k]; //last row with dFr/dPrh partials
         }
         for (long k = 1; k <= layers; k++){
               jmatrix[k][k] = dfrhdprh[k]; //diagonal of dFrh/dPrh partials
         }
         jmatrix[unknowns][unknowns] = dfrdpr; //lower right corner with dFr/dPr partial
         func[unknowns] = frt; //last position in right-hand flow vector
         //'ok, jacobian and righthand vector are loaded
         //'test for total failure
         sum = 0;
         for (long k = 1; k <= layers; k++){
            sum = sum + layer[k];
         }
         if (sum == layers){//total failure
            failspot = "belowground";
            weird = 1; //trigger a restart
         }
         //'test for flow conservation (steady-state)
         threshold = 0;
         for (int k = 1; k <= unknowns; k++){//k = 1 To unknowns
            threshold = threshold + std::abs(func[k]);
         }
         if (ticks == 1){
            initialthreshold = threshold;
         }
         //'remember to replace "n" with "unknowns" in ludcmp and lubksb
         ludcmp(sum,aamax,jmatrix,vv,dum,indx,d,unknowns,imax); //numerical recipe for doing LU decomposition of jacobian prior to solving
         lubksb(sum,jmatrix,indx,func,ii,unknowns,ll); //solves the decomposed jacobian for delta p's
         //'print out solution vector of pressures
         //'revise unknown pressures
         for (int k = 1; k <= layers; k++){//k = 1 To layers
            prh[k] = prh[k] - func[k]; //NOTE lubksb replaces original right-side func()vector with the solution vector
         }
         pr = pr - func[unknowns];
         //'check for jumping lower bound
         for (int k = 1; k <= layers; k++){//k = 1 To layers
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
               weird = 1; //
            }
         }
      } while (!(threshold < 0.02 || weird == 1));//Loop Until threshold < 0.01 Or weird = 1 //weird = 1 restarts NR
   } while (!((threshold < 0.02 && weird == 0) || check > 500));

   if (check > 500){// disable this output if it's causing too much spam
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
   for (int z = 1; z <= layers; z++) {
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

// STEM block
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

// integrates stem elemente weibull
void hydraulics::trapzdwbs(const long &t, double &p1, double &p2, const double& b, const double& c, const double& ksats,double &s, long& it, long& tnm, double& del, double& x, double& sum){ //integrates stem element weibull
   
   if (t == 1) {
      double vg_p1 = get_weibullfit(p1,b,c,ksats);
      double vg_p2 = get_weibullfit(p2,b,c,ksats);
      s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
      it = 1;
   } else {
      tnm = it;
      del = (p2 - p1) / tnm;
      x = p1 + 0.5 * del;
      sum = 0;
      for (long j = 1; j <= it; j++){
         sum = sum + get_weibullfit(x,b,c,ksats);
         x = x + del;
      }
      s = 0.5 * (s + (p2 - p1) * sum / tnm);
      it = 2 * it;
   }
}

void hydraulics::qtrapwbs(double& olds, long& t, const long& f, double &p1, double &p2, const double& b, const double& c, const double& ksats,double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx){ //evaluates accuracy of root element z integration
   
   olds = -1; //starting point unlikely to satisfy if statement below
   
   for (t = 1; t <= f; t++){
      
      trapzdwbs(t, p1, p2, b, c, ksats, s, it, tnm, del, x, sum);
      
      if (std::abs(s - olds) <= epsx * std::abs(olds))
         return;
      olds = s;
   }
}

// stem E(P) curve
void hydraulics::get_stemPcrit(double* es,bool vCurve, double* es_v, double &p1, double &p2, const double& pinc,long& k, double& e,double& olds, long &t, const long& f, const double& b, const double& c, const double& ksats,double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx, double& ksh, const double& kmin, double& pcrits) {
   
   double* es_ptr = es;
   
   if (vCurve) {
      memset(es_v, 0, sizeof(es_v));
      es_ptr = es_v;
   } else {
      memset(es, 0, sizeof(es));
   } //memset(es, 0, sizeof(es));//Erase es //eliminate values from previous calls
   
   p1 = 0;
   k = 1;
   e = 0; //value of integral
   es_ptr[0] = 0;
   
   do {
      p2 = p1 + pinc;
      qtrapwbs(olds, t, f, p1, p2, b, c, ksats, s, it, tnm, del, x, sum, epsx);
      e = e + s;
      es_ptr[k] = e;
      x = p2;
      ksh = get_weibullfit(x,b,c,ksats); //weibull k
      p1 = p2; //reset p1 for next increment
      k = k + 1;
      if (k == 100000){
         break; //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
      }
   } while (!(ksh < kmin)); //Loop Until ksh < kmin
   
   pcrits = p2;//end of line for element z
} //endsub//

// gets stem pressure and conductance from pr and e.
void hydraulics::get_stem(double &p1, const double& pr, const double& pgrav, double& plow, const double& pinc, double* es, double& elow, double& ehigh,
               double& estart, double& efinish, double& e, double& p2, double& ps, double& pcrits,
               long& k, long& test, std::string& failspot){
   //'start with stem
   p1 = pr + pgrav; //add gravity drop before integration
   plow = int(p1 / pinc); //pressure index below target
   k = int(p1 / pinc);
   elow = es[k]; //e below target
   ehigh = es[k + 1]; //e above target
   plow = plow * pinc; //convert index to pressure
   estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //linear interpolation of starting e
   efinish = estart + e;
   
   long j = k;
   //find efinish
   do{
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

// LEAVES block
// gmax per leaf area in mmol m-2s-1
double hydraulics::get_gmaxl(const double& gmax, const double& laperba){
   return gmax * (1.0 / laperba) * (1 / 3600.0) * 55.56 * 1000.0;
}
// Leaf conductance per basal area
double hydraulics::get_leafkmax(const double& ksatp, const double& leafpercent){
   return ksatp * (100.0 / leafpercent);
}
// Leaf specific conductance in kg hr-1m-2MPa-1
double hydraulics::get_LSC(const double& ksatl, const double& laperba){
   //lsc in kg hr-1m-2MPa-1
   //lsc per leaf area in kg hr - 1 
   //Call setValueFromName("o_leafLSC", lsc / (3600 * 0.00001805)) 
   //Cells(2, 9) = lsc / (3600 * 0.00001805) //lsc converted to mmol s - 1m - 2MPa - 1
   return ksatl * 1.0 / laperba; 
}

// integrate leaf element weibull
void hydraulics::trapzdwbl(long &t, double &p1, double &p2, const double& b, const double& c, const double& ksatl, double &s, long& it, long& tnm, double& del, double& x, double& sum) {
   if (t == 1) {
      double vg_p1 = get_weibullfit(p1,b,c,ksatl);
      double vg_p2 = get_weibullfit(p2,b,c,ksatl);
      s = 0.5 * (p2 - p1) * (vg_p1 + vg_p2);
      it = 1;
   } else {
      tnm = it;
      del = (p2 - p1) / tnm;
      x = p1 + 0.5 * del;
      sum = 0;
      for (long j = 1; j <= it; j++) {
         sum = sum + get_weibullfit(x,b,c,ksatl);
         x = x + del;
      }
      s = 0.5 * (s + (p2 - p1) * sum / tnm);
      it = 2 * it;
   }
}

void hydraulics::qtrapwbl(double& olds, long& t, const long& f,double &p1, double &p2, const double& b, const double& c, const double& ksatl,double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx) {
   olds = -1; //starting point unlikely to satisfy if statement below
   for (t = 1; t <= f; t++) {
      trapzdwbl(t, p1, p2, b, c, ksatl, s, it, tnm, del, x, sum);
      if (std::abs(s - olds) <= epsx * std::abs(olds)) {
         return;
      }
      olds = s;
   }
}

// leaves E(P) curve
void hydraulics::get_leafPcrit(double* el, bool vCurve, double* el_v, double &p1, double &p2, const double& pinc,long& k, double& e, double& olds, long& t, const long& f, const double& b, const double& c, const double& ksatl,double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx, double& ksh, const double& kmin, double& pcritl) {
   double* el_ptr = el;
   if (vCurve) { 
      memset(el_v, 0, sizeof(el_v));
      el_ptr = el_v;
   } else {
      memset(el, 0, sizeof(el));
   } //memset(el_ptr, 0, sizeof(el_ptr));//Erase el_ptr //eliminate values from previous calls
   
   p1 = 0;
   k = 1;
   e = 0; //value of integral
   el_ptr[0] = 0;
   do {
      p2 = p1 + pinc;
      qtrapwbl(olds, t, f, p1, p2, b, c, ksatl, s, it, tnm, del, x, sum, epsx);
      e = e + s;
      el_ptr[k] = e;
      x = p2;
      ksh = get_weibullfit(x,b,c,ksatl); //weibull k
      p1 = p2; //reset p1 for next increment
      k = k + 1;
      
      if (k == 100000){
         break;
      } //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
   } while (!(ksh < kmin)); //Loop Until ksh < kmin
   pcritl = p2; //end of line for element z
} //endsub//


// gets leaf pressure from stem pressure
void hydraulics::get_leaf(double &p1, const double& ps, const double& pgrav, double& plow, const double& pinc, double* el, double& elow, double& ehigh,
               double& estart, double& efinish, double& e, double& p2, double& pl, double& pcritl,
               long& k, long& test, std::string& failspot){
   p1 = ps;
   plow = int(p1 / pinc); //pressure index below target
   k = int(p1 / pinc);
   elow = el[k]; //e below target
   ehigh = el[k + 1]; //e above target
   plow = plow * pinc; //convert index to pressure
   estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //linear interpolation of starting e
   efinish = estart + e;
   long j = k;
   //find efinish
   do{
      j = j + 1;
      if (el[j] == 0){//= Empty 
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

