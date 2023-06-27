// Functions to extract weibull parameters and fit weibull function
// TO-DO: merge stem and leaves integration functions as these could be just generic functions
// TO-DO: after merging stem and leaves integration functions check possibility of merging root integration functions
#include "03Hydraulics.h"
#include "05CAssimilation.h"

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
void hydraulics::trapzdwbr(const long& t, double &p1, double &p2, const double& b, const double& c,double* ksatr,
   const int& z, double& s, long& it, long& tnm, double& del, double& x, double& sum) {
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

void hydraulics::qtrapwbr(double& olds, long& t, double& p1, double& p2, const double& b, const double& c,double* ksatr,
   const int& z, double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx){ //'evaluates accuracy of root element z integration
   
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
void hydraulics::get_rootPcrit(double er[][100001],double kr[][100001],long& z, const int& layers,double* ksatr, double& p1,
   double& p2, const double& pinc,long& k,double& e,double& olds, long& t, const double& b, const double& c, double& s,
   long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx, const double& kmin, double* pcritr) { //generates fresh global E(P) curve for the element(erases history)
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

void hydraulics::get_rootPcrit_v(double er_v[][100001],double kr_v[][100001],long& z, const int& layers,double* ksatr,
   double& p1, double& p2, const double& pinc,long& k,double& e,double& olds, long& t, const double& b, const double& c,
   double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx, const double& kmin,
   double* pcritr) { //generates fresh global E(P) curve for the element(erases history)
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
   double kr[6][100001],double& ehigh,  double& khigh, double& estart, double& klower, double& p2, double& efinish,
   double& kupper, double& flow, long& k, long& z){
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
   double& frt, double& dfrdpr, double* pcritrh, double& p1, double& p2, double& plow, double& pinc, double& elow,
   double erh[6][100001], double& klow, double krh[6][100001],double& ehigh, double& khigh, double& estart, double& klower,
   double& efinish, double& kupper, double& flow, double* func, double* dfrhdprh, double* pcritr, double er[6][100001],
   double kr[6][100001], double* dfrhdpr, double* dfrdprh, double& e, double& sum, double& threshold, double& initialthreshold,
   double& aamax, double* vv, double& dum, double* indx, double& pcrits, double& waterold, long* gs_ar_nrFailConverge,
   double* gs_ar_nrFailConverge_Water, double* gs_ar_nrFailConverge_WaterMax, long& weird, long& check, long* layer,
   long& k, long& ticks, long& unknowns, long& d, long& imax, long& ii,long& ll, long& gs_yearIndex, long& dd,
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
void hydraulics::trapzdwbs(const long &t, double &p1, double &p2, const double& b, const double& c, const double& ksats,
   double &s, long& it, long& tnm, double& del, double& x, double& sum){ //integrates stem element weibull
   
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

void hydraulics::qtrapwbs(double& olds, long& t, const long& f, double &p1, double &p2, const double& b, const double& c,
   const double& ksats,double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx){ //evaluates accuracy of root element z integration
   
   olds = -1; //starting point unlikely to satisfy if statement below
   
   for (t = 1; t <= f; t++){
      
      trapzdwbs(t, p1, p2, b, c, ksats, s, it, tnm, del, x, sum);
      
      if (std::abs(s - olds) <= epsx * std::abs(olds))
         return;
      olds = s;
   }
}

// stem E(P) curve
void hydraulics::get_stemPcrit(double* es,bool vCurve, double* es_v, double &p1, double &p2, const double& pinc,long& k,
   double& e,double& olds, long &t, const long& f, const double& b, const double& c, const double& ksats,double &s,
   long& it, long& tnm, double& del, double& x, double& sum, double& epsx, double& ksh, const double& kmin, double& pcrits) {
   
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
void hydraulics::get_stem(double &p1, const double& pr, const double& pgrav, double& plow, const double& pinc, double* es,
   double& elow, double& ehigh, double& estart, double& efinish, double& e, double& p2, double& ps, double& pcrits, long& k,
   long& test, std::string& failspot){
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
void hydraulics::trapzdwbl(long &t, double &p1, double &p2, const double& b, const double& c, const double& ksatl, double &s,
   long& it, long& tnm, double& del, double& x, double& sum) {
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

void hydraulics::qtrapwbl(double& olds, long& t, const long& f,double &p1, double &p2, const double& b, const double& c,
   const double& ksatl,double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx) {
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
void hydraulics::get_leafPcrit(double* el, bool vCurve, double* el_v, double &p1, double &p2, const double& pinc,long& k,
   double& e, double& olds, long& t, const long& f, const double& b, const double& c, const double& ksatl,double &s,
   long& it, long& tnm, double& del, double& x, double& sum, double& epsx, double& ksh, const double& kmin, double& pcritl) {
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
void hydraulics::get_leaf(double &p1, const double& ps, const double& pgrav, double& plow, const double& pinc, double* el,
   double& elow, double& ehigh, double& estart, double& efinish, double& e, double& p2, double& pl, double& pcritl,
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

// WHOLE PLANT block
// composite curves: stores composite E(P)curve and the element conductances
void hydraulics::get_compositecurve(double elayer[6][100001], double* prh, double prhizo[6][100001], double& plow,
   double& p1, double& pinc, double& elow, double er[6][100001], double& klow, double kr[6][100001],double& ehigh,
   double& khigh, double& estart, double& klower, double& p2, double& efinish, double& kupper, double& flow,
   double& pr, double kroot[6][100001], double& x, double* pd, const double& root_b, const double& root_c, double* ksatr,
   double* kminroot, double* pcritrh, double* proot,double* pstem, double& ps, double* pleaf, double& pl, double* kleaf,
   double* kstem, double* kplant, double& kminleaf, double& kminstem, double& kminplant, double& e, double& pgrav,
   const double& leaf_b, const double& leaf_c, double& ksatl, const double& stem_b, const double& stem_c, double& ksats, double* eplant,
   double& einc, double* dedp, double* dedpf, double& pcritsystem, double& ecritsystem, long& k, long& p, long* layer,
   long& test, long& total, const int& layers, const bool& refilling, std::string* layerfailure){
   elayer[0][p] = 0; //'no root mediated flow in topmost layer
   for (long z = 1; z <= layers; z++){//z = 1 To layers
      if (layer[z] == 0){ //if//
         prhizo[z][p] = prh[z];
         p1 = prh[z];
         p2 = pr;
         get_rootflow(plow,p1,pinc,elow,er,klow,kr,ehigh,khigh,estart,klower,p2,efinish,
            kupper,flow,k,z);
         elayer[z][p] = flow; //'flow through layer
         if (flow != 0){
            kroot[z][p] = std::abs(elayer[z][p] / (pr - prhizo[z][p]));
         }
         if (flow == 0) { //if//
            if (refilling == true) { //if// //'for refilling, starting point is always weibull
               x = pd[z];
               kroot[z][p] = get_weibullfitroot(x,root_b,root_c,ksatr,z);
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
            kroot[z][p] = get_weibullfitroot(x,root_b,root_c,ksatr,z); //'root element conductance = instantaneous conductance from weibull curve at pr
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
      if (refilling == true) { //if//
         kleaf[p] = kminleaf;
         kstem[p] = kminstem;
         kplant[p] = kminplant;
      } //End if//
      if (refilling == false) { //if//
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

/* stores historical element e(P) curves*/
void hydraulics::storehistory(double ter[6][100001], double er[6][100001], double kr[6][100001], double tkr[6][100001], double er_v[6][100001],
   double kr_v[6][100001], double* tes, double* es, double* es_v, double* tel, double* el, double* el_v, double* kminroot,
   long& k, long* layer, long* tlayer, int& layers, std::string* layerfailure, std::string* tlayerfailure){
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

   for (long z = 1; z <= layers; z++){
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
   } while (es[k] != 0 || es_v[k] != 0 || tes[k] != 0);
   //Loop Until es(k) = Empty
   k = 0;
   do{
      tel[k] = el[k];
      el[k] = el_v[k];
      k = k + 1;
      if (k == 100000){
         break;
      }
   } while (el[k] != 0 || el_v[k] != 0 || tel[k] != 0);
   //Loop Until el(k) = Empty
   for (long z = 1; z <= layers; z++){
      tlayer[z] = layer[z];
      tlayerfailure[z] = layerfailure[z];
   }
   for (long z = 1; z <= layers; z++){// 'get all layers functioning to start with
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

/* restores historical element e(P) curves*/
void hydraulics::gethistory(double ter[6][100001], double er[6][100001], double kr[6][100001], double tkr[6][100001], double* tes,
   double* es, double* tel, double* el, long& k, long* layer, long* tlayer, int& layers,
   std::string* layerfailure, std::string* tlayerfailure){
   for (long z = 1; z <= layers; z++){
      k = 0;
      do{
         er[z][k] = ter[z][k];
         kr[z][k] = tkr[z][k];
         k = k + 1;
         if (k == 100000)
               break;
      } while (ter[z][k] != 0 || er[z][k] != 0);//Loop Until ter[z][k] = Empty
   }
   k = 0;
   do{
      es[k] = tes[k];
      k = k + 1;
      if (k == 100000)
            break;
   } while (tes[k] != 0 || es[k] != 0);//Loop Until tes[k] = Empty
   k = 0;
   do{
      el[k] = tel[k];
      k = k + 1;
      if (k == 100000)
         break;
   } while (tel[k] != 0 || el[k] != 0);//Loop Until tel[k] = Empty
   //'re-set failure status (this at the critical point)
   for (long z = 1; z <= layers; z++){
      layer[z] = tlayer[z];
      layerfailure[z] = tlayerfailure[z];
   }
}

/* Get canopy pressure*/
void hydraulics::get_canopypressure(const double& ecritsystem, double ter[6][100001], double er[6][100001], double kr[6][100001],
   double tkr[6][100001], double er_v[6][100001], double kr_v[6][100001], double* tes, double* es, double* es_v, double* tel,
   double* el, double* el_v, double* kminroot, double& sum, double* prh, double* pd, double* pcritr, double& pr, double& psynmaxmd,
   double& psynmaxshmd, double& e, double& einc, double& dedplmin, double& ksatp, double jmatrix[7][7], double& frt, double& dfrdpr,
   double* pcritrh, double& p1, double& p2, double& plow, double& pinc, double& elow, double erh[6][100001], double& klow,
   double krh[6][100001],double& ehigh, double& khigh, double& estart, double& klower, double& efinish, double& kupper, double& flow,
   double* func, double* dfrhdprh, double* dfrhdpr, double* dfrdprh, double& threshold, double& initialthreshold, double& aamax,
   double* vv, double& dum, double* indx, double& pcrits, double& waterold, double* gs_ar_nrFailConverge_Water,
   double* gs_ar_nrFailConverge_WaterMax, const double& pgrav, double& ps, double& pl, double& pcritl, double* pleafv, double& predawn,
   double& plold, double& dedplzero, double& dedpl, double* klossv, double& pcritsystem, double& rabs, const double& ssun, const double& sref,
   const double& emiss, const double& la, const double& lg, double& lambda, const double& airtemp, double& grad, double& gha, const double& wind, 
   const double& leafwidth, double& emd, const double& laperba, const double& sbc, double& numerator, double& denominator, const double& sha,
   double& leaftmd, double& lavpdmd, const double& patm, const double& vpd, const double& sshade, double& leaftshmd, double& lavpdshmd,
   double& gcanwmd, const double& gmax, double& gcancmd, double& comp, const double& comp25, const double& gas, const double& svvmax, 
   const double& havmax, const double& hdvmax, const double& vmax25, double& vmax, const double& svjmax, const double& hdjmax,
   const double& hajmax, double& jmax, const double& jmax25, const double& kc25, const double& ko25, double& kc, double& ko,
   double& rday25, double& rdaymd, double& ci, double& jact, const double& qmax, const double& qsl, const double& lightcurv,
   double& je, const double& oa, double& jc, double& var, const double& thetac, const double& ca, double* psynmd, double& cinmd,
   double& marker, double& gcanwshmd, double& gcancshmd, double& rdayshmd, const double& qsh, double* psynshmd, double& cinshmd,
   double& maxkloss, double& dpmax, double& dpamax, double& amaxmax, double& dpamin, double* amaxfrac, double* dpa, double& rmean,
   double& md, double* lavpd, double& dpasun, double& mdsh, double* amaxfracsh, double* lavpdsh, double* pleaf, double* eplantl,
   double& transpiration, double& psynact, double* psyn, double& gcmd, double* gcanw, double& cinc, double* cin, double& transpirationsh,
   double& psynactsh, double* psynsh, double& gcmdsh, double* gcanwsh, double& lavpdmdsh, double& cincsh, double* cinsh, double& kminstem,
   double& kminleaf, double kroot[6][100001], double* kstem, double* kleaf,
   long& k, long* layer, long* tlayer, long& t, long& failure, long& test, long& p, long* gs_ar_nrFailConverge, long& weird, long& check,
   long& ticks, long& unknowns, long& d, long& imax, long& ii,long& ll, long& gs_yearIndex, long& dd, long& totalv,
   const long& runmean, long& total, const double& cutoff, long& halt, long haltsh,
   int& layers, std::string* layerfailure, std::string* tlayerfailure, std::string& failspot, const std::string& night,
   const bool& refilling){
   cassimilation photosynfunctions;
   //computes carbon-based middays; pressures and cost curve from virgin supply function, gas exchange from historical values
   //store history, get MD and cost function from virgin curve
   if (ecritsystem == 0){
      k = 0;
      goto tenMarker; //don't bother to find midday
   }
   storehistory(ter,er,kr,tkr,er_v,kr_v,tes,es,es_v,tel,el,el_v,kminroot,k,layer,tlayer,layers,layerfailure,tlayerfailure); //stores xylem element curves and failure status, resets layers to full functioning
   //rootcurves(); //erases history for md solution
   //stemcurve();
   //leafcurve();
   sum = 0;
   t = 0;
   for (long k = 1; k <= layers; k++){//assign source pressures, set layer participation
      if (layer[k] == 0){
         prh[k] = pd[k]; //initial guess of unknown rhizosphere pressures
         sum = sum + pd[k];
      } else {
         prh[k] = pcritr[k];
         t = t + 1;
      }
   }
   if (t < layers){
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
   do{
      e = e + einc;
      p = p + 1;
      newtonrhapson(kminroot, pr, pd, prh,jmatrix,frt, dfrdpr, pcritrh, p1, p2, plow, pinc,
         elow,erh, klow,krh,ehigh, khigh, estart, klower,efinish, kupper, flow,
         func, dfrhdprh,pcritr,er,kr,dfrhdpr,dfrdprh, e, sum, threshold,
         initialthreshold, aamax, vv, dum, indx, pcrits, waterold,
         gs_ar_nrFailConverge, gs_ar_nrFailConverge_Water, gs_ar_nrFailConverge_WaterMax,
         weird,check,layer,k,ticks,unknowns,d,imax,ii,ll,gs_yearIndex,dd,
         layers,layerfailure,failspot); //pd's already assigned...this solves for p's and e's in fingers of chicken foot

      if (check >= 500){
         //p = p - 1
         //std::cout << "NR failed on virgin curve at e = " << e << " on timestep dd = " << dd << std::endl;
         break; //gone as far as can go
      }
      //gets stem and leaf pressures
      get_stem(p1,pr,pgrav,plow,pinc,es,elow,ehigh,estart,efinish,e,p2,ps,pcrits,k,test,failspot);
      get_leaf(p1,ps,pgrav,plow,pinc,el,elow,ehigh,estart,efinish,e,p2,pl,pcritl,k,test,failspot);
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
      }
      //now get virgin A curve
      photosynfunctions.get_leaftempsmd(rabs,ssun,sref,emiss,la,lg,lambda,airtemp,grad,
         gha,wind,leafwidth,emd,e,laperba,sbc,numerator, denominator,sha,leaftmd,lavpdmd,
         patm,vpd);//gets virgin sun layer leaf temperature from energy balance
      photosynfunctions.get_leaftempsshademd(rabs,sshade,sref,emiss,lg,emd,lambda,airtemp,
         sbc,numerator,denominator,sha,grad,gha,leaftshmd,lavpdshmd,patm,vpd);//gets virgin shade layer leaf temperature
      photosynfunctions.get_assimilationmd(lavpdmd,gcanwmd,gmax,emd,gcancmd,comp,comp25,gas,
         leaftmd,numerator,svvmax,havmax,denominator,hdvmax,vmax25,vmax,svjmax,hdjmax,hajmax,
         jmax,jmax25,kc25,ko25,kc,ko,rday25,rdaymd,ci,jact,qmax,qsl,lightcurv,je,oa,jc,var,
         thetac,ca,psynmd,cinmd,marker,psynmaxmd,p,night); //gets virgin sun layer photosynthesis
      photosynfunctions.get_assimilationshademd(lavpdshmd,gcanwshmd,gmax,emd,gcancshmd,comp,comp25,gas,
         leaftshmd,numerator,svvmax,hdvmax,havmax,denominator,vmax25,vmax,svjmax,hdjmax,hajmax,jmax,
         jmax25,kc25,ko25,kc,ko,rday25,rdayshmd,ci,jact,qmax,qsh,lightcurv,je,oa,jc,var,thetac,ca,psynshmd,
         cinshmd,marker,psynmaxshmd,p,night); //gets virgin shade layer photosynthesis
      //by now we have assigned psynmd(p) and psynshmd(p) and reset psynmaxes
      plold = pl;
      //e = e + einc
      //p = p + 1
   } while (!(test == 1 || p > 99900)); //loop to failure or out of "p"s//Loop Until test = 1 Or p > 99900 'loop to failure or out of "p"s
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
      if (p >= runmean - 1){//get running mean
         sum = 0;
         for (long i = p - runmean + 1; i <= p; i++){
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
   if (psynmaxshmd == 0){//shade layer's below light compensation point, don't open
      mdsh = pleafv[0]; //set midday = predawn
   } else {
      p = -1; //p is still index for virgin curve
      dpmax = 0;
      dpamax = -100;
      amaxmax = 0;
      dpamin = 0.0;
      do{
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
         if (p >= runmean - 1){//get running mean
            sum = 0;
            for (long i = p - runmean + 1; i <= p; i++){
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
   do{
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
   gethistory(ter,er,kr,tkr,tes,es,tel,el,k,layer,tlayer,layers,layerfailure,tlayerfailure); //reinstates historical element curves and failure status prior to updating
   if (refilling == true){//need to record midday kmins, uses sunlit pressures
      for (long z = 1; z <= layers; z++){
         kminroot[z] = kroot[z][halt];
      }
      kminstem = kstem[halt];
      kminleaf = kleaf[halt];
   }
}

/* Resets E(P) element curves*/
void hydraulics::update_curves(double kroot[6][100001],double* kminroot, const double& pinc, double* proot,double er[6][100001],
   double kr[6][100001], double* kstem, double& kminstem, double* pstem, double* es, double* kleaf, double& kminleaf,
   double* pleaf, double* el,long& halt, long& phigh, const int& layers){
   //'if k<kmin, re-assign e//'s on element curve by back-calculating
   for (long z = 1; z <= layers; z++){//z = 1 To layers
      if (true){
         if(kroot[z][halt] < kminroot[z]){
            kminroot[z] = kroot[z][halt];
            phigh = int(proot[halt] / pinc) + 1; //'pressure datum just above the target
            for (long k = phigh; k >= 0; k--){//k = phigh To 0 Step -1 //'back-calculate e//'s
               er[z][k] = er[z][phigh] - kminroot[z] * pinc * (phigh - k);
               kr[z][k] = kminroot[z]; //'back-calculate KR(Z,K) too for roots (not stem or leaves)
            } //EndFor  k
         } //EndIf//
      }
   } //EndFor  z
   
   if (kstem[halt] < kminstem){
      kminstem = kstem[halt];
      phigh = int(pstem[halt] / pinc) + 1;
      for (long k = phigh; k >= 0; k--){//k = phigh To 0 Step -1 //'back-calculate e//'s
         es[k] = es[phigh] - kminstem * pinc * (phigh - k);
      } //EndFor  k
   } //EndIf//
   if (kleaf[halt] < kminleaf){
      kminleaf = kleaf[halt];
      phigh = int(pleaf[halt] / pinc) + 1;
      for (long k = phigh; k >= 0; k--){//k = phigh To 0 Step -1 //'back-calculate e//'s
         el[k] = el[phigh] - kminleaf * pinc * (phigh - k);
      } //EndFor  k
   } //EndIf//
   //'if kplant[halt] < kminplant Then kminplant = kplant[halt]NOTE: kplant CAN go up because of rhizosphere recovery!
}