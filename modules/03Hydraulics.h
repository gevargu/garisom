#ifndef HYDRAULICS
#define HYDRAULICS

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

using namespace std;
class hydraulics
{
public:
        // general
        double get_cweibull(double &p12, double &p50);
        double get_bweibull(double& p12, double& c);
        double get_weibullfit(double& x, double& b, double &c, double &ksat);
        // root block
        double get_rootpercent(double &leafpercent);
        double get_rootkmax(double &rootpercent, double &rsatp);
        double get_weibullfitroot(double &x, double &b, double &c, double (&ksatr)[6], long &z);
        void trapzdwbr(double &p1, double &p2, double &root_b, double &root_c, double (&ksatr)[6],
                double &x, double &sum, double &del, double &s, long &t, long &z, long &it, long &tnm, long &j);
        void qtrapwbr(double &olds, double &p1, double &p2, double &root_b, double &root_c, double (&ksatr)[6],
                double &x, double &sum, double &del, double &s, double &epsx,long &t, long &f, long &z, long &it,
                long &tnm, long &j);
        void get_rootPcrit(double (&kr)[6][100001], double (&ksatr)[6], double &p1, double &e, double (&er)[6][100001],
                double &p2, double &pinc, double &olds, double &root_b, double &root_c, double &x, double &sum, double &del,
                double &s,double &epsx, double &kmin, double (&pcritr)[6], long &z,long &layers, long &k, long &t, long &f, 
                long &it, long &tnm, long &j);
        void get_rootPcrit_v(double (&kr_v)[6][100001], double (&ksatr)[6], double &p1, double &e, double (&er_v)[6][100001],
                double &p2, double &pinc, double &olds, double &root_b, double &root_c, double &x, double &sum, double &del, 
                double &s, double &epsx, double &kmin, double (&pcritr)[6], long &z,long &layers, long &k, long &t, long &f,
                long &it, long &tnm, long &j);
        void get_rhizoflow(double &plow, double &p1, double &pinc, double &elow, double (&erh)[6][100001], double &klow,
                double (&krh)[6][100001], double &ehigh, double &khigh, double &estart, double &klower, double &p2, long &k,
                double &efinish, double &kupper, double &flow, long &z);
        void get_rootflow(double &plow, double &p1, double &pinc, double &elow, double (&er)[6][100001],
                double &klow, double (&kr)[6][100001], double &ehigh, double &khigh, double &estart, double &klower, double &p2,
                double &efinish, double &kupper, double &flow, long &k, long &z);
        void ludcmp(double& sum, double& aamax, double jmatrix[7][7], double* vv, double& dum,
                double* indx, long& d, long& unknowns, long& imax);
        void lubksb(double& sum, double jmatrix[7][7], double* indx, double* func, long& ii, long& unknowns, long& ll);
        void newtonrhapson(double* kminroot, double& pr, double* pd, double* prh, double jmatrix[7][7],
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
                long& layers, std::string* layerfailure, std::string& failspot);
        
        // stem block
        double get_stempercent(double leafpercent);
        double get_stemkmax(double stempercent, double rsatp);
        void trapzdwbs(const long &t, double &p1, double &p2, const double& b, const double& c, const double& ksats,double &s,
                long& it, long& tnm, double& del, double& x, double& sum);
        void qtrapwbs(double& olds, long &t, const long& f, double &p1, double &p2, const double& b, const double& c,
                const double& ksats,double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx);
        void get_stemPcrit(double* es,bool vCurve, double* es_v, double &p1, double &p2, const double& pinc,long& k, double& e,
                double& olds, long &t, const long& f, const double& b, const double& c, const double& ksats,double &s,
                long& it, long& tnm, double& del, double& x, double& sum, double& epsx, double& ksh, const double& kmin,
                double& pcrits);
        void get_stem(double &p1, const double& pr, const double& pgrav, double& plow, const double& pinc, double* es, double& elow, double& ehigh,
                double& estart, double& efinish, double& e, double& p2, double& ps, double& pcrits,
                long& k, long& test, std::string& failspot);
    
        // leaves block
        double get_gmaxl(const double& gmax, const double& laperba);
        double get_leafkmax(const double& ksatp, const double& leafpercent);
        double get_LSC(const double& ksatl, const double& laperba);
        void trapzdwbl(long &t, double &p1, double &p2, const double& b, const double& c, const double& ksatl, double &s, long& it,
                long& tnm, double& del, double& x, double& sum);
        void qtrapwbl(double& olds, long& t, const long& f,double &p1, double &p2, const double& b, const double& c, const double& ksatl,
                double &s, long& it, long& tnm, double& del, double& x, double& sum, double& epsx);
        void get_leafPcrit(double* el, bool vCurve, double* el_v, double &p1, double &p2, const double& pinc,long& k, double& e, double& olds,
                long& t, const long& f, const double& b, const double& c, const double& ksatl,double &s, long& it, long& tnm,
                double& del, double& x, double& sum, double& epsx, double& ksh, const double& kmin, double& pcritl);
        void get_leaf(double &p1, const double& ps, const double& pgrav, double& plow, const double& pinc, double* el, double& elow, double& ehigh,
                double& estart, double& efinish, double& e, double& p2, double& pl, double& pcritl,
                long& k, long& test, std::string& failspot);
        
        // whole plant block
        void get_compositecurve(double elayer[6][100001], double* prh, double prhizo[6][100001], double& plow,
                double& p1, double& pinc, double& elow, double er[6][100001], double& klow, double kr[6][100001],double& ehigh,
                double& khigh, double& estart, double& klower, double& p2, double& efinish, double& kupper, double& flow,
                double& pr, double kroot[6][100001], double& x, double* pd, const double& root_b, const double& root_c, double* ksatr,
                double* kminroot, double* pcritrh, double* proot,double* pstem, double& ps, double* pleaf, double& pl, double* kleaf,
                double* kstem, double* kplant, double& kminleaf, double& kminstem, double& kminplant, double& e, double& pgrav,
                const double& leaf_b, const double& leaf_c, double& ksatl, const double& stem_b, const double& stem_c, double& ksats, double* eplant,
                double& einc, double* dedp, double* dedpf, double& pcritsystem, double& ecritsystem, long& k, long& p, long* layer,
                long& test, long& total, const long& layers, const bool& refilling, std::string* layerfailure);
        void storehistory(double ter[6][100001], double er[6][100001], double kr[6][100001], double tkr[6][100001],
                double er_v[6][100001], double kr_v[6][100001], double* tes, double* es, double* es_v, double* tel,
                double* el, double* el_v, double* kminroot, long& k, long* layer, long* tlayer, long& layers,
                std::string* layerfailure, std::string* tlayerfailure);
        void gethistory(double ter[6][100001], double er[6][100001], double kr[6][100001], double tkr[6][100001], double* tes,
                double* es, double* tel, double* el, long& k, long* layer, long* tlayer, long& layers,
                std::string* layerfailure, std::string* tlayerfailure);
        void get_canopypressure(const double& ecritsystem, double ter[6][100001], double er[6][100001], double kr[6][100001],
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
                long& layers, std::string* layerfailure, std::string* tlayerfailure, std::string& failspot, const std::string& night,
                const bool& refilling);
        void update_curves(double kroot[6][100001],double* kminroot, const double& pinc, double* proot,double er[6][100001],
                double kr[6][100001], double* kstem, double& kminstem, double* pstem, double* es, double* kleaf, double& kminleaf,
                double* pleaf, double* el,long& halt, long& phigh, const long& layers);
};

#endif