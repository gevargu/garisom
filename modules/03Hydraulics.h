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

        // genaral
        double get_cweibull(const double& p12, const double& p50);
        double get_bweibull(const double& p12, const double& c);
        double get_weibullfit(const double& x,const double& b, const double& c, const double& ksat);
        // root block
        double get_rootpercent(const double& leafpercent);
        double get_rootkmax(const double& rootpercent,const double& rsatp);
        double get_weibullfitroot(const double& x, const double& b, const double& c, double* ksatr, const int& z);
        void trapzdwbr(const long& t, double& p1, double& p2, const double& b, const double& c,double* ksatr, const int& z, 
                double& s, long& it, long& tnm, double& del, double& x, double& sum);
        void qtrapwbr(double& olds, long& t, double& p1, double& p2, const double& b, const double& c,double* ksatr, 
                const int& z, double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx);
        void get_rootPcrit(double er[][100001],double kr[][100001],long& z, const int& layers,double* ksatr, double& p1, 
                double& p2, const double& pinc,long& k,double& e,double& olds, long& t, const double& b, const double& c, 
                double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx, 
                const double& kmin,double* pcritr);
        void get_rootPcrit_v(double er_v[][100001],double kr_v[][100001],long& z, const int& layers,double* ksatr, double& p1, 
                double& p2, const double& pinc,long& k,double& e,double& olds, long& t, const double& b, const double& c, 
                double& s, long& it, long& tnm, double& del, double& x, double& sum, const long& f, double& epsx, 
                const double& kmin,double* pcritr);
        void get_rhizoflow(double& plow, double& p1, double& pinc, double& elow, double erh[6][100001], double& klow,
                double krh[6][100001],double& ehigh,  double& khigh, double& estart, double& klower, double& p2,
                double& efinish, double& kupper, double& flow,long& k, long& z);
        void get_rootflow(double& plow, double& p1, double& pinc, double& elow, double er[6][100001], double& klow,
                double kr[6][100001],double& ehigh,  double& khigh, double& estart, double& klower, double& p2,
                double& efinish, double& kupper, double& flow, long& k, long& z);
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
                int& layers, std::string* layerfailure, std::string& failspot);
        
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
};

#endif