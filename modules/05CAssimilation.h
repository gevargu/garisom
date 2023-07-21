#ifndef CASSIMILATION
#define CASSIMILATION

// Necessary Cpp libraries
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

// Define class
class cassimilation
{
    public:
    // irradiance calculations
    void get_solarcalc(double &fet, double &et, const double &pi, double &sm, double &longitude, double &tsn, double &tsncorr,
    double &sindec, double &dec, double &cosdec, double &tim, double &tod, double &coszen, double &lat, double &zen, double &cosaz, double &az,
    double &m, double &patm, double &sp, const double &solar, double &tau, double &sb, double &sd, double &st, double &cloud, double &obssolar,
    double &fcd, double &xang, double &kbe, double &kbezero, double &mleafang, double &rad, double &sum, double &k1, double &t1, double &lai,
    double &told, double &t2, double &kd, double &qd, double &qds, const double &abspar, double &qdt, double &qb, double &qbt, double &qsc,
    double &qsh, double &qsl, double &laisl, double &laish, double &parsh, double &parsl, double &parbottom, const double &absnir,
    const double &absolar, double &nirsh, double &nirsl, double &sshade, double &ssun, double &sbottom, double &ssunb, double &ssund,
    double &sref, double &par, double &ppfd, double &ea, double &maxvpd, double &vpd, double &eac, double &airtemp, double &la, const double &sbc,
    double &lg,long &jd);

    // leaf energy balance
    // sun layer
    void get_leaftemps(double& rabs, const double& ssun, const double& sref, const double& emiss,
        const double& la, const double& lg, double& lambda, const double& airtemp, double& grad,
        double& gha, const double& wind, const double& leafwidth, const double& laperba, double* eplantl,
        double* eplant, double& numerator, double& denominator, const double& sbc, const double& sha,
        double* leaftemp, double* lavpd, const double& patm, const double& vpd,
        long& p);
    // virgin
    void get_leaftempsmd(double& rabs, const double& ssun, const double& sref, const double& emiss, const double& la,
        const double& lg, double& lambda, const double& airtemp, double& grad, double& gha, const double& wind,
        const double& leafwidth, double& emd, const double& e, const double& laperba, const double& sbc,
        double& numerator, double& denominator, const double& sha, double& leaftmd, double& lavpdmd,
        const double& patm, const double& vpd);
    
    // shade layer
    void get_leaftempsshade(double& rabs, const double& sshade, const double& sref, const double& emiss,
        const double& lg, double& numerator, const double& sbc, const double& airtemp, const double& lambda,
        double* eplantl, double& denominator, const double& sha, const double& grad, const double& gha,
        double* leaftempsh, double* lavpdsh, const double& patm, const double& vpd,
        long& p);
    
    // virgin
    void get_leaftempsshademd(double& rabs, const double& sshade, const double& sref, const double& emiss, const double& lg,
        double& emd, double& lambda, const double& airtemp, const double& sbc, double& numerator, double& denominator,
        const double& sha, const double& grad, const double& gha, double& leaftshmd, double& lavpdshmd, const double& patm,
        const double& vpd);
    
    // carbon assimilation
    // sun layer
    void get_assimilation(double* lavpd, double* gcanw, const double& gmax, double* eplantl, double* gcanc,
        const double& comp25, double& comp, const double& gas, double* leaftemp, double& numerator,
        const double& svvmax, const double& hdvmax, const double& havmax, double& denominator,
        const double& vmax25, double& vmax, const double& svjmax, const double& hdjmax, const double& hajmax,
        double& jmax, const double& jmax25, const double& kc25, const double& ko25, double& kc, double& ko,
        double& rday25, double* rday, double& ci, double& jact, const double& qmax, const double& qsl,
        const double& lightcurv, double& je, const double& oa, double& jc, double& var, const double& thetac,
        const double& ca, double* psyn, double* cin, double& marker, double& psynmax,
        long& p, const std::string& night);
    
    // virgin
    void get_assimilationmd(double& lavpdmd, double& gcanwmd, const double& gmax, double& emd, double& gcancmd, double& comp,
        const double& comp25, const double& gas, const double& leaftmd, double& numerator, const double& svvmax,
        const double& havmax, double& denominator,const double& hdvmax, const double& vmax25, double& vmax, const double& svjmax,
        const double& hdjmax, const double& hajmax, double& jmax, const double& jmax25, const double& kc25, const double& ko25,
        double& kc, double& ko, double& rday25, double& rdaymd, double& ci, double& jact, const double& qmax, const double& qsl,
        const double& lightcurv, double& je, const double& oa, double& jc, double& var, const double& thetac, const double& ca,
        double* psynmd, double& cinmd, double& marker, double& psynmaxmd, long& p, const std::string& night);
    
    // shade layer
    void get_assimilationshade(double* lavpdsh, double* gcanwsh, const double& gmax, double* eplantl, double* gcancsh,
        double& comp, const double& comp25, const double& gas, double* leaftempsh, double& numerator, const double& svvmax,
        const double& hdvmax, const double& havmax, double& denominator, const double& vmax25, double& vmax, const double& svjmax,
        const double& hdjmax, const double& hajmax, double& jmax, const double& jmax25, const double& kc25, const double& ko25,
        double& kc, double& ko, double& rday25, double* rdaysh, double& ci, double& jact, const double& qmax, const double& qsh,
        const double& lightcurv, double& je, const double& oa, double& jc, double& var, const double& thetac, const double& ca,
        double* psynsh, double* cinsh, double& marker, double& psynmaxsh, long& p, const std::string& night);
    
    // virgin
    void get_assimilationshademd(double& lavpdshmd, double& gcanwshmd, const double& gmax, double& emd, double& gcancshmd,
        double& comp, const double& comp25, const double& gas, double& leaftshmd, double& numerator, const double& svvmax,
        const double& hdvmax, const double& havmax, double& denominator, const double& vmax25, double& vmax, const double& svjmax,
        const double& hdjmax, const double& hajmax, double& jmax, const double& jmax25, const double& kc25, const double& ko25,
        double& kc, double& ko, double& rday25, double& rdayshmd, double& ci, double& jact, const double& qmax, const double& qsh,
        const double& lightcurv, double& je, const double& oa, double& jc, double& var, const double& thetac, const double& ca,
        double* psynshmd, double& cinshmd, double& marker, double& psynmaxshmd, long& p, const std::string& night);
};

#endif