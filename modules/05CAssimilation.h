#ifndef CASSIMILATION
#define CASSIMILATION

// Necessary Cpp libraries
#include <cmath> // math utility functions
#include <iostream>

// Define class
class cassimilation
{
    public:
    // irradiance calculations
    void get_solarcalc(double& fet, double& et, double& sm, const double& longitude, double& tsn,
                        const double& pi, double& tsncorr, double& sindec, double& dec, double& cosdec,
                        double& tim, double& tod, double& coszen, const double& lat, double& zen,
                        double& cosaz, double& az, double& patm, double& m, const double& solar,
                        double& tau, double& sp, double& sb, double& sd, double& st, double& cloud,
                        double& obssolar, double& fcd, const double& xang, double& kbe, double& kbezero,
                        double& mleafang, double& rad, double& sum, double& k1, double& t1, double& lai,
                        double& told, double& t2, double& kd, double& qd, double& qds, const double& abspar,
                        double& qdt, double& qb, double& qbt, double& qsc, double& qsh, double& qsl,
                        double& laisl, double& laish, double& parsh, double& parsl, double& parbottom,
                        const double& absnir, double& nirsh, double& nirsl, double& sshade, double& ssun,
                        double& sbottom, double& ssunb, double& ssund, double& sref, const double& absolar,
                        double& par, double& ppfd, double& ea, double& eac, double& la, double& lg,
                        double& vpd, double& airtemp, double& maxvpd, const double& sbc,
                        long& jd);
};

#endif