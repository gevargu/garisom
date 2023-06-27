#include "05CAssimilation.h"

/* Solar irradiance function*/
//'gets radiative terms for energy balance and assimilation
void cassimilation::get_solarcalc(double& fet, double& et, double& sm, const double& longitude, double& tsn,
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
                                  long& jd) {
    fet = 279.575 + 0.9856 * jd; //'jd is julian day, fet is factor for CN eqn 11.4
    fet = fet * pi / 180.0; //'convert to radians
    et = (-104.7 * sin(fet) + 596.2 * sin(2 * fet) + 4.3 * sin(3 * fet) - 12.7 * sin(4 * fet) - 429.3 * cos(fet) - 2 * cos(2 * fet) + 19.3 * cos(3 * fet)) / 3600.0; //'"equation of time" in fraction of HOURS C&N 11.4
    sm = 15 * int(longitude / 15.0); //'standard meridian east of longitude
    //'lc = 0.0666667 * (sm - longitude) //'CN 11.3 for longitude correction in fraction of hours
    tsn = 12 - tsncorr - et; //'time of solar noon in hour fraction from midnight, CN 11.3
    sindec = pi / 180.0 * (356.6 + 0.9856 * jd);
    sindec = sin(sindec);
    sindec = pi / 180.0 * (278.97 + 0.9856 * jd + 1.9165 * sindec);
    sindec = 0.39785 * sin(sindec); //'sine of the angle of solar declination
    dec = atan(sindec / pow((-sindec * sindec + 1), 0.5)); //'arcsin of sindec in radians, dec is solar declination
    cosdec = cos(dec);
    tim = 15 * (tod - tsn);
    tim = pi / 180.0 * tim; //'convert to radians
    coszen = sin(lat) * sindec + cos(lat) * cosdec * cos(tim); //'cos of zenith angle of sun from overhead, CN11.1
    zen = atan(-coszen / pow((-coszen * coszen + 1), 0.5)) + 2 * atan(1); //'zenith in radians
    cosaz = -(sindec - cos(zen) * sin(lat)) / (cos(lat) * sin(zen)); //'cos of azimuth angle measured counterclockwize from due south CN 11.5
    if (cosaz < -1){
        cosaz = -1; //'keep it in limits
    }
    if (cosaz > 1){
        cosaz = 1; //'ditto
    }
    if (cosaz == 1 || cosaz == -1) { //'keeps stupid acos eqn from crashing
        if (cosaz == 1){
            az = 0;
        }
        if (cosaz == -1){
            az = 3.14159; //'180 in radians
        }
    } else {
        az = atan(-cosaz / pow((-cosaz * cosaz + 1), 0.5)) + 2 * atan(1); //'solar az in radians
    } //Endif//
    
    if (tod > tsn){//'correct for afternoon hours!
        az = 6.28319 - az;
    }
    //'dayl = (-Sin(lat) * sindec) / (Cos(lat) * Cos(dec))
    //'dayl = Atn(-dayl / Sqr(-dayl * dayl + 1)) + 2 * Atn(1)
    //'dayl = 2 * dayl / 15
    if (zen * 180 / pi < 90) { //'sun//'s up: calculate solar radiation
        //'night = "n" //'its officially day
        m = patm / (101.3 * cos(zen)); //'CN 11.12
        sp = solar * pow(tau, m); //'direct beam irradiance Wm-2
        sb = sp * cos(zen); //'direct beam irradiance on horizontal surface
        sd = 0.3 * (1 - pow(tau, m)) * solar * cos(zen); //'clear sky diffuse radiation
        st = sd + sb; //'total horizontal irradiance from sun (w/o reflected radiation)
        cloud = solar * pow(0.4, m) * cos(zen); //'overcast threshold
        if (obssolar > 0) { //'we//'ve got solar data
            if (obssolar < st) { //'we//'ve got clouds
                if (obssolar > cloud) { //'we//'ve got partial clouds
                    fcd = 1 - (obssolar - cloud) / (st - cloud); //'fraction for converting beam to diffuse
                    sd = sd / st + fcd * sb / st; //'diffuse/total rato
                    sd = sd * obssolar; //'multiply ratio by total observed to get total diffuse
                    sb = obssolar - sd; //'leftover beam
                    st = obssolar; //'reset to stobs
                } else { //'its all clouds
                    sd = obssolar;
                    sb = 0;
                    st = obssolar;
                } //Endif//
            } //Endif// //'if no clouds, everything//'s already set
        } //Endif// //'if no solar data, we assume no clouds
        //'calculate reflected light as if it is equal to light at bottom of canopy
        //'if zen < 1.57 { //'sun//'s up:
        double evenMoreTest = tan(zen);
        evenMoreTest = pow((tan(zen)), 2.0);
        evenMoreTest = pow(xang, 2.0);
        evenMoreTest = pow((pow(xang, 2.0) + pow((tan(zen)), 2.0)), 2.0);
        double testNum = pow((pow(xang, 2.0) + pow((tan(zen)), 2.0)), 2.0);
        double testDenom = (xang + 1.774 * pow((xang + 1.182), -0.733));
        kbe = pow((pow(xang, 2.0) + pow((tan(zen)), 2.0)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733)); //'beam extinction coefficient CN 15.4
                                                                                                                 //kbe = Sqr(xang ^ 2 + (Tan(zen)) ^ 2) / (xang + 1.774 * (xang + 1.182) ^ -0.733) 'beam extinction coefficient CN 15.4
        kbezero = xang / (xang + 1.774 * pow((xang + 1.182), -0.733)); //'beam extinction for zen=0(overhead
        mleafang = atan(-kbezero / pow((-kbezero * kbezero + 1), 0.5)) + 2 * atan(1); //'mean leaf angle in radians
        //'gdbeamf = Exp(-Sqr(abssol) * kbe * lai) //'CN 15.6, fraction of solar beam radiation reaching ground
        //'now solve for kd (diffuse extinction) by integrating beam over all possible zenith angles from 0 to 90, CN 15.5
        rad = 0;
        sum = 0;
        k1 = pow((pow(xang, 2) + pow((tan(rad)), 2)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733));  //'beam extinction coefficient CN 15.4
        t1 = exp(-k1 * lai); //'transmittance CN 15.1
        told = t1 * sin(rad) * cos(rad); //'integral function
        do {
            rad = rad + 0.015708; //'0.9 degree intervals, zenith angle in radians
            k1 = pow((pow(xang, 2) + pow((tan(rad)), 2)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733));  //'beam extinction coefficient CN 15.4
            t1 = exp(-k1 * lai); //'transmittance CN 15.1
            t2 = t1 * sin(rad) * cos(rad); //'integral function
            sum = sum + (t2 + told) / 2.0 * 0.015708; //'integral sum
            told = t2; //'reset
        } while (!(rad > 1.5708)); //Loop Until rad > 1.5708 //'loop until 90 degrees
        sum = sum * 2; //'complete summing
        kd = -log(sum) / lai; //'extinction coefficient for diffuse radiation
        //'now...compute shaded leaf ppfd, q denotes a ppfd
        qd = 0.45 * sd * 4.6; //'converts total solar diffuse to PPFD diffuse in umol m-2 s-1
        qds = qd * (1 - exp(-(pow(abspar, 0.5) * kd * lai))) / (pow(abspar, 0.5) * kd * lai); //'mean diffuse irradiance for shaded leaves CN p. 261
        qdt = qd * exp(-(pow(abspar, 0.5) * kd * lai)); //'diffuse irradiance at bottom of canopy, CN p. 255, Eqn 15.6
        qb = 0.45 * sb * 4.6; //'converts total solar beam to PPFD
        qbt = qb * exp(-(pow(abspar, 0.5) * kbe * lai)); //'direct AND downscattered PPFD at bottom of canopy, CN 15.6
        qb = qb * exp(-(kbe * lai)); //'direct PPFD at bottom of canopy,CN 15.1
        qsc = (qbt - qb) / 2.0; //'average backscattered beam on shaded leaves
        qsh = qds + qsc; //'average PPFD incident (not absorbed!) on shaded leaves
        qb = 0.45 * sb * 4.6; //'re-set qb to top of canopy
        //'now get sunlit leaf ppfd
        qsl = kbe * qb + qsh; //'average PPFD incident (not absorbed) on sunlit leaves
        //'now get sunlit vs. shaded lai
        laisl = (1 - exp(-(kbe * lai))) / kbe; //'sunlit lai
        laish = lai - laisl; //'shaded lai
        //'now get sun and shade PAR, NIR, & longwave in preparation for energy balance
        //'in par range
        parsh = qsh / 4.6; //'incident PAR, Wm-2, shaded leaves, 100% diffuse
        parsl = qsl / 4.6;//'incident PAR, sunlit leaves
        parbottom = (qbt + qdt) / 4.6; //'PAR making it through the canopy
        //'in near-infrared range (assume same equations as for PAR, but different absorptances and incoming fluxes in Wm-2
        qd = 0.55 * sd; //'converts total solar diffuse to NIR diffuse in Wm-2
        qds = qd * (1 - exp(-(pow(absnir, 0.5) * kd * lai))) / (pow(absnir, 0.5) * kd * lai); //'mean diffuse nir irradiance for shaded leaves CN p. 261
        qdt = qd * exp(-(pow(absnir, 0.5) * kd * lai)); //'diffuse NIR irradiance at bottom of canopy, CN p. 255, Eqn 15.6
        qb = 0.55 * sb; //'converts total solar beam to NIR
        //'qb = 1600
        qbt = qb * exp(-(pow(absnir, 0.5) * kbe * lai)); //'direct AND downscattered NIR at bottom of canopy, CN 15.6
        qb = qb * exp(-(kbe * lai)); //'direct NIR at bottom of canopy,CN 15.1, using same extinction coefficient as for PAR
        qsc = (qbt - qb) / 2.0; //'average backscattered beam on shaded leaves
        nirsh = (qds + qsc); //'incident NIR on shaded leaves, 100% diffuse
        qb = 0.55 * sb; //'re-set nir qb to top of canopy to get sunlit leaves
        nirsl = kbe * qb + nirsh; //'average incident NIR on sunlit leaves
        sshade = parsh + nirsh; //'total solar incident on shaded leaves, 100% diffuse
        ssun = parsl + nirsl; //'total solar incident on sunlit leaves
        sbottom = parbottom + qdt + qbt; //'total solar at bottom of canopy
        ssunb = sb / st * ssun; //'beam solar on sunlit (an approximation)
        ssund = sd / st * ssun;//'diffuse solar on sunlit (approximation)
        //'abssolar = Cells(8, 17) //'absorptivity of leaves for total solar (0.5)
        sref = (1 - absolar) * sshade; //'reflected light...this used for sun/shade dichotomy
        sref = (1 - absolar) * sbottom; //'reflected light for monolayer version
        //'these below are used for monolayer version:
        par = 0.45 * st; //'wm-2 in par wavelength...45% of total solar
        ppfd = par * 4.6; //'assumes 4.6 moles photons per Joule conversion factor
    } else { //'sun//'s down
        sp = 0; sb = 0; sd = 0; st = 0; sref = 0; par = 0; ppfd = 0; sref = 0; ssun = 0; sshade = 0;
        qsh = 0; qsl = 0; ssunb = 0; ssund = 0; laisl = 0; laish = 0; sbottom = 0; //'sun//'s down
        //'night = "y" //'it//'s officially night
    } //Endif//
    //'now compute long wave irradiance
    ea = patm * (maxvpd - vpd); //'vapor pressure in kPa
    eac = 1.72 * pow((ea / (airtemp + 273.15)), (1.0 / 7.0)); //'emissivity of clear sky CN  10.10
    la = eac * sbc * pow((airtemp + 273.15), 4); //'long wave irradiance from clear sky
    lg = 0.97 * sbc * pow((airtemp + 273.15), 4); //'long wave irradiance from ground...assumes equilibrium with air temp
}

/* Leaf energy balance functions*/
// gets sun layer leaf temp and leaf-to-air vpd from E, Campbell and Norman book
void cassimilation::get_leaftemps(double& rabs, const double& ssun, const double& sref, const double& emiss,
    const double& la, const double& lg, double& lambda, const double& airtemp, double& grad,
    double& gha, const double& wind, const double& leafwidth, const double& laperba, double* eplantl,
    double* eplant, double& numerator, double& denominator, const double& sbc, const double& sha,
    double* leaftemp, double* lavpd, const double& patm, const double& vpd,
    long& p){
    rabs = 0.5 * (0.5 * ssun + 0.5 * sref) + emiss * (0.5 * la + 0.5 * lg); //'total absorbed radiation for sun leaves; CN 11.14
    lambda = -42.9143 * airtemp + 45064.3; //'heat of vaporization for water at air temp in J mol-1
    grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * pow(airtemp, 2); //'radiative conductance (long wave) at air temp in mol m-2 s-1
    gha = 1.4 * 0.135 * pow((wind / leafwidth), 0.5); //'heat conductance in mol m-2s-1
    eplantl[p] = eplant[p] * (1.0 / laperba) * (1.0 / 3600.0) * 55.4; //'convert to E per leaf area in mol m-2s-1
    numerator = rabs - emiss * sbc * pow((airtemp + 273.2), 4) - lambda * eplantl[p] / 2.0; //'divide E by 2 because energy balance is two sided.
    denominator = sha * (grad + gha);
    leaftemp[p] = airtemp + numerator / denominator; //'leaf temp for supply function
    lavpd[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    lavpd[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftemp[p])) - lavpd[p] + vpd; //'leaf-to-air vpd
    if (lavpd[p] < 0){
        lavpd[p] = 0; //'don//'t allow negative lavpd
    }               //'eplantl[p] = eplantl[p] * 1000 //'convert to mmol m-2 s-1
}

//'gets virgin sun layer leaf temp and leaf-to-air vpd from E, Campbell and Norman book
void cassimilation::get_leaftempsmd(double& rabs, const double& ssun, const double& sref, const double& emiss, const double& la,
    const double& lg, double& lambda, const double& airtemp, double& grad, double& gha, const double& wind,
    const double& leafwidth, double& emd, const double& e, const double& laperba, const double& sbc,
    double& numerator, double& denominator, const double& sha, double& leaftmd, double& lavpdmd,
    const double& patm, const double& vpd){
    rabs = 0.5 * (0.5 * ssun + 0.5 * sref) + emiss * (0.5 * la + 0.5 * lg); //'total absorbed radiation for sun leaves; CN 11.14
    lambda = -42.9143 * airtemp + 45064.3; //'heat of vaporization for water at air temp in J mol-1
    grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * pow(airtemp, 2); //'radiative conductance (long wave) at air temp in mol m-2 s-1
    gha = 1.4 * 0.135 * pow((wind / leafwidth), 0.5); //'heat conductance in mol m-2s-1
    emd = e * (1 / laperba) * (1.0 / 3600.0) * 55.4; //'convert to E per leaf area in mol m-2s-1
    numerator = rabs - emiss * sbc * pow((airtemp + 273.2), 4) - lambda * emd / 2.0; //'divide E by 2 because energy balance is two sided.
    denominator = sha * (grad + gha);
    leaftmd = airtemp + numerator / denominator; //'leaf temp for supply function
    lavpdmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    lavpdmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftmd)) - lavpdmd + vpd; //'leaf-to-air vpd
    if (lavpdmd < 0){
        lavpdmd = 0;//'don't allow negative lavpd
    }
    //'eplantl(p) = eplantl(p) * 1000 //'convert to mmol m-2 s-1
}

// gets shade layer leaf temp and leaf-to-air vpd from E, Campbell and Norman book   
void cassimilation::get_leaftempsshade(double& rabs, const double& sshade, const double& sref, const double& emiss,
    const double& lg, double& numerator, const double& sbc, const double& airtemp, const double& lambda,
    double* eplantl, double& denominator, const double& sha, const double& grad, const double& gha,
    double* leaftempsh, double* lavpdsh, const double& patm, const double& vpd,
    long& p){
    rabs = 0.5 * (0.5 * sshade + 0.5 * sref) + emiss * lg; //'total absorbed radiation for shaded leaves
    //'lambda = -42.9143 * airtemp + 45064.3 //'heat of vaporization for water at air temp in J mol-1
    //'grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * airtemp ^ 2 //'radiative conductance (long wave) at air temp in mol m-2 s-1
    //'gha = 1.4 * 0.135 * (wind / leafwidth) ^ 0.5 //'heat conductance in mol m-2s-1
    //'eplantl[p] = eplant[p] * (1 / laperba) * (1 / 3600) * 55.4 //'convert to E per leaf area in mol m-2s-1
    numerator = rabs - emiss * sbc * pow((airtemp + 273.2), 4) - lambda * eplantl[p] / 2.0; //'divide E by 2 because energy balance is two sided.
    denominator = sha * (grad + gha);
    leaftempsh[p] = airtemp + numerator / denominator; //'leaf temp for supply function
    lavpdsh[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    lavpdsh[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftempsh[p])) - lavpdsh[p] + vpd; //'leaf-to-air vpd
    if (lavpdsh[p] < 0){
        lavpdsh[p] = 0; //'don//'t allow negative lavpd
    }
    eplantl[p] = eplantl[p] * 1000; //'convert to mmol m-2 s-1
}

//'gets virgin shade layer leaf temp and leaf-to-air vpd from E, Campbell and Norman book   
void cassimilation::get_leaftempsshademd(double& rabs, const double& sshade, const double& sref, const double& emiss, const double& lg,
    double& emd, double& lambda, const double& airtemp, const double& sbc, double& numerator, double& denominator,
    const double& sha, const double& grad, const double& gha, double& leaftshmd, double& lavpdshmd, const double& patm, const double& vpd){
    rabs = 0.5 * (0.5 * sshade + 0.5 * sref) + emiss * lg; //'total absorbed radiation for shaded leaves
    numerator = rabs - emiss * sbc * pow((airtemp + 273.2), 4) - lambda * emd / 2.0; //'divide E by 2 because energy balance is two sided.
    denominator = sha * (grad + gha);
    leaftshmd = airtemp + numerator / denominator; //'leaf temp for supply function
    lavpdshmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
    lavpdshmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftshmd)) - lavpdshmd + vpd; //'leaf-to-air vpd
    if (lavpdshmd < 0){
        lavpdshmd = 0; //'don't allow negative lavpd
    }
    emd = emd * 1000; //'convert to mmol m-2 s-1
}

/* Photosynthesis fcuntions*/
// sun layer photosynthesis
void cassimilation::get_assimilation(double* lavpd, double* gcanw, const double& gmax, double* eplantl, double* gcanc,
    const double& comp25, double& comp, const double& gas, double* leaftemp, double& numerator,
    const double& svvmax, const double& hdvmax, const double& havmax, double& denominator,
    const double& vmax25, double& vmax, const double& svjmax, const double& hdjmax, const double& hajmax,
    double& jmax, const double& jmax25, const double& kc25, const double& ko25, double& kc, double& ko,
    double& rday25, double* rday, double& ci, double& jact, const double& qmax, const double& qsl,
    const double& lightcurv, double& je, const double& oa, double& jc, double& var, const double& thetac,
    const double& ca, double* psyn, double* cin, double& marker, double& psynmax,
    long& p, const std::string& night){
    //'get g from D and e
    if (lavpd[p] == 0) { //if//
        gcanw[p] = gmax; //'maximum g if no lavpd
    } else {
        gcanw[p] = eplantl[p] / lavpd[p]; //'gcanopy in mmol m-2s-1 (leaf area)
    } //End if//
    //'gcanw[p] = eplantl[p] / lavpd[p] //'gcanopy in mmol m-2s-1 (leaf area)
    gcanc[p] = (gcanw[p] / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
    //'adjust photosynthetic inputs for Tleaf
    comp = comp25 * exp((37830 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
    numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (gas * 298.2))) * exp((havmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftemp[p])));
    denominator = 1 + exp((svvmax * (273.2 + leaftemp[p]) - hdvmax) / (gas * (273.2 + leaftemp[p])));
    vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
    numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (gas * 298.2))) * exp((hajmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftemp[p])));
    denominator = 1 + exp((svjmax * (273.2 + leaftemp[p]) - hdjmax) / (gas * (273.2 + leaftemp[p])));
    jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
    kc = kc25 * exp((79430 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
    ko = ko25 * exp((36380 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
    rday25 = vmax25 * 0.01; //'from Medlyn 2002
    rday[p] = rday25 * pow(2, ((leaftemp[p] - 25) / 10.0));
    rday[p] = rday[p] * pow((1 + exp(1.3 * (leaftemp[p] - 55))), -1); //'high temp inhibition, collatz
    if (night == "n" && gcanc[p] > 0) { //if// //'solve for A and ci
        if (p == 1) { //if// //'stomata have just opened: find the mitochondrial compensation point
            ci = comp - 0.00000001; //'start at photorespiratory compensation point
            do{//'loop through A-ci curve to find mitochondrial compensation point
                ci = ci + 0.00000001;
                jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + oa / ko);
                jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rday[p]; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));//Loop Until var >= 0 Or ci >= ca //'there//'s positive Anet or not enough light
            psyn[p] = var; //'always start with zero or above
            cin[p] = ci; //'the dark compensation point...not predicting negative Anet
        } //End if// //'p=1 if
        
        if (p > 1) { //if//
            ci = cin[p - 1]; //'start from previous ci
            do{//'loop through A-ci curve to find ci
                ci = ci + 0.00000001;
                marker = gcanc[p] * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
                jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + oa / ko);
                jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rday[p]; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));//Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that//'s the right ci...
            psyn[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psyn[p] > psynmax){
                psynmax = psyn[p];
            }
            cin[p] = ci; //'store ci
        } //End if// //'p>1 if
    } else { //'it//'s night or stomata are closed
        if (night == "y"){
            psyn[p] = 0 - rday[p];
            psynmax = 0;
            cin[p] = ca; //'respiration accounted for at night
        }
        if (night == "n") {
            psyn[p] = 0;
            psynmax = 0;
            cin[p] = ca; //'it//'s day and p=0, g=0
        }
    } //End if// //'night if
}

//'gets virgin assimilation for sun leaves
void cassimilation::get_assimilationmd(double& lavpdmd, double& gcanwmd, const double& gmax, double& emd, double& gcancmd, double& comp,
    const double& comp25, const double& gas, const double& leaftmd, double& numerator, const double& svvmax,
    const double& havmax, double& denominator,const double& hdvmax, const double& vmax25, double& vmax, const double& svjmax,
    const double& hdjmax, const double& hajmax, double& jmax, const double& jmax25, const double& kc25, const double& ko25,
    double& kc, double& ko, double& rday25, double& rdaymd, double& ci, double& jact, const double& qmax, const double& qsl,
    const double& lightcurv, double& je, const double& oa, double& jc, double& var, const double& thetac, const double& ca,
    double* psynmd, double& cinmd, double& marker, double& psynmaxmd, long& p, const std::string& night){
    //'get g from D and e
    if (lavpdmd == 0){
        gcanwmd = gmax; //'maximum g if no lavpd
    } else{
        gcanwmd = emd / lavpdmd; //'gcanopy in mmol m-2s-1 (leaf area)
    }
    //'gcanw(p) = eplantl(p) / lavpd(p) //'gcanopy in mmol m-2s-1 (leaf area)
    gcancmd = (gcanwmd / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
    //'adjust photosynthetic inputs for Tleaf
    comp = comp25 * exp((37830 * ((leaftmd + 273.15) - 298.15)) / (298.15 * gas * (leaftmd + 273.15))); //'Bernacchi via Medlyn
    numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (gas * 298.2))) * exp((havmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftmd)));
    denominator = 1 + exp((svvmax * (273.2 + leaftmd) - hdvmax) / (gas * (273.2 + leaftmd)));
    vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
    numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (gas * 298.2))) * exp((hajmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftmd)));
    denominator = 1 + exp((svjmax * (273.2 + leaftmd) - hdjmax) / (gas * (273.2 + leaftmd)));
    jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
    kc = kc25 * exp((79430 * ((leaftmd + 273.15) - 298.15)) / (298.15 * gas * (leaftmd + 273.15))); //'Bernacchi via Medlyn
    ko = ko25 * exp((36380 * ((leaftmd + 273.15) - 298.15)) / (298.15 * gas * (leaftmd + 273.15)));//'Bernacchi via Medlyn
    rday25 = vmax25 * 0.01; //'from Medlyn 2002
    rdaymd = rday25 * pow(2, ((leaftmd - 25) / 10.0));
    rdaymd = rdaymd * pow((1 + exp(1.3 * (leaftmd - 55))), -1); //'high temp inhibition, collatz
    if (night == "n" && gcancmd > 0){ //'solve for A and ci
        if (p == 1){//'stomata have just opened: find the mitochondrial compensation point
            ci = comp - 0.00000001; //'start at photorespiratory compensation point
            do{//'loop through A-ci curve to find mitochondrial compensation point
                ci = ci + 0.00000001;
                jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + oa / ko);
                jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdaymd; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there's positive Anet or not enough light
            psynmd[p] = var; //'always start with zero or above
            cinmd = ci; //'the dark compensation point...not predicting negative Anet
        } //'p=1 if
        if (p > 1){
            ci = cinmd - 0.00000001; //'cin(p - 1) //'start from previous ci
            do{//'loop through A-ci curve to find ci
                ci = ci + 0.00000001;
                marker = gcancmd * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
                jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + oa / ko);
                jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdaymd; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that's the right ci...
            psynmd[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynmd[p] > psynmaxmd){
                psynmaxmd = psynmd[p];
            }
            cinmd = ci; //'store ci
        } //'p>1 if
    } else {//'it's night or stomata are closed
        if (night == "y"){
            psynmd[p] = 0 - rdaymd;
            psynmaxmd = 0; //'respiration accounted for at night
        }
        if (night == "n"){
            psynmd[p] = 0;
            psynmaxmd = 0; //'it's day and p=0, g=0
        }
    } //'night if
    //'Cells(16 + p, 65) = psynsh(p)
}

// shade layer photosynthesis
void cassimilation::get_assimilationshade(double* lavpdsh, double* gcanwsh, const double& gmax, double* eplantl, double* gcancsh,
    double& comp, const double& comp25, const double& gas, double* leaftempsh, double& numerator, const double& svvmax,
    const double& hdvmax, const double& havmax, double& denominator, const double& vmax25, double& vmax, const double& svjmax,
    const double& hdjmax, const double& hajmax, double& jmax, const double& jmax25, const double& kc25, const double& ko25,
    double& kc, double& ko, double& rday25, double* rdaysh, double& ci, double& jact, const double& qmax,  const double& qsh,
    const double& lightcurv, double& je, const double& oa, double& jc, double& var, const double& thetac, const double& ca,
    double* psynsh, double* cinsh, double& marker, double& psynmaxsh, long& p, const std::string& night){
    //'get g from D and e
    if (lavpdsh[p] == 0){ //if//
        gcanwsh[p] = gmax; //'set to maximum if no vpd
    } else {
        gcanwsh[p] = eplantl[p] / lavpdsh[p]; //'gcanopy in mmol m-2s-1 (leaf area)
    } //End if//
    //'gcanwsh[p] = eplantl[p] / lavpdsh[p] //'gcanopy in mmol m-2s-1 (leaf area)
    gcancsh[p] = (gcanwsh[p] / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
    //'adjust photosynthetic inputs for Tleaf
    comp = comp25 * exp((37830 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
    numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (gas * 298.2))) * exp((havmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftempsh[p])));
    denominator = 1 + exp((svvmax * (273.2 + leaftempsh[p]) - hdvmax) / (gas * (273.2 + leaftempsh[p])));
    vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
    numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (gas * 298.2))) * exp((hajmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftempsh[p])));
    denominator = 1 + exp((svjmax * (273.2 + leaftempsh[p]) - hdjmax) / (gas * (273.2 + leaftempsh[p])));
    jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
    kc = kc25 * exp((79430 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
    ko = ko25 * exp((36380 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
    rday25 = vmax25 * 0.01; //'from Medlyn 2002
    rdaysh[p] = rday25 * pow(2, ((leaftempsh[p] - 25) / 10.0));
    rdaysh[p] = rdaysh[p] * pow((1 + exp(1.3 * (leaftempsh[p] - 55))), -1); //'high temp inhibition, collatz
    if (night == "n" && gcancsh[p] > 0) { //if// //'solve for A and ci
        if (p == 1) { //if// //'first find the mitochondrial compensation point
            ci = comp - 0.00000001; //'start at photorespiratory compensation point
            do{//'loop through A-ci curve to find mitochondrial compensation point
                ci = ci + 0.00000001;
                jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + oa / ko);
                jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdaysh[p]; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there//'s positive Anet
            psynsh[p] = var; //'always start with zero or above
            cinsh[p] = ci; //'the dark compensation point...not predicting negative Anet
        } //End if// //'p=0 if
        if (p > 1) { //if//
            ci = cinsh[p - 1]; //'start from previous ci
            do{//'loop through A-ci curve to find mitochondrial compensation point
                ci = ci + 0.00000001;
                marker = gcancsh[p] * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
                jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + oa / ko);
                jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdaysh[p]; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that//'s the right ci...
            psynsh[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynsh[p] > psynmaxsh){
                psynmaxsh = psynsh[p];
            }
            cinsh[p] = ci; //'store ci
        } //End if// //'p>1 if
    } else { //'it//'s night
        if (night == "y") {
            psynsh[p] = 0 - rdaysh[p];
            psynmaxsh = 0;
            cinsh[p] = ca; //'respiration accounted for at night
        }
        if (night == "n") {
            psynsh[p] = 0;
            psynmaxsh = 0;
            cinsh[p] = ca; //'it//'s day and p=0, g=0
        }
    } //End if// //'night if
}

// gets virgin assimilation for shade leaves for midday solution
void cassimilation::get_assimilationshademd(double& lavpdshmd, double& gcanwshmd, const double& gmax, double& emd, double& gcancshmd,
    double& comp, const double& comp25, const double& gas, double& leaftshmd, double& numerator, const double& svvmax,
    const double& hdvmax, const double& havmax, double& denominator, const double& vmax25, double& vmax, const double& svjmax,
    const double& hdjmax, const double& hajmax, double& jmax, const double& jmax25, const double& kc25, const double& ko25,
    double& kc, double& ko, double& rday25, double& rdayshmd, double& ci, double& jact, const double& qmax,  const double& qsh,
    const double& lightcurv, double& je, const double& oa, double& jc, double& var, const double& thetac, const double& ca,
    double* psynshmd, double& cinshmd, double& marker, double& psynmaxshmd, long& p, const std::string& night){//'get g from D and e
    if (lavpdshmd == 0){
        gcanwshmd = gmax; //'set to maximum if no vpd
    } else {
        gcanwshmd = emd / lavpdshmd; //'gcanopy in mmol m-2s-1 (leaf area)
    }
    //'gcanwsh(p) = eplantl(p) / lavpdsh(p) //'gcanopy in mmol m-2s-1 (leaf area)
    gcancshmd = (gcanwshmd / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
    //'adjust photosynthetic inputs for Tleaf
    comp = comp25 * exp((37830 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * gas * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
    numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (gas * 298.2))) * exp((havmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftshmd)));
    denominator = 1 + exp((svvmax * (273.2 + leaftshmd) - hdvmax) / (gas * (273.2 + leaftshmd)));
    vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
    numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (gas * 298.2))) * exp((hajmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftshmd)));
    denominator = 1 + exp((svjmax * (273.2 + leaftshmd) - hdjmax) / (gas * (273.2 + leaftshmd)));
    jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
    kc = kc25 * exp((79430 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * gas * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
    ko = ko25 * exp((36380 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * gas * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
    rday25 = vmax25 * 0.01; //'from Medlyn 2002
    rdayshmd = rday25 * pow(2, ((leaftshmd - 25) / 10.0));
    rdayshmd = rdayshmd * pow((1 + exp(1.3 * (leaftshmd - 55))), -1); //'high temp inhibition, collatz
    if (night == "n" && gcancshmd > 0){ //'solve for A and ci
        if (p == 1) { //'first find the mitochondrial compensation point
            ci = comp - 0.00000001; //'start at photorespiratory compensation point
            do{//'loop through A-ci curve to find mitochondrial compensation point
                ci = ci + 0.00000001;
                jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + oa / ko);
                jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdayshmd; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there's positive Anet
            psynshmd[p] = var; //'always start with zero or above
            cinshmd = ci; //'the dark compensation point...not predicting negative Anet
        } //'p=1 if
        if (p > 1){
            ci = cinshmd - 0.00000001; //'start from previous ci backed off a bit
            do{//'loop through A-ci curve to find mitochondrial compensation point
                ci = ci + 0.00000001;
                marker = gcancshmd * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
                jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
                je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
                numerator = vmax * (ci - comp);
                denominator = ci + kc * (1 + oa / ko);
                jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
                var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
                var = var - rdayshmd; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that's the right ci...
            psynshmd[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynshmd[p] > psynmaxshmd){
                psynmaxshmd = psynshmd[p];
            }
            cinshmd = ci; //'store ci
        } //'p>1 if
    } else {//'it's night
        if (night == "y"){
            psynshmd[p] = 0 - rdayshmd;
            psynmaxshmd = 0; //'respiration accounted for at night
        }
        if (night == "n"){
            psynshmd[p] = 0;
            psynmaxshmd = 0; //'it's day and p=0, g=0
        }
    } //'night if
    //'Cells(16 + p, 65) = psynsh(p)
}
