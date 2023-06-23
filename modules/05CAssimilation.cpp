#include "05CAssimilation.h"

/* Solar irradiance function*/
//'gets radiative terms for energy balance and assimilation
void cassimilation::get_solarcalc() {
    //'j = Cells(14 + i, 3) //'julian day
    fet = 279.575 + 0.9856 * jd; //'jd is julian day, fet is factor for CN eqn 11.4
    fet = fet * pi / 180.0; //'convert to radians
    et = (-104.7 * sin(fet) + 596.2 * sin(2 * fet) + 4.3 * sin(3 * fet) - 12.7 * sin(4 * fet) - 429.3 * cos(fet) - 2 * cos(2 * fet) + 19.3 * cos(3 * fet)) / 3600.0; //'"equation of time" in fraction of HOURS C&N 11.4
                                                                                                                                                                       //'long = Cells(11, 4) //'longitude in degree fraction W
    sm = 15 * int(longitude / 15.0); //'standard meridian east of longitude
    //'lc = 0.0666667 * (sm - longitude) //'CN 11.3 for longitude correction in fraction of hours
    tsn = 12 - tsncorr - et; //'time of solar noon in hour fraction from midnight, CN 11.3
    //'Cells(14 + i, 5) = tsn //'output solar noon
    sindec = pi / 180.0 * (356.6 + 0.9856 * jd);
    sindec = sin(sindec);
    sindec = pi / 180.0 * (278.97 + 0.9856 * jd + 1.9165 * sindec);
    sindec = 0.39785 * sin(sindec); //'sine of the angle of solar declination
    //'lat = Cells(10, 4) //'latitude, fraction of degrees N
    dec = atan(sindec / pow((-sindec * sindec + 1), 0.5)); //'arcsin of sindec in radians, dec is solar declination
    //'Cells(14 + i, 6) = dec * 180 / pi //'output solar declination
    cosdec = cos(dec);
    //'timst = Cells(14 + i, 4) //'local standard time, hour fraction from midnight
    tim = 15 * (tod - tsn);
    tim = pi / 180.0 * tim; //'convert to radians
    coszen = sin(lat) * sindec + cos(lat) * cosdec * cos(tim); //'cos of zenith angle of sun from overhead, CN11.1
    zen = atan(-coszen / pow((-coszen * coszen + 1), 0.5)) + 2 * atan(1); //'zenith in radians
    //'if zen < 1.57 { Cells(14 + i, 7) = zen * 180 / pi //'output when sun//'s up (zen<90)
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
    //'if zen < 1.57 { Cells(14 + i, 8) = az * 180 / pi//'output azimuth during day
    //'dayl = (-Sin(lat) * sindec) / (Cos(lat) * Cos(dec))
    //'dayl = Atn(-dayl / Sqr(-dayl * dayl + 1)) + 2 * Atn(1)
    //'dayl = 2 * dayl / 15
    //'Cells(14 + i, 9) = dayl * 180 / pi
    
    if (zen * 180 / pi < 90) { //'sun//'s up: calculate solar radiation
        //'pa = Cells(12, 4) //'atmospheric p in kpa
        //'Cells(16 + dd, 8) = zen * 180 / pi //'zenith angle in degrees
        //'Cells(16 + dd, 9) = lat
        //'night = "n" //'its officially day
        m = patm / (101.3 * cos(zen)); //'CN 11.12
        //'spo = Cells(9, 4) //'solar constant
        //'tau = Cells(8, 4) //'transmittance, CN 11.11
        sp = solar * pow(tau, m); //'direct beam irradiance Wm-2
        sb = sp * cos(zen); //'direct beam irradiance on horizontal surface
        //'Cells(14 + i, 11) = sb
        sd = 0.3 * (1 - pow(tau, m)) * solar * cos(zen); //'clear sky diffuse radiation
        //'Cells(14 + i, 12) = sd
        st = sd + sb; //'total horizontal irradiance from sun (w/o reflected radiation)
        cloud = solar * pow(0.4, m) * cos(zen); //'overcast threshold
        //'stobs = Cells(13, 4) //'observed solar radiation on the horizontal Wm-2
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
        //'xang = Cells(12, 11) //'leaf angle parameter, CN 15.4
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
        //'Cells(14 + i, 20) = kbe
        //'Cells(14 + i, 21) = mleafang * 180 / pi //'mean leaf angle in degrees
        //'lai = Cells(13, 11) //'canopy leaf area index
        //'abspar = Cells(7, 17) //'absorptivity for PAR of leaves
        //'abssol = Cells(8, 17) //'absorptivity for total solar of leaves
        //'gdbeamf = Exp(-Sqr(abssol) * kbe * lai) //'CN 15.6, fraction of solar beam radiation reaching ground
        //'Cells(14 + i, 22) = gdbeamf
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
        //'Cells(14 + i, 23) = kd //'output
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
        //'absnir = Cells(9, 17) //'absorptivity of leaves to NIR
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
    //'ta = Cells(3, 19) + 273.15 //'air temp in K
    eac = 1.72 * pow((ea / (airtemp + 273.15)), (1.0 / 7.0)); //'emissivity of clear sky CN  10.10
    //'boltz = Cells(9, 11) //'boltzman constant
    la = eac * sbc * pow((airtemp + 273.15), 4); //'long wave irradiance from clear sky
    lg = 0.97 * sbc * pow((airtemp + 273.15), 4); //'long wave irradiance from ground...assumes equilibrium with air temp
    //'Cells(14 + i, 16) = la
    //'Cells(14 + i, 17) = lg
}