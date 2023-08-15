// Functions to calculate root morphological traits
#include "04Morphology.h"

/* Root morphological parameters*/

double morphology::get_rootbeta(double rootdepth)
  {
    /* Root biomass distribution is allocated based on the equation reported in
    Love et al (2019): 
    M = 1 - Beta^d, where M is the fraction of biomass above depth d
    expressed in cm. We find the Beta that provides an M of 0.995 for the
    maximum rooting depth.*/
    double d;
    double beta;
    d = rootdepth*100; // Converts bedrock rooting depth from m to cm
    beta = pow((1-0.995),(1/d)); // Calculates Beta
    return beta;
  }
