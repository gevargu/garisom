#include "06MiscFunctions.h"

// Calculate atmospheric pressure using altitude
// T = 15 C, average sealevel patm; approximation
double MiscFunctions::get_patm(double alt)
{
    return 101.325 * pow((1 - 0.0065 * alt / (288.15 + 0.0065 * alt)), 5.257);
}

// Convert conductance to resistance
double MiscFunctions::get_resist(double conductance)
{
    return 1.0 / conductance;
}