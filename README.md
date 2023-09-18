# Carbon Gain vs Hydraulic Risk Stomatal Optimization Model V 2.0

__Coding language:__ C++

__Authors:__ German Vargas G. & William R.L. Anderegg

__Contact:__ german.vargas@utah.edu

------------

## Introduction:

The model uses a stomatal carbon gain vs. hydraulic risk optimization (Sperry et al. 2017 -- see references below) combined with a soil water budget to run a continuous growing season simulation, producing stand outputs such as net carbon assimilation (Anet), internal [CO2] (Ci), transpiration (E), total evapotranspiration (ET), and element conductances (k) on an hourly and summary basis (see output details below).

------------

## Basic usage instructions:

First, fork the repository to your own account. If you don't have a github account, just download the repository. Provide plant and site traits through the parameters file, then supply hourly weather data to drive model. The model will expect these files to be located in the current working directory. See the included examples and the details below for formatting and units. These example files contain all inputs necessary to test-run the model immediately after building.

__a)__ To run, build and execute the model program (no command line arguments -- runs from files in the working directory). Building requires c++11, GNU example:

```{}
g++ -std=c++11
```

The -O3 and -ffast-math optimizations are recommended with GNU compilers:

```{}
g++ -std=c++11 -O3 -ffast-math
```

__b)__ To build and run from the terminal in a Mac OS system:

__b.1)__ Clone the repository in a desired location in your system.

```
git clone https://github.com/gevargu/garisom.git
```

__b.2)__ Navigate to the repository by usin the command cd:

```
cd garisom
```

__b.3)__ Run the following command will compile the code and build an executable:

```
make
```

__b.4)__ Run this program from the same folder with this command:

```
./garisom
```

__b.5)__ Press <kbd>Command</kbd> + <kbd>C</kbd> if you want to stop the model before it completes (<kbd>Ctrl</kbd> + <kbd>C</kbd> in Windows/Linux)

This version has also been tested with Visual Studio 2017's compiler with similar optimizations (floating point mode fast, maximum optimization preferring speed). The following files (included in this repository) should be located in the working directory (normally the same directory as the executable) before running:
	
  - __parameters_2.0.0.csv__ (site, atmospheric, soils, stand, plant, hydraulics, and carbon assimilation parameters).
  - __configuration_2.0.0.csv__ (model controls)
  - __dataset.csv__ (hourly weather drivers).
  - __dataheader.csv__ (a header row for the hourly data output).
  - __sumheader.csv__ (a header row for the summary data output).
  - __seasonlimits_2.0.0.csv__ (growing season limits and yearly atmospheric CO2, only required if using "sequential year mode" described below)

Upon completion, two output files are produced:
	
  - The __timesteps.csv__ output contains all of the hourly model outputs corresponding to the weather inputs.
  - The __summary.csv__ output includes various total values for each year (for example, net growing season productivity Anet, net transpiration E, etc.)

------------

## Plant and Stand Parameters:

Configure plant traits and other parameters in __parameters.csv__ (expected input units are indicated). There is also an Excel .xslx version of this file included which highlights the inputs used in yellow and includes additional comments. The "parameters" sheet from this workbook can be exported as __parameters.csv__ after editing, or the __parameters.csv__ file can be edited directly. Ideally we will use the parameter build function in R to produce the __paramaters.csv__ file.

__Noteworthy plant traits:__

- Whole plant kMax (saturated whole-plant conductance).
- Percent of resistance in leaves (determines how tree conductance is partitioned to woody vs. leaf elements).
- Vulnerability curves (in the form of P12 and P50).
- Basal area/ground area (BA:GA, tree density or BAI).
- Leaf area/basal area (together LA:BA and BA:GA determine LAI).
- Leaf width.
- Root depth in m.
- Maximum carboxylation rate at 25C, Vcmax25 (and associated maximum electron transport rate Jmax25, assumed to be Vmax25 * 1.67).

__Important environment or site traits:__

- Ambient [CO2] Ca (input as ppm).
- Soil hydraulic parameters.
- Elevation.
- Lat/lon.
- Solar noon correction (offset between hour 12 in weather data and actual solar noon at this location)
- Atmospheric "clear sky" transmittance (tau). Calibrates the amount of observed solar radiation considered to be "clear sky" (no clouds), generally between 0.6-0.75. See the equations in the "solarcalc" function if you would like to back-calculate transmittance from a observed clear sky data point.
- The soil parameters we used for many soil types can be found in the __Common Soil Types__ sheet of __parameters - inputs worksheet.xlsx__. __Note:__ Currently only supports using the same soil type for all active layers.
- Soil layers count (Default: 5) can be up to 5. A higher number of soil layers provides a more robust soil water budget simulation, while fewer soil layers may improve performance slightly.

__Model Configuration:__

| Group    	| Model Control       	    | Description          |
| ------------- | ------------------------- | -------------------- |
| Soil     	| __igWaterEnable__   	    | Turns on/off groundwater flow. Values: n (off); y (on). It provides an unlimited source of water at a set potential and distance below the root layers. This water will flow up into the soil layers, and potentially allow layers to fill above field capacity (from the bottom layer up). When disabled (default), the only sources of water input will be the initial fraction of field capacity and observed rainfall (and any water over field capacity will become "drainage").	|
| Soil     	| __i_soilRedEnable__  	    | Turns on/off soil redistribution routine. Values: n (off); y (on). It allows water to flow between soil layers.	|
| Soil     	| __i_soilEvapEnable__ 	    | Turns on/off soil evaporation routine. Values: n (off); y (on). It enables simulation of water evaporation from the surface soil layer.	|
| Climate  	| __i_rainEnable__  	    | Turns on/off rain inputs. Values: n (off); y (on). It allows for precipitation events. Weather data rainfall will be ignored if disabled.	|
| Climate  	| __i_useGSDataStress__     | Turns on/off growing season data for multiple year modeling. Vakyes: n (off); y (on). If enabled, multiple years will be run "sequentially" with on and off seasons defined in seasonlimits_2.0.0.csv and a continuous water budget. When disabled (default), all weather timesteps provided are treated as part of the growing season and the user is expected to truncate individual years to their start/end days. Water budget is reset between years when disabled, treating years as totally independent.	|
| Climate	| __i_useGSDataOpt__ 	    | Turns on/off growing season data for multiple year modeling during BAGA optimization. Values: n (off); y (on).	|
| Hydraulics  	| __i_refilling__ 	    | Turns on/off xylem refilling within a growing season. Values: n (off); y (on). It allows trees to restore lost conductance, however the refilling model is not sufficient to simulate authentic xylem refilling behavior and has __not been thoroughly tested in the current version of the code__.	|
| Hydraulics  	| __i_predawnsMode__ 	    | Turns on/off if model should consider measured pre-dawn water potential values. Values: n (off); y (on). If set to 'y', disables soil simulation and runs from hourly inputs of canopy predawn water potential. These are read from the "rain" column (rain is not used with the soil sim disabled) in MPa. See "dataset - predawns example.csv". This mode is especially useful when comparing to other models which run from canopy predawn measurements, and generally runs significantly faster as it does not need to solve for root layer pressures.	|
| Hydraulics  	| __i_cavitFatigue__ 	    | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism	|
| Hydraulics  	| __i_stemOnly__ 	    | Turns on/off xylem stress hysteresis only in stem xylem. Values: n (off); y (on). When disabled it allows for a weighted estimation of xylem vulnerability to embolism for both stem and roots	|
| BAGA  	| __i_iter_gwEnable__ 	    | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism	|
| BAGA  	| __i_iter_ffcEnable__ 	    | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism	|
| BAGA  	| __i_iter_bagaEnable__     | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism	|
| BAGA  	| __i_iter_useAreaTable__   | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism	|
| BAGA  	| __i_iter_yearsAsCount__   | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism	|
| BAGA  	| __i_iter_runSupplyCurve__ | Turns on/off xylem stress hysteresis to carry effects from previous growing season. Values: n (off); y (on). It allows for a weighted estimation of xylem vulnerability to embolism	|
| Community	| __i_multipleSP__	    | Turns on/off whether our model configuration has 1 species per site (monodominant) or multiple species per site (diverse). Values: n (off); y (on). __NOT READY FOR MULTIPLE SPECIES/PFTs in version 2.0.1__.	|
| Community	| __i_speciesN__	    | Nnumber of species/PFT to run the model.	|
| Forcing files	| __i_ClimateData__	    | Path to file with climate forcing variables __dataset.csv__	|
| Forcing files	| __i_GSData__		    | Path to file with growing season data __seasonlimits.2.0.0.csv__	|

------------

## Weather Data:

See example data for formatting. Weather drivers should be in hourly timesteps and can include multiple years of data. Note that while year values are arbitrary, they should be sequential. For example, if running data for the years 1997 and 2005 these should be numbered sequentially as 1 and 2 (or 1997 and 1998, etc). Inputs:

- Year.
- Julian Day (1-366).
- Hour (0-23).
- Obs. Solar (W m-2).
- Rain (mm).
- Wind (m s-1).
- Tair (C).
- Tsoil (C, if not available substitute air temp).
- D (kPa).

------------

## Outputs:

- Autosave is always enabled regardless of the setting, as this version of the model has no alternative output method. Output files will be generated in the working directory when the run completes.

-Hourly Outputs (see dataheader.csv for full list):
	-Pressures (predawn soil layer pressures, sun and shade "mid-day" canopy pressures, MPa), 
	-Water flows (mmol m-2s-1), 
	-PS assimilation (A, umol s-1m-2 (leaf area)), 
	-Gain-risk optimized stomatal conductance to water (Gw, mmol m-2s-1),
	-Element and whole plant conductances, (k, kghr-1m-2), 
	-Water content and deltas (mm), 
	-Ci

-Summary Outputs (per year, see sumheader.csv for full list):
	-Total Anet (mmol yr-1 m-2(leaf area)),
	-Total E (mm = mm3/mm2(ground area)), 
	-Minimum whole plant conductance during the growing season (kghr-1m-2), 
	-Percent Loss Conductance (PLC, percent, relative to a reference conductance at field capacity),
	-Mean Ci/Ca (+ A weighted Ci/Ca), 
	-Water summary (start/end content, total growing season input (mm)).

------------

### Independent year mode (default)

- Set "Use GS Data" to "n" under "Program Options"
- Use growing season trimmed data (see the example: "dataset.csv").
- Ensure that the growing season limits are defined in "seasonlimits.csv"

The default setting is to reset the tree hydraulics and reset the soil water content to the specified percent of field capacity every year. The years are completely independent, only run in a single dataset for convenience. In this mode, weather data should be trimmed to __only the growing season days__ as in the included dataset.csv (so that the last day of one growing season is followed immediately by the first day of the next). Year values are used for output, and each must be unique and optimally sequential (to determine when new years begin), but the values are otherwise unused by the model and thus do not need to be meaningful values. When running in this mode the growing season limits (seasonlimits.csv) will not be used; All days in the dataset will be considered to be in the growing season.

### Sequential year mode:

In this mode, plant hydraulics will reset between seasons and plant transpiration/productivity will be disabled during the off-season, but soil water budget will continue to be computed. Soil may or may not be refilled to field capacity depending on the availability of off-season precipitation.
- Set "Use GS Data" to "y" under "Program Options".
- Use full-year data (see the example: "dataset - full year example.csv").
- Ensure that the growing season limits are defined in "seasonlimits.csv".
  - __Note:__ The year values in "seasonlimits.csv" are for reference only; The first row of start/end days will be used for the first year of data, etc.
- Soil surface evaporation will also be disabled during the off-season. This is not particularly realistic, but the functionality was intended to answer the question: Is there at minimum enough recorded rainfall to refill the soil? A more robust off-season water simulation would require additional data (snow pack) and simulation of soil behavior under snow and is not provided here.

-------------

## References:

Describing the version 1.0 of the C++ code:
- Venturas MD, JS Sperry, DM Love, EH Frehner, MG Allred, Y Wang, and WRL Anderegg. (2017). A Stomatal Control Model Based on Optimization of Carbon Gain versus Hydraulic Risk Predicts Aspen Sapling Responses to Drought. New Phytologist 220: 836–50.

Describing the gain/risk algorithm used in the model:
- Sperry JS, MD Venturas, WRL Anderegg, M Mencucinni, DS Mackay, Y Wang, and DM Love. (2017). Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40: 816-830

Describing the original hydraulic model the gain-risk optimization was based on:
- Sperry JS, and DM Love (2015) Tansley Review: What plant hydraulics can tell us about plant responses to climate-change droughts. New Phytologist 207: 14-17.
- Sperry JS, Y Wang, BT Wolfe, DS Mackay, WRL Anderegg, NG McDowell, and WT Pockman. (2016). Pragmatic hydraulic theory predicts stomatal responses to climatic water deficits. New Phytologist 212: 577-589
-------------
