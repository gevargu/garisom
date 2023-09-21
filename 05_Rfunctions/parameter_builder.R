##' @name garisompars
##' @title garisompars
##' @author German Vargas G.
##' @export
##' 
##' @description
##' A short description...
##' 
##' 
##' @param parameter_data data.frame with for each parameter of interest in long format with a column for site_codes and a column for species/PFT
##' @param site_codes character vector that contains each site for which we will create the parameter and configuration files
##' @param data_path the path to the folder where to print the parameters, it can be overwritten if included as a column in the main data.frame
##' @param weather_data
##' 
##'
##' name# Parameters file builder function----------------
#'
#' Author: German Vargas G., Sept - 2023
#' 
#' Instructions: please create two folders in your working directory (WD) as: ~WD/scripts/functions
#' 
param_file = function(parameter_data,#data frame with the parameter info, this will change in the PECAN framework
                       site_codes,# list of sites for which we have parameters
                       data_path,# the path to the folder where to print the parameters file
                       region, # where are the plots located
                       weather_data=FALSE,
                       weather_folder,#directory with the weather data
                       model_version="C++",#version of the model
                       co2=410,#atmospheric CO2
                       yeari=1980,# initial year in the weather and GS data, default is 1980
                       yearf=2020,# last year to be simulated in the weather and GS data
                       groundwater_d = 1,#distance to ground water source, default 1 m
                       groundwater = FALSE,#Do plants have access to ground water, default to FALSE
                       soil_redistribution = TRUE,#vertical soil water redistribution
                       soil_evaporation = TRUE,#bare soil evaporation
                       refilling = FALSE,# xylem refilling during the growing season
                       rainfall = TRUE,# rainfall events
                       GS = FALSE,# sequential years mode on
                       predawn_mode = FALSE){
  
  print(x = "Creating files for GARISOM v 2.0.1 Cpp")
  
  #' load miscellaneous functions
  source(paste(getwd(),"/functions/support_functions.R",sep = ""))
  
  #' select the plots to work with
  parameter_data = parameter_data[parameter_data$site %in% site_codes,]
  
  if(model_version=="R"){
    # R version here---------------
    #'
    #' This will write the data using a template csv file for the R model
    #' 
    print(x = "Creating file for the R model")
    parameters_template = read.csv(file = paste(getwd(),"/scripts/functions/R_parameters_template.csv",sep = ""))
    parameters_template[,4] = round(parameters_template[,4],digits = 4)
    
    # Duplicate the default column to the number of sites
    for(p in 1:length(site_codes)){
      column_name = paste("value",p,sep = "")
      parameters_template$x1 = parameters_template[,4]
      colnames(parameters_template)[p+3] = column_name
    }
    
    ## Site parameters
    parameters_template[60,c(4:(length(site_codes)+3))] = site_codes#plot name
    
    site_index = map(site_codes, function(x) which(sites$plot == x)) %>% unlist
    parameters_template[1,c(4:(length(site_codes)+3))] = sites$lat[site_index]#latitude
    parameters_template[2,c(4:(length(site_codes)+3))] = sites$long[site_index]#longitude
    parameters_template[3,c(4:(length(site_codes)+3))] = sites$elev_m[site_index]#altitude
    parameters_template[4,c(4:(length(site_codes)+3))] = sites$slope_perc[site_index]#slope
    parameters_template[5,c(4:(length(site_codes)+3))] = asp_fun(asp_card = sites$aspect[site_index])#slope aspect
    parameters_template[6,c(4:(length(site_codes)+3))] = sites$total_ba_ga_m2ha[site_index]# total basal area
    
    if(is.na(LAI_data)){
      print(x = "No LAI data, inputing 2 for each plot")
      parameters_template[7,c(4:(length(site_codes)+3))] = rep(x = 2,lenght.out=length(site_codes)) ## No data
    }else{
      #We need to collect LAI data this summer during peak canopy cover, or estimate it with satellite images...
    }
    
    parameters_template[8,c(4:(length(site_codes)+3))] = rep(x = co2,lenght.out=length(site_codes))#atmospheric co2
    #parameters_template[9,c(4:(length(site_codes)+3))] = 0.65#, atmospheric transmittance, using default values for now
    
    SN = sapply(X = as.numeric(parameters_template[2,c(4:(length(site_codes)+3))]),FUN = function(x) {
      SN= Solar_noon_correction(longitude = x,UTC.weather = 8)
      return(SN)})
    parameters_template[10,c(4:(length(site_codes)+3))] = SN[1,]# solar noon correction
    
    #parameters_template[11,c(4:(length(site_codes)+3))] = rep(x = 0.97, lenght.out=length(site_codes))#emissivity, using default values for now
    #parameters_template[12,c(4:(length(site_codes)+3))] = rep(x = 1, lenght.out=length(site_codes))#understory heigh above soil, using default values
    
    ## Soil parameters
    soil_index = map(site_codes, function(x) which(soils$Site == x)) %>% unlist
    parameters_template[13,c(4:(length(site_codes)+3))] = soils$Texture_mean[soil_index]# soil texture
    parameters_template[14,c(4:(length(site_codes)+3))] = rep(x = 5, lenght.out=length(site_codes))# number of soil layers, using default of 5 now
    parameters_template[15,c(4:(length(site_codes)+3))] = soils$Rock_fraction[soil_index]#rock fraction
    parameters_template[16,c(4:(length(site_codes)+3))] = rep(x = 50, lenght.out=length(site_codes))#whole plant water resistance in rhizosfere
    parameters_template[17,c(4:(length(site_codes)+3))] = rep(x = 100, lenght.out=length(site_codes))#percent field capacity for starting the season
    parameters_template[18,c(4:(length(site_codes)+3))] = rep(x = 0.5, lenght.out=length(site_codes))#fraction that field capacity is of saturation
    parameters_template[19,c(4:(length(site_codes)+3))] = rep(x = 0.94, lenght.out=length(site_codes))#absorptivity of soil surface for solar
    
    ## Tree parameters
    
    SPCODE = substring(parameters_template[60,c(4:(length(site_codes)+3))],first = 1,4)
    
    for (sp in 1:length(site_codes)){# loop over site_codes, over species
      parameters_template[20,3+sp] = 1#max radius of root system per max depth, for now using default
      parameters_template[21,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"Rooting_depth"]#maximum rooting depth
      parameters_template[22,3+sp] = 1#leaf angle parameter; CN 15.4
      parameters_template[23,3+sp] = as.numeric(parameters_template[7,3+sp])/as.numeric(parameters_template[6,3+sp])*1e4#leaf area per basal area
      parameters_template[24,3+sp] = 25#saturated leaf resistance
      parameters_template[25,3+sp] = Calculate_Kmax(LABA = as.numeric(parameters_template[23,3+sp]),LSC = 10,Rleaf = as.numeric(parameters_template[24,3+sp]))#kmax of tree in kg hr-1 m-2 MPa-1 per basal area
      parameters_template[26,3+sp] = height_data[height_data$Species == SPCODE[sp] & height_data$plots == site_codes[sp],"h_m"]# height per species per plot
      parameters_template[27,3+sp] = 0.04096#leaf width, using default for now...
      parameters_template[28,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"b"]*0.5#weibull b for each root element
      parameters_template[29,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"c"]*0.5#weibull c for each root element
      parameters_template[30,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"b"]#weibull b for stem
      parameters_template[31,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"c"]#weibull c for stem
      parameters_template[32,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"b"]*1.5#weibull b for each root element
      parameters_template[33,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"c"]*1.5#weibull c for each root element
      parameters_template[34,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"Vcmax"]#umol m-2 s-1; maximum carboxylation rate (vmax) at 25C
      parameters_template[35,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"Vcmax"]*1.67#umol m-2 s-1; maximum electron transport rate (jmax) at 25C
      parameters_template[66,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"sapwood_t"]#sapwood turnover rate
      parameters_template[67,3+sp] = traits[traits$SpeciesID == SPCODE[sp],"conduit_d"]#Xylem conduit diameter
      
      #weather file name per plot
      parameters_template[62,3+sp] = paste("weather_plot_",site_codes[sp],".csv",sep = "")
      parameters_template[63,3+sp] = paste("GS_plot_",site_codes[sp],".csv",sep = "")
    }
    
    # Model control
    parameters_template[51,c(4:(length(site_codes)+3))] = model_controls[parameters_template[51,"parameter"]]#perform xylem refilling during GS? yes = 1, no = 0
    parameters_template[52,c(4:(length(site_codes)+3))] = model_controls[parameters_template[52,"parameter"]]#turns on/off groundwater flow? yes = 1, no = 0
    parameters_template[53,c(4:(length(site_codes)+3))] = model_controls[parameters_template[53,"parameter"]]#turns on/off soil redistribution routine? yes = 1, no = 0
    parameters_template[54,c(4:(length(site_codes)+3))] = model_controls[parameters_template[54,"parameter"]]#turns on/off soil evaporation routine? yes = 1, no = 0
    parameters_template[55,c(4:(length(site_codes)+3))] = model_controls[parameters_template[55,"parameter"]]#turns on/off rain inputs? yes = 1, no = 0
    parameters_template[56,c(4:(length(site_codes)+3))] = model_controls[parameters_template[56,"parameter"]]#turns on/off drainage blockage allowing saturation? yes = 1, no = 0
    parameters_template[57,c(4:(length(site_codes)+3))] = model_controls[parameters_template[57,"parameter"]]#turns on/off resetting kmax at the beginning of GS? yes = 1, no = 0
    parameters_template[58,c(4:(length(site_codes)+3))] = model_controls[parameters_template[58,"parameter"]]#turns on/off computing water budget during nonGS? yes = 1, no = 0
    parameters_template[59,c(4:(length(site_codes)+3))] = model_controls[parameters_template[59,"parameter"]]#turns on/off soil evaporation during non growing season? yes = 1, no = 0
    parameters_template[61,c(4:(length(site_codes)+3))] = weather_folder
    parameters_template[64,c(4:(length(site_codes)+3))] = model_controls[parameters_template[64,"parameter"]]#perform cavitation fatigue routine? yes = 1, no = 0
    parameters_template[65,c(4:(length(site_codes)+3))] = model_controls[parameters_template[65,"parameter"]]#update pcrit with CF? yes = 1, no = 0
    
    return(parameters_template)
    }else{
      # C++ version-----
      #'
      #' This will write the data using a template csv file for the C++ model
      #' 
      #' Some important notes:
      #' This function runs in parallel, so add as many sites as you want. It will produce 5 files per plot/site:
      #' parameters.csv (plant, site params and program options)
      #' nametable.csv (maps parameter names to row/col locations in the parameters.csv sheet)
      #' dataset.csv (hourly weather drivers)
      #' dataheader.csv (a header row for the hourly data output)
      #' sumheader.csv (a header row for the summary data output)
      #' seasonlimits.csv (growing season limits, only required if using "sequential year mode"
      #' 
      
      
      # Load packages and functions to run function in parallel
      require(parallel)
      require(doParallel)
      registerDoParallel(cores = detectCores())
      require(foreach)
      
      # Create list files for storing dataframes and temporal objects
      parameters_list = list()
      nametable_list = list()
      dataset_list = list()
      dataheader_list = list()
      sumheader_list = list()
      seasonlimits_list = list()
      
      SN = list()# solar noon correction list
      SPAR = list()# soil parameters list
      
      # Run the loop in parallel for creating each parameter file
      n_sites = length(site_codes)
      foreach(p = 1:n_sites) %do% {#parallelization over folders with files
        print(x = paste("Extracting information for ",site_codes[p],sep = ""))
        
        #' 
        #' Create a directory at the plot level
        #' 
        
        dir.create(path = paste(data_path,"/",site_codes[p],sep = ""),showWarnings = F)
        
        #'
        #' Read parameter files and save them as objects of a list
        #'
        parameters_list[[p]] = read.csv(file = paste(getwd(),"/functions/cpp_parameters_template.csv",sep = ""),
                                         header = FALSE)
        
        #'
        #'
        #'
        # General information----------
        parameters_list[[p]][17,2] = parameter_data[parameter_data[,"site"] == site_codes[p],"species"]#write species name or code (e.g., potr)
        parameters_list[[p]][18,2] = region#write region name (e.g., Mountain West)
        parameters_list[[p]][19,2] = site_codes[p]#write site or plot name (e.g., POTR-M-E)
        
        # Stand parameters----------
        parameters_list[[p]][42,2] = parameter_data[parameter_data$site == site_codes[p],"lat"]#latitude
        parameters_list[[p]][43,2] = parameter_data[parameter_data$site == site_codes[p],"lon"]#longitude
        parameters_list[[p]][44,2] = ifelse(test = is.na(parameter_data[parameter_data$site %in% site_codes[p],"slope_topo"]),#slope inclination
                                            yes = median(parameter_data[,"slope_topo"],na.rm = TRUE),
                                            no = parameter_data[parameter_data$site %in% site_codes[p],"slope_topo"])
        parameters_list[[p]][45,2] = asp_fun(asp_card = parameter_data[parameter_data$site == site_codes[p],"aspect"])#slope aspect
        parameters_list[[p]][46,2] = 1#leaf angle parameter; CN 15.4
        parameters_list[[p]][47,2] = parameter_data[parameter_data$site == site_codes[p],"LAI"]
        parameters_list[[p]][48,2] = 1 #vegetation height above soil
        parameters_list[[p]][49,2] = parameter_data[parameter_data$site == site_codes[p],"elevation"] # site elevation
        SN[[p]] = Solar_noon_correction(longitude = as.numeric(parameters_list[[p]][43,2]),UTC.weather = 8)#solar noon correction function     
        parameters_list[[p]][50,2] = SN[[p]][1]#solar noon correction
        parameters_list[[p]][51,2] = parameter_data[parameter_data$site == site_codes[p],"basal_area"]# extract basal area per ground area of a given plot
        
        # Soil parameters--------
        parameters_list[[p]][54,2] = 5# number of soil layers, using default of 5 for now
        parameters_list[[p]][55,2] = 0.5#fraction of field capacity
        parameters_list[[p]][56,2] = 100#initial field capacity
        parameters_list[[p]][57,2] = parameter_data[parameter_data$site %in% site_codes[p],"rock_fraction"]# fraction of soils that is rocks 
        parameters_list[[p]][58,2] = 0.94#soil absorptivity of solar radiation, using default
        SPAR[[p]] = soil_parameters_from_texture(texture = parameter_data[parameter_data$site %in% site_codes[p],"t_usda_tex"])#texture parameters function
        parameters_list[[p]][61:65,3:6] = rep(x = SPAR[[p]],each=5)# soil texture parameters function
        parameters_list[[p]][69,2] = 0.00075#pressure increment for curve generation, using default
        
        # Plant parameters------
        parameters_list[[p]][26,3] = parameter_data[parameter_data$site %in% site_codes[p],"weibull_b"]#stem element b
        parameters_list[[p]][27,3] = parameter_data[parameter_data$site %in% site_codes[p],"weibull_c"]#stem element c
        parameters_list[[p]][26,2] = as.numeric(parameters_list[[p]][26,3])*0.5#root element b
        parameters_list[[p]][27,2] = as.numeric(parameters_list[[p]][27,3])*0.5#root element c
        parameters_list[[p]][26,4] = as.numeric(parameters_list[[p]][26,3])*1.5#leaf element b
        parameters_list[[p]][27,4] = as.numeric(parameters_list[[p]][27,3])*1.5#leaf element c
        parameters_list[[p]][28,4] = parameter_data[parameter_data$site %in% site_codes[p],"leaf_resist"]#leaf resistance at Ksat, for now 25%
        parameters_list[[p]][33,2] = parameter_data[parameter_data$site %in% site_codes[p],"leaf_cond"]#leaf conductance
        parameters_list[[p]][34,2] = parameter_data[parameter_data$site %in% site_codes[p],"canopy_height"]#Average tree height
        parameters_list[[p]][35,2] = parameter_data[parameter_data$site %in% site_codes[p],"LA_BA"]#leaf area per basal area
        parameters_list[[p]][32,2] = Calculate_Kmax(LABA = as.numeric(parameters_list[[p]][35,2]),#Kmax of plant
                                                     LSC = as.numeric(parameters_list[[p]][33,2]),
                                                     Rleaf = as.numeric(parameters_list[[p]][28,4]))#kmax of tree in kg hr-1 m-2 MPa-1 per basal area
        parameters_list[[p]][36,2] = parameter_data[parameter_data$site %in% site_codes[p],"leaf_width"]
        parameters_list[[p]][37,2] = 50#whole plant water resistance in rhizosfere, using default
        parameters_list[[p]][38,2] = 1#Root aspect ratio, using default
        parameters_list[[p]][39,2] = root.beta(rootdepth = parameter_data[parameter_data$site %in% site_codes[p],"rooting_depth"])#root beta
        parameters_list[[p]][78,2] = parameter_data[parameter_data$site %in% site_codes[p],"Vcmax25"]#umol m-2 s-1; maximum carboxylation rate (vmax) at 25C
        parameters_list[[p]][79,2] = parameter_data[parameter_data$site %in% site_codes[p],"Jmax25"]#umol m-2 s-1; maximum electron transport rate (jmax) at 25C
        
        # Atmospheric parameters------
        parameters_list[[p]][76,2] = co2# atmospheric CO2, this got me to think!!!!!
        parameters_list[[p]][95,2] = 0.97# sky emissivity, using default values
        parameters_list[[p]][96,2] = 0.75# atm transmittance, using default values
        parameters_list[[p]][98,2] = 0.21# O2 Mole fraction, using default values
        
        # Model controls--------
        parameters_list[[p]][24,9] = 0# water potential of groundwater
        parameters_list[[p]][25,9] = groundwater_d# distance to ground water
        parameters_list[[p]][26,9] = ifelse(test = groundwater == TRUE,yes = "y",no = "n")# access to ground water
        parameters_list[[p]][27,9] = ifelse(test = soil_redistribution == TRUE,yes = "y",no = "n")# vertical soil water redistribution
        parameters_list[[p]][28,9] = ifelse(test = soil_evaporation == TRUE, yes = "y",no = "n")# soil water evaporation
        parameters_list[[p]][29,9] = ifelse(test = refilling == TRUE,yes = "y",no = "n")# xylem refilling during the growing season
        parameters_list[[p]][30,9] = ifelse(test = rainfall == TRUE,yes = "y",no = "n")# rainfall events
        parameters_list[[p]][34,9] = ifelse(test = GS == TRUE,yes = "y",no = "n")# is the water budget sequential according to the growing season
        #parameters_list[[p]][35,9] = ifelse(test = autosave == TRUE, yes = "y", no = "n")# Autosave is always enabled
        parameters_list[[p]][36,9] = ifelse(test = predawn_mode == TRUE,yes = "y",no = "n")# To run from hourly predawn canopy water potentials (psoil + pgrav, MPa) -- predawn input read from "rain" column
        
        # Exporting files----------
        #'
        #' Export file with the parameters for the species in the given plot directory by:
        #' 
        #' a) Creates a species directory within the plot directory
        #' 
        #dir.create(path = paste(data_path,"/",site_codes[p],"/",sp_code[s],sep = ""),showWarnings = F)
        #' b) exporting the parameters.csv file
        try(write.table(x = parameters_list[[p]],
                        file = paste(data_path,"/",site_codes[p],"/","parameters.csv",sep = ""),
                        row.names = FALSE,col.names = FALSE,sep = ",",quote = FALSE))
        #'
        #' c) exporting the nametable.csv
        #'
        nametable_list[[p]] = read.csv(file = paste(getwd(),"/functions/nametable.csv",sep = ""),
                                        header = FALSE)
        try(write.table(x = nametable_list[[p]],
                        file = paste(data_path,"/",site_codes[p],"/","nametable.csv",sep = ""),
                        row.names = FALSE,col.names = FALSE,sep = ",",quote = FALSE))
        
        if(weather_data == TRUE){
          #'
          #' d) exporting the dataset.csv (hourly weather drivers)
          #' 
          dataset_list[[p]] = read.csv(paste(weather_folder,"/","weather_plot_",site_codes[p],".csv",sep = ""))
          dataset_list[[p]] = dataset_list[[p]][dataset_list[[p]]$Year %in% c(yeari:yearf),]
          try(write.table(x = dataset_list[[p]],
                          file = paste(data_path,"/",site_codes[p],"/","dataset.csv",sep = ""),
                          row.names = FALSE,col.names = TRUE,sep = ",",quote = FALSE))
        }else{
          # do nothing for now
        }
        
        #' 
        #' e) exporting the dataheader.csv (a header row for the hourly data output)
        #'
        dataheader_list[[p]] = read.csv(file = paste(getwd(),"/functions/dataheader.csv",sep = ""),
                                         header = FALSE)
        try(write.table(x = dataheader_list[[p]],
                        file = paste(data_path,"/",site_codes[p],"/","dataheader.csv",sep = ""),
                        row.names = FALSE,col.names = FALSE,sep = ",",quote = FALSE))
        #' 
        #' f) exporting the sumheader.csv (a header row for the summary data output)
        #'
        sumheader_list[[p]] = read.csv(file = paste(getwd(),"/functions/sumheader.csv",sep = ""),
                                        header = FALSE)
        try(write.table(x = sumheader_list[[p]],
                        file = paste(data_path,"/",site_codes[p],"/","sumheader.csv",sep = ""),
                        row.names = FALSE,col.names = FALSE,sep = ",",quote = FALSE))
        
        if(weather_data == TRUE){
          #'
          #' g) Exporting the growing season data
          #'
          seasonlimits_list[[p]] = read.csv(paste(weather_folder,"/","GS_plot_",site_codes[p],".csv",sep = ""))
          seasonlimits_list[[p]] = seasonlimits_list[[p]][seasonlimits_list[[p]]$Year %in% c(yeari:yearf),]
          try(write.table(x = seasonlimits_list[[p]],
                          file = paste(data_path,"/",site_codes[p],"/","seasonlimits.csv",sep = ""),
                          row.names = FALSE,col.names = TRUE,sep = ",",quote = FALSE))
        } else{
          #Do nothing for now!
        }

        #'
        #'
        #'
      }
      # END of C++ version-------
  }
}
#' 
#' 
# END--------------
#' 
#' 