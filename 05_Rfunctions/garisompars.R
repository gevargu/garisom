##' @name garisompars
##' @title garisompars
##' @author German Vargas G.
##' @export
##' 
##' @description
##' Finds every parameter value for a given site within a region and creates the parameters_2.0.0.csv and the configuration_2.0.0.csv files to be used in GARISOM > v 2.0.0.
##' This function runs in parallel, so add as many sites as you want.
##' 
##' 
##' @param parameter_data Data.frame with for each parameter of interest in long format with a column for site_codes and a column for species/PFT.
##' @param i_site Character vector that contains each site for which we will create the parameter and configuration files.
##' @param i_region The broader geographical area containing the sites we want to model.
##' @param data_path The path to the folder where to print the parameters.
##' @param co2 Atmospheric CO2 at the beginning of the simulation.
##' @param i_gWaterEnable Turns on/off groundwater flow. Default to FALSE.
##' @param i_soilRedEnable Turns on/off soil redistribution routine. Default to TRUE.
##' @param i_soilEvapEnable Turns on/off soil evaporation routine. Default to TRUE.
##' @param i_rainEnable Turns on/off precipitation inputs. Default to TRUE.
##' @param i_useGSDataStress Turns on/off growing season data for multiple year modeling. Default to TRUE (sequential year mode).
##' @param i_useGSDataOpt Turns on/off growing season data for multiple year modeling during BAGA optimization. Default to TRUE.
##' @param i_refilling Turns on/off xylem refilling within a growing season. Default to FALSE.
##' @param i_predawnsMode Turns on/off if model should consider measured pre-dawn water potential values. Default to FALSE.
##' @param i_cavitFatigue Turns on/off xylem stress hysteresis to carry effects from previous growing season. Default to FALSE.
##' @param i_stemOnly Turns on/off xylem stress hysteresis only in stem xylem. Default to TRUE.
##' @param i_iter_gwEnable Not tested in version 2.0.1. Default to FALSE.
##' @param i_iter_ffcEnable Not tested in version 2.0.1. Default to FALSE.
##' @param i_iter_bagaEnable Not tested in version 2.0.1. Default to FALSE.
##' @param i_iter_useAreaTable Not tested in version 2.0.1. Default to FALSE.
##' @param i_iter_yearsAsCount Not tested in version 2.0.1. Default to FALSE.
##' @param i_iter_runSupplyCurve Not tested in version 2.0.1. Default to FALSE.
##' @param i_multipleSP Turns on/off whether our model configuration has 1 species per site (monodominant) or multiple species per site (diverse). Not ready for multiple species in version 2.0.1. Default to FALSE.
##' @param i_speciesN Number of species/PFT to run the model.
##' @param i_ClimateData Path to file with climate forcing variables dataset.csv
##' @param i_GSData Path to file with growing season data seasonlimits.2.0.0.csv
##'
garisompars <- function(parameter_data,#data frame with the parameter info, this will change in the PECAN framework
                      i_site,# list of sites for which we have parameters
                      data_path,# the path to the folder where to print the parameters file
                      i_gWaterEnable = FALSE,
                      i_soilRedEnable = TRUE,
                      i_soilEvapEnable = TRUE,
                      i_rainEnable = TRUE,
                      i_useGSDataStress = TRUE,
                      i_useGSDataOpt = TRUE,
                      i_refilling = FALSE,
                      i_predawnsMode = FALSE,
                      i_cavitFatigue = FALSE,
                      i_stemOnly = TRUE,
                      i_iter_gwEnable = FALSE,
                      i_iter_ffcEnable = FALSE,
                      i_iter_bagaEnable	= FALSE,
                      i_iter_useAreaTable	= FALSE,
                      i_iter_yearsAsCount	= FALSE,
                      i_iter_runSupplyCurve	= FALSE,
                      i_multipleSP	= FALSE,
                      i_speciesN	= 1){
  cat(x = "Loading default values\n")# for now based on ASPEN values, will based on Venturas et al. 2021
  default.parameters <- matrix(nrow = 68,ncol = 1)
  colnames(default.parameters) <- "PFT_1" 
  row.names(default.parameters) <- c("i_sp","i_region","i_site","i_latitude","i_longitude","i_elevation","i_slopeI","i_slopeA","i_gWaterP","i_gWaterDist","i_atmTrans","i_solarNoon","i_emiss","i_co2AmbPPM","i_layers","i_fieldCapFrac","i_fieldCapPercInit","i_rockFrac","i_soilAbsSol","i_rhizoPer","i_texture","i_baperga","i_leafAreaIndex","i_soilXHeight","i_height","i_treeToPhotoLAI","i_leafPerBasal","i_leafWidth","i_leafAngleParam","i_aspect","i_rootDepth","i_leafPercRes","i_kmaxTree","i_pinc","i_rootP12","i_rootP50","i_stemP12","i_stemP50","i_leafP12","i_leafP50","i_sapwoodT","i_conduitDiam","i_qMax","i_vmax25","i_jmax25","i_kc25","i_ko25","i_comp25","i_thetaC","i_havmax","i_hdvmax","i_svvmax","i_lightCurv","i_lightComp","i_hajmax","i_hdjmax","i_svjmax","i_iter_gwInc","i_iter_gwStart","i_iter_gwEnd","i_iter_ffcInc","i_iter_ffcStart","i_iter_ffcEnd","i_iter_bagaInc","i_iter_bagaStart","i_iter_bagaEnd","i_iter_bagaRef","i_iter_bagaCutoff")
  default.parameters[,"PFT_1"] <- c("Populus tremuloides","Rocky Mountains","POTR-PFT","37.47543","-108.16058","3071","9.6","337.5","0","1","0.75","-15.21070533","0.97","390","5","0.5","100","0.24555","0.94","50","loam","88.55597","1.990195159","1","22.60892422","3.82","1736.379839","0.044555168","1","1","1.15","25","102.5965238","0.00075","-0.501002","-0.803847","-1.002004","-1.607694","-1.503006","-2.411541","0.1316129","23.55","0.3","52.02396278","86.88001785","0.0004048","0.27833131","0.00004275","0.98","73637","149252","486","0.9","30","50300","152044","495","-1.2","12","4","0.2","0.1","1","5","5","150","2","0.5")
  
  cat(x = "Creating files for GARISOM v > 2.0.0 Cpp")
  
  # select the sites for which we will create the files
  parameter_data <- parameter_data[parameter_data$i_site %in% i_site,]
  # number of sites
  n_sites <- length(i_site)
  
  for(s in 1:n_sites){
    cat(x = paste("\nExtracting information for: ",i_site[s],sep = ""))
    
    # Create a directory at the site level
    dir.create(path = paste(data_path,"/",i_site[s],sep = ""),showWarnings = F)
    
    ## Configuration file: configuration_2.0.0.csv ----------------
    config.file <- matrix(nrow = 3,ncol = 23)
    config.file[1,] <- c("soil","soil","soil","climate","climate","climate","hydraulics","hydraulics","hydraulics","hydraulics","baga","baga","baga","baga","baga","baga","community","community","forcing_files","forcing_files","forcing_files","forcing_files","NULL")
    config.file[2,] <- c("i_gWaterEnable","i_soilRedEnable","i_soilEvapEnable","i_rainEnable","i_useGSDataStress","i_useGSDataOpt","i_refilling","i_predawnsMode","i_cavitFatigue","i_stemOnly","i_iter_gwEnable","i_iter_ffcEnable","i_iter_bagaEnable","i_iter_useAreaTable","i_iter_yearsAsCount","i_iter_runSupplyCurve","i_multipleSP","i_speciesN","i_ClimateData","i_GSData","i_dataheader","i_sumheader","NULL")
    config.file[c(1,2,3),23] <- "NULL"# last column should say NULL, this will avoid problems in Cpp
    model_controls <- config.file[2,]
    N <- length(model_controls)-1
    #### Inputting the model configurations ------------------------------------
    for (c in 1:N) {
      pos.config <- as.vector(which(x = config.file==model_controls[c],arr.ind = TRUE))
      pos.config[1] <- pos.config[1]+1# get the right row
      if(model_controls[c] %in% c("i_ClimateData","i_GSData","i_dataheader","i_sumheader")){
        config.file[pos.config[1],pos.config[2]] <- parameter_data[parameter_data$i_site == i_site[s],model_controls[c]]
      } else {
        if(isTRUE(get(model_controls[c]))){
          config.file[pos.config[1],pos.config[2]] <- "y"# assign a y value
        } else {
          config.file[pos.config[1],pos.config[2]] <- "n"# assign a n value
        }
      }
    }
    rm(N,pos.config)
    # exporting configuration file ---------------------------------------------
    cat(x = paste("\nSUCCESS! Saving model configuration for: ",i_site[s]))
    try(write.table(x = config.file,
                    file = paste(data_path,"/",i_site[s],"/","configuration_2.0.0.csv",sep = ""),
                    row.names = FALSE,col.names = FALSE,sep = ",",quote = FALSE))
    
    ## Parameters file: paremeters_2.0.0.csv ----------------
    #### 1. create parameter matrix where each line is a species/pft -----------
    param.file <- matrix(nrow = 2+i_speciesN,ncol = 69)
    param.file[1,] <- c("site","site","site","site","site","site","site","site","site","site","atmosphere","atmosphere","atmosphere","atmosphere","soil","soil","soil","soil","soil","soil","soil","stand","stand","stand","stand","stand","stand","tree","tree","tree","tree","hydraulics","hydraulics","hydraulics","hydraulics","hydraulics","hydraulics","hydraulics","hydraulics","hydraulics","hydraulics","hydraulics","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","photosynthesis","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","BAGA_optimizer","NULL")
    param.file[2,] <- c("i_sp","i_region","i_site","i_latitude","i_longitude","i_elevation","i_slopeI","i_slopeA","i_gWaterP","i_gWaterDist","i_atmTrans","i_solarNoon","i_emiss","i_co2AmbPPM","i_layers","i_fieldCapFrac","i_fieldCapPercInit","i_rockFrac","i_soilAbsSol","i_rhizoPer","i_texture","i_baperga","i_leafAreaIndex","i_soilXHeight","i_height","i_treeToPhotoLAI","i_leafPerBasal","i_leafWidth","i_leafAngleParam","i_aspect","i_rootDepth","i_leafPercRes","i_kmaxTree","i_pinc","i_rootP12","i_rootP50","i_stemP12","i_stemP50","i_leafP12","i_leafP50","i_sapwoodT","i_conduitDiam","i_qMax","i_vmax25","i_jmax25","i_kc25","i_ko25","i_comp25","i_thetaC","i_havmax","i_hdvmax","i_svvmax","i_lightCurv","i_lightComp","i_hajmax","i_hdjmax","i_svjmax","i_iter_gwInc","i_iter_gwStart","i_iter_gwEnd","i_iter_ffcInc","i_iter_ffcStart","i_iter_ffcEnd","i_iter_bagaInc","i_iter_bagaStart","i_iter_bagaEnd","i_iter_bagaRef","i_iter_bagaCutoff","NULL")
    param.file[3,69] <- "NULL"# last column should say NULL, this will avoid problems in Cpp
    params <- param.file[2,]
    Npar <- length(params)-1
    #### 2. loop over parameters and save them ---------------------------------
    # NOTE: if parameter value is empty or NA or does not exist in the df we input a default
    # TODO: add default values for generic PFTs according to the literature
    
    for(i in 1:i_speciesN){
      nsp <- i# account for the number of species
      for(p in 1:Npar){
        pos.params <- as.vector(which(x = param.file==params[p],arr.ind = TRUE))
        pos.params[1] <- pos.params[1] + nsp # the parameter row plus the species
        if(params[p] %in% colnames(parameter_data)){# user values
          param.file[pos.params[1],pos.params[2]] <- parameter_data[parameter_data$i_site %in% i_site[1],params[p]]
        } else { # default values
          param.file[pos.params[1],pos.params[2]] <- default.parameters[params[p],"PFT_1"]
        }
      }
    }
    rm(Npar,pos.params)
    
    #### 3. exporting parameters file ------------------------------------------
    cat(x = paste("\nSUCCES! Saving model parameters for: ",i_site[s]))
    
    try(write.table(x = param.file,
                    file = paste(data_path,"/",i_site[s],"/","parameters_2.0.0.csv",sep = ""),
                    row.names = FALSE,col.names = FALSE,sep = ",",quote = FALSE))
  }
  cat(x = paste("\nDONE"))
}