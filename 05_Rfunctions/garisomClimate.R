##' @name garisomClimate
##' @title garisomClimate
##' @author German Vargas G.
##' @export
##' 
##' @description
##' This function will trimm the climate data and merge the NOAAA atmospheric CO2 trend data with the growing season data. It creates two files: dataset.csv and seasonlimits_2.0.0.csv
##' 
##' 
##' @param climate_folder The path to the folder containing the climate forcing data.
##' @param data_path The path to the folder where to print the climate forcing data.
##' @param i_site Character vector that contains each site for which we will create the dataset.csv and seasonlimits_2.0.0.csv files. It has to match climate forcing data file names. If NULL, it will create a single pair of files for all sites.
##' @param yeari Initial year
##' @param yearf Final year
##' @param noaa_ca URL to NOAA annual atmospheric Carbon Dioxide concentration records, check here: https://gml.noaa.gov/ccgg/trends/gl_data.html
##'
##'
garisomClimate <- function(climate_folder,data_path,i_site,yeari,yearf,noaa_ca){
  #### download NOAA yearly CO2 inputs ----
  if(RCurl::url.exists(url = noaa_ca)){
    cat("SUCCESS: Fetcing NOAA CO2 data\n")
    co2 <- read.table(noaa_ca)
    colnames(co2) <- c("Year", "ca_ppm", "sd_year")
  } else {
    cat("ERROR: unoable to fetch NOAA CO2 data\n Check url\n Check internet connection")
  }
  co2 <- co2[co2$Year %in% c(yeari:yearf),]
  
  #### subset and export climate forcing drivers ------
  # get a vector of file names
  climate_files <- list.files(climate_folder)
  # get the number of sites
  i_siteN <- length(i_site)
  # create a matrix to store path to files
  climate_paths <- matrix(nrow = i_siteN,ncol=3)
  colnames(climate_paths) <- c("i_site","i_ClimateData","i_GSData")
  # loop over each site to extract the climate forcing data for a given time period
  for(s in 1:i_siteN){
    climate_paths[s,"i_site"] <- i_site[s]# store site s id
    
    # create a directory at the site level
    dir.create(path = paste(data_path,"/",i_site[s],sep = ""),showWarnings = F)
    cat(paste("SUCCESS: Created a directory for",i_site[s],"in:\n",sep = " "))
    cat(paste(data_path,"/",i_site[s],"\n",sep = ""))
    
    # find the files for a given site
    files_per_site <- climate_files[grep(i_site[s],climate_files)]# look for site s
    
    #### hourly climate data --------
    weather_dat_file <- files_per_site[grep("weather",files_per_site)] # select the files that contain hourly weather data
    weather_dat <- read.csv(paste(climate_folder,"/",weather_dat_file,sep = "")) # read hourly weather data
    weather_dat <- weather_dat[weather_dat$Year %in% c(yeari:yearf),] # trim for the year interval
    climate_paths[s,"i_ClimateData"] <- paste(data_path,"/",i_site[s],"/","dataset.csv",sep = "")# save path to the file
    # save file
    try(write.table(x = weather_dat, file = climate_paths[s,"i_ClimateData"], row.names = FALSE,col.names = TRUE,sep = ",",quote = FALSE))
    cat(paste("SUCCESS: Saved hourly climate data for",i_site[s],"in:\n",sep = " "))
    cat(paste(climate_paths[s,"i_ClimateData"],"\n",sep = ""))
    
    rm(weather_dat,weather_dat_file)# delete unnecessary objects
    
    #### season limits data ---------
    GS_dat_file <- files_per_site[grep("GS",files_per_site)] # select the files that contain yearly growing season data
    GS_dat <- read.csv(paste(climate_folder,"/",GS_dat_file,sep=""))# read yearly growing season data
    GS_dat <- GS_dat[GS_dat$Year %in% c(yeari:yearf),]# trim for the year interval
    GS_dat <- merge(x = GS_dat,y = co2[,c("Year","ca_ppm")],by = "Year") # merge with NOAA CO2 data
    climate_paths[s,"i_GSData"] <- paste(data_path,"/",i_site[s],"/","seasonlimits_2.0.0.csv",sep = "")# save path to the file
    # save file
    try(write.table(x = GS_dat,file = climate_paths[s,"i_GSData"],row.names = FALSE,col.names = TRUE,sep = ",",quote = FALSE))
    cat(paste("SUCCESS: Saved growing season data for",i_site[s],"in:\n",sep = " "))
    cat(paste(climate_paths[s,"i_GSData"],"\n",sep = ""))
    rm(GS_dat,GS_dat_file)# delete unnecessary objects
  }
  cat(paste("SUCCESS: DONE\n"))
  return(data.frame(climate_paths))# return a data frame with i_site, i_ClimateData, i_GSData
}