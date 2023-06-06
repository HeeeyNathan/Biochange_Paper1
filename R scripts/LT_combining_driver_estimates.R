### whole time series ####

#combine driver estimates

DriverDir <- "Outputs/Drivers_all_sites"

#make sure it is just the rds trend files
DriverFiles <- list.files(DriverDir)[!grepl("txt",list.files(DriverDir))]
DriverFiles <- list.files(DriverDir)[grepl(".RDS",list.files(DriverDir))]

DriverEsts <- lapply(DriverFiles,function(x){
  
  temp <- readRDS(paste(DriverDir,x,sep="/"))
  
  #add on response from file name
  temp$Response <- strsplit(as.character(x),"__")[[1]][2]
  return(temp)
  
})

DriverEsts <- do.call(rbind,DriverEsts)
names(DriverEsts)[which(names(DriverEsts)=="covariate")] <- "driver"
DriverEsts
saveRDS(DriverEsts,file="outputs/glsTrends_drivers_site_level.rds")
# write.csv(DriverEsts,file="outputs/glsTrends_drivers_site_level_all.rds")

#check we have all data

#get lists of tasks
TaskID <- read.csv("data/LT_ResponseTrends_TaskIDs.csv",as.is=T)
sort(unique(TaskID$Response))
sort(unique(DriverEsts$Response))
#yes!!

##### CLEAN UP --------------------
library(pacman)
# Clear data
rm(list = ls())  # Removes all objects from environment
# Clear packages
p_unload(all)  # Remove all contributed packages
# Clear plots
graphics.off()  # Clears plots, closes all graphics devices
# Clear console
cat("\014")  # Mimics ctrl+L
# Clear mind :)
