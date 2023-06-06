### whole time series ####

#combine trends across all countries afterwards

trendsDir <- "Outputs/Site_trends"

#make sure it is just the rds trend files
trendsFiles <- list.files(trendsDir)[!grepl("txt",list.files(trendsDir))]
trendsFiles <- list.files(trendsDir)[grepl(".RDS",list.files(trendsDir))]

siteTrends <- lapply(trendsFiles,function(x){
  
  temp <- readRDS(paste(trendsDir,x,sep="/"))
  
  #add on response from file name
  temp$Response <- strsplit(as.character(x),"__")[[1]][2]
  return(temp)
  
})

siteTrends <- do.call(rbind,siteTrends)
names(siteTrends)[which(names(siteTrends)=="siteID")] <- "site_id"
siteTrends
saveRDS(siteTrends,file="outputs/glsTrends_site_level.rds")
# write.csv(siteTrends,file="outputs/glsTrends_site_level_all.csv")

#check we have all data

#get lists of tasks
TaskID <- read.csv("data/LT_ResponseTrends_TaskIDs.csv",as.is=T)
sort(unique(TaskID$Response))
sort(unique(siteTrends$Response))
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
