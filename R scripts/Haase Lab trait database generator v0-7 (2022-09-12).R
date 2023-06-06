###---v0-7 produced on 2022.09.12
###---Code authored by James S. Sinclair for use by the Haase Lab, Division of River Ecology and Conservation, Senckenberg Research Institute

##Load required libraries (install if necessary)
install.packages("Hmisc")
install.packages("stringr")
library(Hmisc)
library(stringr)

##Load your taxa list with verified freshwaterecology.info (FWE) taxa names in the first column and AQEM codes in the second
#make sure the list was verified using v8 of FWE
taxa<-read.csv(file.choose(), header=TRUE)
#taxa<-read.delim(file.choose(), header=TRUE, sep=";"), use this if you have a semicolon separated file

str(taxa) #make sure everything reads in correctly
Encoding(taxa[,1])<-"latin1" #fix encoding and convert to uppercase
taxa[,1]<-toupper(as.character(taxa[,1]))

##Load the FWE taxa dictionary titled "FWE taxa dictionary v8 - w traits and DISPERSE.csv"
dict<-read.csv(file.choose(), header=TRUE) #trait categories must have unique names and modalities should follow those in FWE
#note all Tubificidae have been re-classified to Naididae
#any taxa in the DISPERSE dataset have also had those values used for the Tachet traits for body size, 
#life duration, reproductive cycles, and dispersal
#also manually added Nematoda Gen. sp. (8813), Coenagrioninae Gen. sp. (AQEM 19850), Mermithidae Gen. sp. (AQEM 9249), and Crangonyx/Niphargus sp. (AQEM 20204)
#manually averaged traits for Oligochaeta Gen. sp. and Crangonyx/Niphargus sp.

##Convert names to uppercase
Encoding(dict$name)<-"latin1"
Encoding(dict$fam)<-"latin1"
Encoding(dict$name1)<-"latin1"
Encoding(dict$name2)<-"latin1"
dict$name<-toupper(as.character(dict$name))
dict$fam<-toupper(as.character(dict$fam))
dict$name1<-toupper(as.character(dict$name1))
dict$name2<-toupper(as.character(dict$name2))

##--Get the unique trait names from the dictionary
traits<-unique(gsub("_.*", "", x=colnames(dict)[6:length(colnames(dict))]))

##--Create a new dataframe to hold all the trait values for each taxon from your taxa list
taxa2<-taxa[,1:2]
colnames(taxa2)[1:2]<-c("name", "aqem")
taxa2$taxon<-NA
taxa2$aqem2<-NA
taxa2$gen<-NA
taxa2$sp<-NA

##--Put in the genus and species names, NAs will be placed in these columns for AQEMs above genus level
#this step let you know which entry from the dictionary you will be pulling values from
for (i in 1:length(taxa2$aqem)) {
  sub<-dict[which(dict$aqem==taxa2$aqem[i]),]
  name1<-sub$name1
  name2<-sub$name2
  if (length(sub$aqem)==0) { #if it cant find an AQEM code, usually only happens with Ad., Lv., Agg., or Gr. entries
    name0<-strsplit(taxa2$name[i], split=" ")[[1]] #split the name at the blanks
    sub<-dict[grep(paste(name0[1], name0[2], sep=" "), dict$name)[1],] #then search using the first two parts and take the first entry
    name1<-sub$name1
    name2<-sub$name2
  }
  if (name2 %in% "GEN. SP.") {
    name1<-NA
    name2<-NA
  }
  taxa2$taxon[i]<-sub$name
  taxa2$aqem2[i]<-sub$aqem #this can be different if it had to search using the taxon name
  taxa2$gen[i]<-name1
  taxa2$sp[i]<-name2
  if (taxa2$name[i]=="SPHAERIUM CORNEUM/OVALE") { #a specific case in which the AQEM value in FWE has changed recently
    taxa2$taxon[i]<-"SPHAERIUM CORNEUM/OVALE"
    taxa2$aqem2[i]<-22427
    taxa2$gen[i]<-"SPHAERIUM"
    taxa2$sp[i]<-"CORNEUM/OVALE"
  }
  if (taxa2$name[i]=="CYPHON SP.") {
    taxa2$taxon[i]<-"CYPHON SP. AD."
    taxa2$aqem2[i]<-14016
    taxa2$gen[i]<-"CYPHON"
    taxa2$sp[i]<-"SP. AD."
  }
}

##Check if any taxa are still missing codes, these will have to be added manually to the database
taxa2[which(taxa2$aqem2 %in% NA),]

##Check which AQEM values had to be altered, i.e., those in which taxa2$aqem2 is different from taxa2$aqem
sub<-taxa2[which(match(taxa2$aqem, taxa2$aqem2) %in% NA),]
View(sub)#can open sub to look at these entries

##Create empty columns to hold the trait values
taxa2[, colnames(dict)[6:length(colnames(dict))]]<-NA

##Create an empty dataframe to track the level at which each trait value was assigned
#1 is used when values are pulled directly from the dictionary, 0 is used when values are produced by averaging
cover<-taxa2[,c(3:4)]
for (i in 1:length(traits)) {
  cover[, (paste(traits[i], c("aqem", "subspecies", "genus", "family"), sep="_"))]<-NA
}


###---AQEM level

##Check for AQEM-level trait values first, this includes the species-level data if available
count<-0.1 #used to track progress through the loop
for (i in 1:length(taxa2$aqem2)) {
  sub<-dict[which(dict$aqem==taxa2$aqem2[i]), 6:length(colnames(dict))] #find the associated AQEM code in the database
  #If you have an -Agg. or -Gr. taxon with little data then use the associated species entry, comment this IF statement out to skip this step 
  if (length(grep("-AGG\\.", taxa2$taxon[i]))>0 | length(grep("-GR\\.", taxa2$taxon[i]))>0) {
    if (sum(unlist(sub), na.rm=TRUE)<11) { #if the AQEM code has little data then do this, many have only feeding traits, i.e., sum=10
      name1<-str_remove(taxa2$taxon[i], "-AGG.") #remove -Agg. text
      name1<-str_remove(name1, "-GR.") #remove -Gr. text
      sub<-dict[grep(name1, dict$name)[1],] #search for the name that remains and take the first entry
    }
  }
  for (j in 1:length(traits)) { #for each trait
    sub2<-sub[,grep(traits[j], colnames(sub))] #get the associated trait values
    if (sum(unlist(sub2), na.rm=TRUE)>0) { #if data is available
      taxa2[i , grep(traits[j], colnames(taxa2))]<-as.numeric(sub2) #then put that data in
      cover[i, grep(traits[j], colnames(cover))][1]<-1 #and set the coverage value to 1
    }
  }
  if (i/length(taxa2$aqem2)>=count) { #report progress
    cat(paste(count*100,"% ", sep=""))
    count<-count+0.1
  } 
}

###---Subspecies level

##Use subspecies averages for each subspecies if they have no data but others from the same subspecies do
nums<-numeric()
count<-1

for (i in 1:length(taxa2$aqem2)) {
  sub<-strsplit(taxa2$sp[i], " ")
  if (length(sub[[1]])>1 & sub[[1]][1] != "GEN." & sub[[1]][2] != "AD." & sub[[1]][2] != "LV." & length(grep("\\(", sub[[1]][1]))==0 & length(grep("\"", sub[[1]][1]))==0 & length(grep("/", sub[[1]][1]))==0 & taxa2$aqem[i] != 7456) {
    #this IF statement isolates entries with more than 2 names, but excludes those that have more than 2 because they are Gen. sp., Ad., Lv., multiple species (/), etc
    #also excludes Rhyacophila s. str. sp. (AQEM 7456) to ensure its values are assigned at the genus level
    nums[count]<-i
    count<-count+1
  }
}
#nums holds the positions of all subspecies in the dataset

##These are the subspecies in your original taxalist, check to make sure these are indeed subspecies
taxa2$name[nums]

for (i in nums) { #for each subspecies
  sub<-dict[which(dict$aqem==taxa2$aqem2[i]),] #get the associated AQEM entry from the dictionary
  for (j in 1:length(traits)) { #for each trait
    if (sum(taxa2[i , grep(traits[j], colnames(taxa2))], na.rm=TRUE)==0) { #if the trait has no values
      sub1<-grep(taxa2$gen[i], dict$name1) #then search for the genus
      sub2<-grep(strsplit(taxa2$sp[i], split=" ")[[1]][1], dict$name2) #and the first part of the subspecies name
      sub3<-dict[sub1[which(sub1 %in% sub2)], grep(traits[j], colnames(dict))] #take only members of the genus that have the first subspecies name
      if (sum(unlist(sub3), na.rm=TRUE)>0) { #if trait values are available
        taxa2[i , grep(traits[j], colnames(taxa2))]<-colMeans(sub3, na.rm=TRUE) #then use averaged values across the subspecies
        cover[i, grep(traits[j], colnames(cover))][2]<-0 #store a 0 for the cover value since averages have been used
      }
    }
  }
}


###---Genus level

##Genus - if still missing data get genus level data

count<-0.1 #used to track progress through the loop
for (i in 1:length(taxa2$aqem2)) { #for each taxon
  for (j in 1:length(traits)) { #for each trait
    if (sum(taxa2[i, grep(traits[j], colnames(taxa2))], na.rm=TRUE)==0) { #if it has no trait data
      sub<-dict[which(dict$aqem==taxa2$aqem2[i]),]
      if (length(grep("GEN. SP\\.", sub$name2))==0) { #if it has a genus (i.e., not family or subfamily)
        sub1<-subset(dict, dict$name1==sub$name1) #get all data for that genus
        nums1<-which(sub1$name2 == "SP.") #find if there is a genus-level entry
        if (length(grep("AD\\.", sub$name2[i]))>0) { #if it has an Adult specification
          sub1<-sub1[grep("AD\\.", sub1$name2),] #then get only the Adult entries
          nums1<-which(sub1$name2 == "SP. AD.")
        }
        if (length(grep("LV\\.", sub$name2[i]))>0) { #if it has a Larvae specification
          sub1<-sub1[grep("LV\\.", sub1$name2),] #then get only the Larvae entries
          nums1<-which(sub1$name2 == "SP. LV.")
        }
        nums2<-grep("SP\\.", sub1$name2) #find all the species-level entries
        sub2<-sub1[nums1, grep(traits[j], colnames(sub1))] #get the genus-level entry
        sub3<-sub1[-nums2, grep(traits[j], colnames(sub1))] #get only species-level entries for that genus
        if (length(nums1)==0) { #if there is no genus-level entry
          sub3<-sub1[, grep(traits[j], colnames(sub1))] #then sub3 is just sub1
        }
        if (sum(unlist(sub2), na.rm=TRUE)>0) { #if there is genus-level data
          taxa2[i , grep(traits[j], colnames(taxa2))]<-colMeans(sub2, na.rm=TRUE) #use those values, colMeans is used to deal with instances of multiple genus-level entries, generally only one of which has data
          cover[i, grep(traits[j], colnames(cover))][3]<-1 #and set the cover value to 1
        }
        if (sum(unlist(sub2), na.rm=TRUE)==0) { #if no genus-level data is available
          if (sum(unlist(sub3), na.rm=TRUE)>0) { #and if species-level values are available
              taxa2[i , grep(traits[j], colnames(taxa2))]<-colMeans(sub3, na.rm=TRUE) #then average across all species within the genus
              cover[i, grep(traits[j], colnames(cover))][3]<-0 #and set the cover value to 0
          }
        }
      }
    }
  }
  if (i/length(taxa2$aqem2)>=count) { #report progress
    cat(paste(count*100,"% ", sep=""))
    count<-count+0.1
  } 
}


###---Family level

count<-0.1 #used to track progress through the loop
for (i in 1:length(taxa2$aqem2)) { #for each taxon
  for (j in 1:length(traits)) { #for each trait
    if (sum(taxa2[i, grep(traits[j], colnames(taxa2))], na.rm=TRUE)==0) { #if it has no trait data
      sub<-dict[which(dict$aqem==taxa2$aqem2[i]),]
      if (sub$fam %nin% NA) { #if it has a family
        sub1<-subset(dict, dict$fam==sub$fam) #get all data for the associated family
        nums1<-grep(paste(str_to_title(sub$fam), "GEN. SP.", sep=" "), sub1$name) #find if there is a family-level entry
        if (length(grep("AD\\.", sub$name2[i]))>0) { #if it has an Adult specification
          sub1<-sub1[grep("AD\\.", sub1$name2),] #then get only the Adult entries
          nums1<-grep(paste(str_to_title(sub$fam), "GEN. SP. AD.", sep=" "), sub1$name) #find if there is a family-level entry
        }
        if (length(grep("LV\\.", sub$name2[i]))>0) { #if it has a Larvae specification
          sub1<-sub1[grep("LV\\.", sub1$name2),] #then get only the Larvae entries
          nums1<-grep(paste(str_to_title(sub$fam), "GEN. SP. LV.", sep=" "), sub1$name) #find if there is a family-level entry
        }
        sub2<-sub1[nums1, grep(traits[j], colnames(sub1))] #get the family-level entry
        sub3<-sub1[-nums1, grep(traits[j], colnames(sub1))] #get all non-family-level for that family
        if (length(nums1)==0) {
          sub3<-sub1[, grep(traits[j], colnames(sub1))]
        }
        if (sum(unlist(sub2), na.rm=TRUE)>0) { #if there is family-level data
          taxa2[i , grep(traits[j], colnames(taxa2))]<-colMeans(sub2, na.rm=TRUE) #use those values
          cover[i, grep(traits[j], colnames(cover))][4]<-1 #and set the cover value to 1
        }
        if (sum(unlist(sub2), na.rm=TRUE)==0) { #if there is no family-level data
          if (sum(unlist(sub3), na.rm=TRUE)>0) { #and if taxa values are available for the family
            taxa2[i , grep(traits[j], colnames(taxa2))]<-colMeans(sub3, na.rm=TRUE) #then average across all taxa within the family
            cover[i, grep(traits[j], colnames(cover))][4]<-0 #and set the cover value to 0
          }
        }
      }
    }
  }
  if (i/length(taxa2$aqem2)>=count) { #report progress
    cat(paste(count*100,"% ", sep=""))
    count<-count+0.1
  } 
}


##--Percent and total coverage for each trait at each level
cover<-rbind(cover, NA, NA, NA)
cover$taxon[length(cover$taxon)-2]<-"Percent using AQEM"
cover$taxon[length(cover$taxon)-1]<-"Percent using Averages"
cover$taxon[length(cover$taxon)]<-"Total percent"

for (i in 3:length(colnames(cover))) {
  cover[length(cover$taxon)-2, i]<-length(which(cover[(1:(length(cover$taxon)-3)),i]==1))/length(taxa2$aqem2)
  cover[length(cover$taxon)-1, i]<-length(which(cover[(1:(length(cover$taxon)-3)),i]==0))/length(taxa2$aqem2)
  cover[length(cover$taxon), i]<-cover[length(cover$taxon)-2, i]+cover[length(cover$taxon)-1, i]
}

##--Save the trait and coverage data to your working directory
write.csv(taxa2, "Outputs/Trait database.csv", row.names=FALSE)
write.csv(cover, "Outputs/Trait coverage.csv", row.names=FALSE)

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
