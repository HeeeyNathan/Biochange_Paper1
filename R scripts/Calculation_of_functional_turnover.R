##### LOAD LIBRARIES --------------------
# Load the required packages
library(FD) # functional diversity 
library(pacman) # package manager
library(vegan) # community ecology, taxonomic diversity
library(BAT) 

##### LOAD THE DATA --------------------
#### TAXONOMIC DATA
### Full site-by-species matrix for taxonomic analysis
comm_raw <- read.csv("Data/LT_taxalist_2010-2020_short_transposed.csv", header = TRUE) # Add "LT_taxalist_2010-2020_short_final_transposed.csv"
comm_raw <- comm_raw[, -21] # remove 'Ametropus_sp.' because it has no trait info
comm <- comm_raw[, -1] # remove some columns that are needed now
rownames(comm) <- comm_raw$sample_id

#### FUNCITONAL DATA
### Full species-by-trait matrix for functional analysis
traits_raw <- read.csv("Data/LT_traits_2010-2020_short.csv", header = TRUE) # Add "LT_taxalist_2010-2020_tachet_bio_traits.csv"
traits_raw <- traits_raw[-20, ] # remove 'Ametropus_sp.' because it has no trait info
traits <- traits_raw[, -c(1:5)] # remove some columns that are needed now
rownames(traits) <- traits_raw$taxon_name # make row names the taxon names

## Changes trait data to proportions (between 0 & 1)
# Trait category blocks
traits <- prep.fuzzy.var(traits, c(4, 4, 8, 9, 2, 8, 7, 3, 8, 5, 5)) # These numbers show the number of traits in each trait group
rowSums(traits) #rowsums should sum to the number of trait groups you included

### functional site-by-species matrix
identical(row.names(traits), colnames(comm)) # check if names are the same

##### FUNCTIONAL TURNOVER --------------------
# Load required packages
library(codyn) #for functional turnover - stops turnover from becoming masked by base turnover function

# species-by-trait matrix
comm_traits <- traits
#### NOTE: all these files must be in order by year (smallest to largest) within site
# Site-by-species matrix
comm_abund <- comm
#Site-Year factors
facto <- read.csv("Data/LT_siteYr_Factors.csv", h = T, sep = ",", stringsAsFactors = FALSE)
head(facto)
#traits long form
cwm_l <- read.csv("Data/LT_CWMs_long.csv", h = T, sep = ",", stringsAsFactors = FALSE)
head(cwm_l)

############################# 1 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR1", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur_o <- cbind(fact_site, FC_biol_ord)

############################# 2 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR105", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 3 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR11", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 4 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR127", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 5 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR1282", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 6 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR13", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 7 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR1301", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 8 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR1319", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 9 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR133", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 10 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR1348", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 11 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR137", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 12 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR1438", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 13 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR1462", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 14 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR161", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 15 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR175", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 16 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR192", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 17 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR231", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 18 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR245", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 19 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR26", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 20 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR268", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 21 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR325", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 22 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR327", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 23 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR33", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 24 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR387", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 25 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR40", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 26 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR401", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 27 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR41", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 28 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR43", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 29 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR450", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 30 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR50", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 31 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR612", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 32 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR62", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 33 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR65", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 34 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR70", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 35 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR78", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 36 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR787", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 37 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR79", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 38 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR82", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 39 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR86", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 40 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR88", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

############################# 41 ##################################
# Create site subset of community abundances
site_sub <- comm_abund[facto$site_id == "LTR92", ]
site_s <- site_sub[, colSums(site_sub) > 0] # only keeps columns with column sums > 0

# Create site subset of traits
spp_trait = (row.names(comm_traits) %in% colnames(site_s)) #make a list of T or F is there is a match in spp
traits_site = comm_traits[spp_trait, ] #subset by this list

# Create site subset of CWMs long
spp_cwm_long = (cwm_l$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
cwms_site_long = cwm_l[spp_cwm_long, ] #subset by this list

# Create site subset of factor
spp_fact = (facto$sample_id %in% row.names(site_s)) #make a list of T or F is there is a match in spp
fact_site = facto[spp_fact, ] #subset by this list

# Check if names are the same
identical(row.names(traits_site), colnames(site_s))

# Functional Traits
# Principal coordinates analysis (PCoA) to return PCoA axes 
FC_biol_ord <- dbFD(traits_site, site_s, calc.FRic = T, m = 6, stand.FRic = T, 
                    scale.RaoQ = T, calc.FDiv = T, asym.bin = 5:10, calc.CWM = FALSE)

# Functional turnover 
# Calculate relative total turnover within replicates
s1_tot_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "total") #replicate.var = "trait",
fact_site$F_to <- c("NA", s1_tot_FuncTu$total) # Turnover per yr And adds this col to fact_site

# Calculate relative species appearances within replicates
s1_app_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "appearance")
fact_site$F_to_app <- c("NA", s1_app_FuncTu$appearance) 

# Calculate relative species disappearances within replicates
s1_dis_FuncTu <- turnover(cwms_site_long, time.var = "year", species.var = "modality", abundance.var = "CWM", metric = "disappearance")
fact_site$F_to_disap <- c("NA", s1_dis_FuncTu$disappearance) 

## Add results to main dataset
tur <- cbind(fact_site, FC_biol_ord)
tur_o <- rbind(tur, tur_o)

#### Save output
write.csv(tur_o, "Outputs/LT_siteYr_Functurnover.csv") # keep row names for later ordering

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