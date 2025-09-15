#Attempting Geange Analysis on the pools data, utilizing Julian_Day_Wrap as 
#a continuous variable and pool number as a categorical variable 

#SECTION 1: DATA PREP, IMPORT, AND CLEANING-------------------------------------

#import packages and data for analysis

#Install packages
install.packages('tidyverse')
install.packages('vegan')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('ggridges') 
install.packages('readxl')
remove.packages("hrbrthemes")
remove.packages("systemfonts")
install.packages("hrbrthemes")
install.packages("systemfonts")
install.packages("gridExtra")
install.packages("grid")
install.packages("abind")
install.packages("here")

#call packages
library(tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(ggridges) 
library(readxl) 
library(systemfonts)
#library(hrbrthemes)
library(gridExtra)
library(grid)
library(abind) 
library(here)

#import and prep data 
#ensure wd is correct
setwd(here())
#import data
WideData <- read_excel(here("PoolsData_Final.xlsx"))
#flip data
LongData <- gather(data  = WideData, 
                   key   = "Species", 
                   value = "Count", "Ostracod":"Springtail") 


#turn NA strings into real R na values
LongData <- LongData |> 
  mutate(Count            = na_if(Count, "NA")) |> 
  mutate(Depth_cm         = na_if(Depth_cm, "NA")) |> 
  mutate(Temperature_C    = na_if(Temperature_C, "NA")) |> 
  mutate(Dissolved_Oxygen = na_if(Dissolved_Oxygen, "NA"))

#check to make sure there aren't any NAs remaining in other columns
LongDataCheck1 <- LongData[, !(names(LongData) %in% c("Count", "Depth_cm", "Temperature_C", "Dissolved_Oxygen", "NOTES"))]
sum(is.na(LongDataCheck1))
#we good if you see 0

#turn Count column into numeric
LongData$Count <- as.numeric(LongData$Count)
#note there are now going to be some NAs where data was missing 

#remove any sampling sequence other than the first of each month 
#define sequences to remove
other_sequences = c("2nd Feb", "2nd Mar", "3rd Mar", "4th Mar", "2nd Apr", "2nd May")
#use filter function to remove them
LongData <- LongData |>  
  filter(!(Sampling_Sequence %in% other_sequences))
#same process to remove pond data 
LongData <- LongData |> 
  filter(!(Pool_Number %in% "Pond"))

#Expand counts, so each observation is its own line
  #convert count to integer 
  LongData$Count <- as.integer(LongData$Count)
  #convert NAs in count  
  LongData$Count[is.na(LongData$Count)] <- 0
  #expand counts 
  LongData <- uncount(LongData, weights = Count)

#Load functions from Geange files
source(here("ExampleCfiles", "MEE3_070_sm_NicheFunctions.txt"))

#add id column and move it and species to the front
LongData <- LongData |> 
  mutate(id = paste0("id", row_number())) |> 
  select(id, Species, everything())

#START GEANGE ANALYSIS----------------------------------------------------------

##JULIAN DAY (COUNTINUOUS VARIABLE)---------------------------------------------

### get data prepped -----------------------------------------------------------

# read in the individual data file
B.df <- select(LongData, id, Species, Julian_Day)


# Ensure the first two column names are "id" and "species".
colnames(B.df)[1] <- "id"
colnames(B.df)[2] <- "species"

# Ensure that the first 2 cols are factors.
B.df$id      <- as.factor(B.df$id)
B.df$species <- as.factor(B.df$species)

# Store some vectors of names:
spnames   <- sort(unique(as.character(B.df$species)))
no.spp    <- length(spnames)

varnames <- colnames(B.df)[-(1:2)]    
no.vars  <- length(varnames)  

#make a vector of variable types 
vartypes <- c("cts") 
#check they are correctly labeled: 
cbind(varnames,vartypes)

# Set up a list of objects which are NULL if this is not
# a resource selection variable, and with the availability
# vector if it is resource selection.
avail.list        <- vector("list",no.vars)
names(avail.list) <- varnames
avail.list

### Set up R objects to store results ------------------------------------------

# alpha.list

# The object alpha.list has one component per variable.
# The components are NULL for ordinary variables.

Balpha.list        <- vector("list",no.vars)
names(Balpha.list) <- varnames

for (vv in 1:no.vars) if (vartypes[vv]=="rsel")
{
  choices                     <- unique(B.df[,vv+2])
  no.ch                       <- length(choices)
  Balpha.list[[vv]]           <- matrix(NA,no.spp,no.ch)
  dimnames(Balpha.list[[vv]]) <- list(spnames,choices)
}

# no.array

# Set up an array of niche overlaps.
# The object no.array is an array of niche overlaps.
# It is a 3-D array, with rows and columns being species 
# (a square symmetric matrix for pairwise niche overlaps), 
# and the layers are the dimensions for the multivariate 
# niche overlap measure (one dimension per variable).
# Rows and columns are species, layers are variables.


Bno.array           <- array(1,c(no.spp,no.spp,no.vars))
dimnames(Bno.array) <- list(spnames,spnames,varnames) 

# Run through each variable in turn, identify its type,
# calculate the appropriate NO matrix and store it in
# the right layer of the no.array. 

for (vv in 1:no.vars)
{
  #slight change to the Geange code here: 
  y <- B.df[[varnames[vv]]]
  #this ensures that y is a vector, and not a 1 column tibble 
  #the latter happened with my data and not the example dataset 
  #no idea why but this seems to work 
  
  #adding prints in here for trouble shooting:
  print(paste("vv =", vv))
  print(str(y))
  print(paste("vartype =", vartypes[vv])) 
  
  if (vartypes[vv] == "bin")
    Bno.array[,,vv] <- no.bin.fn(B.df$species,y)
  if (vartypes[vv] == "cat")
    Bno.array[,,vv] <- no.cat.fn(B.df$species,y)
  if (vartypes[vv] == "count")
    Bno.array[,,vv] <- no.count.fn(B.df$species,y)
  if (vartypes[vv] == "cts")
    Bno.array[,,vv] <- no.cts.fn(B.df$species,y)
  if (vartypes[vv] == "meas")
    Bno.array[,,vv] <- no.cts.fn(B.df$species,log(y))
  if (vartypes[vv] == "pcent")
    Bno.array[,,vv] <- no.cts.fn(B.df$species,
                                 log(y/(100-y)))
  if (vartypes[vv] == "propn")
    Bno.array[,,vv] <- no.cts.fn(B.df$species,
                                 log(y/(1-y)))
  if (vartypes[vv] == "rsel")
  {
    
    # Do Manly's alpha calculations, store.
    avail.vect       <- avail.list[[vv]]
    alpha.mat        <- alpha.fn(B.df$species,y,avail.vect)
    alpha.list[[vv]] <- alpha.mat         
    
    # Do niche overlaps, as proportions in categories:
    Bno.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
  }
}

### Permutation Testing -------------------------------------------------------- 
# Permutation of the species labels would give data 
# satisfying the null model of complete niche overlap, 
# i.e. that none of the variables 
# serves to differentiate species into different niches.

# Hence for each replication, permute the species labels
# and run through all the calculations above.
# Stor NOs in an array with one extra dimension, one
# layer for each replication.
# Then the null distributions are all stored.
# Can use the original availability data, but need a new 
# alpha list each time.

# Choose no. of replications.
# Start low, eg. with 10 reps, to check it is working.
# Then do more reps, e.g. 1000 reps for 3 decimal places in p-values.
replic <- 10

Bpseudo.no.array <- array(1,c(no.spp,no.spp,no.vars,replic))
dimnames(Bpseudo.no.array) <- list(spnames,spnames,varnames,NULL)

# Set a temporary data frame, which will change each time
# through the cycle by having its species column permuted.
temp.df <- B.df


# For each replication, permute the species labels, run the
# niche overlap calculations, and store the results in the
# pseudo NO array
for (rr in 1:replic)
{
  
  # Permute the species labels in the temporary dataframe:
  temp.df$species <- sample(temp.df$species)
  for (vv in 1:no.vars)
  {
    
    # Read out the column from this variable:
      #same change to Geange code as before:
    y <- temp.df[[varnames[vv]]]
    
    # Run through the variable types, do appropriate analyses:
    if (vartypes[vv] == "bin")
      Bpseudo.no.array[,,vv,rr] <- no.bin.fn(temp.df$species,y)
    if (vartypes[vv] == "cat")
      Bpseudo.no.array[,,vv,rr] <- no.cat.fn(temp.df$species,y)
    if (vartypes[vv] == "count")
      Bpseudo.no.array[,,vv,rr] <- no.count.fn(temp.df$species,y)
    if (vartypes[vv] == "cts")
      Bpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,y)
    if (vartypes[vv] == "meas")
      Bpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,log(y))
    if (vartypes[vv] == "pcent")
      Bpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(100-y)))
    if (vartypes[vv] == "propn")
      Bpseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(1-y)))
    if (vartypes[vv] == "rsel")
    {
      
      # Do Manly's alpha calculations, store.
      avail.vect <- avail.list[[vv]]
      alpha.mat  <- alpha.fn(temp.df$species,y,avail.vect)
      
      # Do niche overlaps, as proportions in categories:
      Bpseudo.no.array[,,vv,rr] <- no.rsel.cat.fn(alpha.mat)
    }
  }
  print(paste("Rep",rr,"done"))
}

## POOL #/CREEK (CATEGORICAL DATA)---------------------------------------------- 
