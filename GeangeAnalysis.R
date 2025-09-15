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
    install.packages("openxlsx")

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
    library(openxlsx)

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
    source(here("ExampleAfiles", "MEE3_070_sm_NicheFunctions.txt"))
  
  #add id column and move it and species to the front
    LongData <- LongData |> 
      mutate(id = paste0("id", row_number())) |> 
      select(id, Species, everything())

#SECTION 2: GEANGE ANALYSIS-----------------------------------------------------

##2.1: JULIAN DAY (COUNTINUOUS VARIABLE) AND POOL # (CATEGORICAL VARIABLE)------

### get data prepped -----------------------------------------------------------

# read in the individual data file
  LD.df <- select(LongData, id, Species, Julian_Day, Pool_Number)


# Ensure the first two column names are "id" and "species".
  colnames(LD.df)[1] <- "id"
  colnames(LD.df)[2] <- "species"

# Ensure that the first 2 cols are factors.
  LD.df$id      <- as.factor(LD.df$id)
  LD.df$species <- as.factor(LD.df$species)

# Store some vectors of names:
  spnames   <- sort(unique(as.character(LD.df$species)))
  no.spp    <- length(spnames)

  varnames <- colnames(LD.df)[-(1:2)]    
  no.vars  <- length(varnames)  

#make a vector of variable types 
  vartypes <- c("cts", "cat") 
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

  alpha.list        <- vector("list",no.vars)
  names(alpha.list) <- varnames
  
  for (vv in 1:no.vars) if (vartypes[vv] == "rsel")
  {
    choices                    <- unique(LD.df[,vv + 2])
    no.ch                      <- length(choices)
    alpha.list[[vv]]           <- matrix(NA,no.spp,no.ch)
    dimnames(alpha.list[[vv]]) <- list(spnames,choices)
  }

# no.array
  
  # Set up an array of niche overlaps.
  # The object no.array is an array of niche overlaps.
  # It is a 3-D array, with rows and columns being species 
  # (a square symmetric matrix for pairwise niche overlaps), 
  # and the layers are the dimensions for the multivariate 
  # niche overlap measure (one dimension per variable).
  # Rows and columns are species, layers are variables.
  
  
  no.array           <- array(1,c(no.spp,no.spp,no.vars))
  dimnames(no.array) <- list(spnames,spnames,varnames) 

# Run through each variable in turn, identify its type,
# calculate the appropriate NO matrix and store it in
# the right layer of the no.array. 
  
  for (vv in 1:no.vars)
  {
    #slight change to the Geange code here: 
    y <- LD.df[[varnames[vv]]]
    #this ensures that y is a vector, and not a 1 column tibble 
    #the latter happened with my data and not the example dataset 
    #no idea why but this seems to work 
    
    #adding prints in here for trouble shooting:
    print(paste("vv =", vv))
    print(str(y))
    print(paste("vartype =", vartypes[vv])) 
    if (vartypes[vv] == "bin")
      no.array[,,vv] <- no.bin.fn(LD.df$species,y)
    if (vartypes[vv] == "cat")
      no.array[,,vv] <- no.cat.fn(LD.df$species,y)
    if (vartypes[vv] == "count")
      no.array[,,vv] <- no.count.fn(LD.df$species,y)
    if (vartypes[vv] == "cts")
      no.array[,,vv] <- no.cts.fn(LD.df$species,y)
    if (vartypes[vv] == "meas")
      no.array[,,vv] <- no.cts.fn(LD.df$species,log(y))
    if (vartypes[vv] == "pcent")
      no.array[,,vv] <- no.cts.fn(LD.df$species,
                                  log(y/(100 - y)))
    if (vartypes[vv] == "propn")
      no.array[,,vv] <- no.cts.fn(LD.df$species,
                                  log(y/(1 - y)))
    if (vartypes[vv] == "rsel")
    {
      
      # Do Manly's alpha calculations, store.
      no.choices       <- length(avail.list[[vv]])
      choicenames      <- names(avail.list[[vv]])
      avail.vect       <- avail.list[[vv]]
      alpha.mat        <- alpha.fn(LD.df$species,y,avail.vect)
      alpha.list[[vv]] <- alpha.mat         
      
      # Do niche overlaps, as proportions in categories:
      no.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
    }
  }

#also calculate overall NO measures, averaged over dimensions
no.overall.mat    <- apply(no.array,c(1,2),mean)
no.overall.mat.sd <- apply(no.array,c(1,2),sd)


##2.2: PERMUTATION TESTING -----------------------------------------------------

### set up species label permutation -------------------------------------------

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

  pseudo.no.array           <- array(1,c(no.spp,no.spp,no.vars,replic))
  dimnames(pseudo.no.array) <- list(spnames,spnames,varnames,NULL)

# Set a temporary data frame, which will change each time
# through the cycle by having its species column permuted.
temp.df <- LD.df


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
        pseudo.no.array[,,vv,rr] <- no.bin.fn(temp.df$species,y)
      if (vartypes[vv] == "cat")
        pseudo.no.array[,,vv,rr] <- no.cat.fn(temp.df$species,y)
      if (vartypes[vv] == "count")
        pseudo.no.array[,,vv,rr] <- no.count.fn(temp.df$species,y)
      if (vartypes[vv] == "cts")
        pseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,y)
      if (vartypes[vv] == "meas")
        pseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,log(y))
      if (vartypes[vv] == "pcent")
        pseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                              log(y/(100 - y)))
      if (vartypes[vv] == "propn")
        pseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                              log(y/(1 - y)))
      if (vartypes[vv] == "rsel")
      {
        
        # Do Manly's alpha calculations, store.
        no.choices  <- length(avail.list[[vv]])
        choicenames <- names(avail.list[[vv]])
        avail.vect  <- avail.list[[vv]]
        alpha.mat   <- alpha.fn(temp.df$species,y,avail.vect)
        
        # Do niche overlaps, as proportions in categories:
        pseudo.no.array[,,vv,rr] <- no.rsel.cat.fn(alpha.mat)
      }
    }
    print(paste("Rep",rr,"done"))
  }


### null model analysis --------------------------------------------------------

# Calculate p values for each pair of species 
# separately for each variable.
  sep.pvals           <- array(1,c(no.spp,no.spp,no.vars))
  dimnames(sep.pvals) <- list(spnames,spnames,varnames)

  for (spa in 1:(no.spp - 1)) for (spb in (spa + 1):no.spp)
    for (vv in 1:no.vars)   
    {
      pseudo.nos            <- pseudo.no.array[spa,spb,vv,]
      data.no               <- no.array[spa,spb,vv]
      sep.pvals[spa,spb,vv] <- mean(pseudo.nos < data.no) 
      length(pseudo.nos[data.no < pseudo.nos])
      sep.pvals[spb,spa,vv] <- sep.pvals[spa,spb,vv] 
    }

# Also find p value for overall NO measure.
  overall.pvals           <- matrix(1,no.spp,no.spp)
  dimnames(overall.pvals) <- list(spnames,spnames)

  for (spa in 1:(no.spp - 1)) for (spb in (spa + 1):no.spp)
  {
    temp.mat               <- pseudo.no.array[spa,spb,,]
    pseudo.nos             <- apply(temp.mat,2,mean)
    data.no                <- no.overall.mat[spa,spb]
    overall.pvals[spa,spb] <- mean(pseudo.nos < data.no) 
    length(pseudo.nos[data.no < pseudo.nos])
    overall.pvals[spb,spa] <- overall.pvals[spa,spb] 
  }

#Null model analysis to determine if distribution of species is more 
#differentiated or more clustered than expected 

#reformat observed data to derive matrix of niche overlaps with one row per 
#species, and one column for each niche dimension 
  VV <- ncol(LD.df[,-c(1:2)])
  RR <- replic   # Number of replications.

  #making a slight adjustment to Geange code because it doesn't work with the  
  #number of dimensions our array has 
  
  no.pairs <- no.spp * (no.spp - 1) / 2
  no.mat   <- matrix(NA, nrow = no.pairs, ncol = VV)
  
  for (vv in 1:VV)
    no.mat[, vv] <- as.vector(as.dist(no.array[, , vv]))

# Next, reformat the pseudo data to derive a matrix of niche overlaps
# with one row per species, and one column for each niche dimension,
# with one extra dimension, one layer for each replication
  
  #applying the same adjustment here as before because of dimensional issues
  pseudo.mat <- array(NA, dim = c(no.pairs, VV, RR))
  
  for (vv in 1:VV) for (rr in 1:RR) {
    pseudo.mat[, vv, rr] <- as.vector(as.dist(pseudo.no.array[,, vv, rr]))
  }

# For each niche dimension, calculate mean and variance over the species
# pairs, and hence the test statistic ch = coefficient of heterogeneity.
# Note: Need to use variance formula based on n, not n-1.
  
  KK <- ncol(no.mat)      # Number of niche dimensions
  SS <- nrow(no.mat)      # Number of species pairs
  RR <- replic            # Number of replications.
  
  
  data.ch   <- rep(NA,KK)
  pseudo.ch <- matrix(NA,RR,KK)
  
  for (kk in 1:KK)
  {
    # Calculate data test statistic:
    x <- mean(no.mat[,kk])
    v <- var(no.mat[,kk])*(SS - 1)/SS # Adjust for denom n, not n-1
    data.ch[kk] <- v/x/(1 - x)
    
    # Calculate test stats for all pseudo-data:
    for (rr in 1:RR)
    {
      x <- mean(pseudo.mat[,kk,rr])
      v <- var(pseudo.mat[,kk,rr])*(SS - 1)/SS
      pseudo.ch[rr,kk] <- v/x/(1 - x)
    }
  }

# For each niche dimension, see if data more differentiated than random.
  p.dims.diff <- rep(NA,KK)
  for (kk in 1:KK)
    p.dims.diff[kk]  <- mean(data.ch[kk] > pseudo.ch[,kk])
  names(p.dims.diff) <- paste("diff.dim",sort(varnames))

# For each niche dimension, see if data more clustered than random.
  p.dims.clus <- rep(NA,KK)
  for (kk in 1:KK)
    p.dims.clus[kk]  <- mean(data.ch[kk] < pseudo.ch[,kk])
  names(p.dims.clus) <- paste("clus.dim",sort(varnames))
  
# For average niche overlap, calculate mean and variance over the species
# pairs, and hence the test statistic ch = coefficient of heterogeneity.
# Note: Need to use variance formula based on n, not n-1.
  
  overall.data.ch   <- mean(data.ch)
  overall.pseudo.ch <- apply(pseudo.ch,1,mean)
  
  # Test if this community is more differentiated than random:
  p.all.diff <- mean(overall.data.ch > overall.pseudo.ch)
  
  # Test if this community is more clustered than random:
  p.all.clus <- mean(overall.data.ch < overall.pseudo.ch)

#save results! 
  NO.results <- list(
    info = list(variables = cbind(varnames,vartypes),
                perm.reps = replic),
    NOestimates = no.array,
    separate.pvalues = sep.pvals,
    separate.cluster.pvalues = p.dims.clus,
    separate.differentiated.pvalues = p.dims.diff,
    ests.overall = no.overall.mat,
    ests.overall.sd = no.overall.mat.sd,
    overall.pvalues = overall.pvals,
    overall.cluster.pvalues = p.all.clus,
    overall.differentiated.pvalues = p.all.diff)

#save results as xlsx file

#save NOestimates as a xlsx with 2 sheets, one for each niche variable
  wb <- createWorkbook()
  
  for (i in 1:dim(NO.results$NOestimates)[3]) {
    df         <- as.data.frame(NO.results$NOestimates[, , i])
    sheet_name <- paste0("Layer_", i)
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, df, rowNames = TRUE)
  }
  
  # Save workbook
  saveWorkbook(wb, file = "NOestimates.xlsx", overwrite = TRUE)

#save ests.overall and overall.pvalues as one xlsx
  wb <- createWorkbook()
  
  # Add first worksheet and write ests.overall
    addWorksheet(wb, "Ests_Overall")
    writeData(wb, sheet = "Ests_Overall", as.data.frame(NO.results$ests.overall), rowNames = TRUE)
    
  # Add second worksheet and write matrix2
    addWorksheet(wb, "p_values")
    writeData(wb, sheet = "p_values", as.data.frame(NO.results$overall.pvalues), rowNames = TRUE)
    
  # Save workbook to a file
    saveWorkbook(wb, file = "EstsAndpValues.xlsx", overwrite = TRUE)
  