#This analysis iterates on GeangeAnalysis01, adding temperature, oxygen, and
#recruitment amount as niche axes. 
#it will also differentiate between pool and creek data 

#SECTION 1: DATA PREP, IMPORT, AND CLEANING ------------------------------------
  
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
    install.packages("tictoc")
  
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
    library(tictoc)
  
  #import and prep data 
  #ensure wd is correct
    setwd(here())
  #import data
    WideData <- read_excel(here("PoolsDataWReg.xlsx"))
  #flip data
    LongData <- gather(data  = WideData, 
                      key   = "Species", 
                      value = "Count", "Ostracod":"Springtail") 
  
  #turn NA strings into real R na values
   LongData <- LongData |> 
      mutate(Count            = na_if(Count, "NA")) |> 
      mutate(Depth_cm         = na_if(Depth_cm, "NA")) |> 
      mutate(Temperature_C    = na_if(Temp_Reg, "NA")) |> 
      mutate(Dissolved_Oxygen = na_if(O2_Reg, "NA"))
  
  #check to make sure there aren't any NAs remaining in other columns
    sum(is.na(LongData[, !(names(LongData) %in% c("Count", "Depth_cm", 
                                                                    "Temperature_C", 
                                                                    "Dissolved_Oxygen",
                                                                    "NOTES"))]))
    #we good if you see 0
  
  #turn Count column into numeric
    LongData$Count <- as.numeric(LongData$Count)
   #note there are now going to be some NAs where data was missing 
  
  #remove any sampling sequence other than the first of each month 
  #define sequences to remove
    other_sequences = c("2nd Feb", "2nd Mar", "3rd Mar", "4th Mar", "2nd Apr", 
                        "2nd May")
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
  
 #make new LongData with just the pools (removing the creek) 
   LongDataPools <- LongData |> 
     filter(!(Pool_Number %in% "Creek"))
      
      
  #remove uncommon groups from analysis (uncommon = prevelance under 0.5%)
      #summarise data and sort by numerical prevelance
      LD_counts <- LongDataPools |> 
        group_by(Species) |> 
        summarise(count = n()) |> 
        mutate(percentage = (count/sum(count))*100)
      
      #7 groups with percentage > 1%
      #10 groups with percentage > 0.5%
      #18 groups with percentage > 0.1%
      
      species05 <- c("Cladoceran", "Copepod", "Diptera_Pupae", 
                     "Midge_Larvae", "Mites", "Mosquito_Larvae", "Ostracod", 
                     "Roundworm", "Springtail") 
      
      LongDataPools <- LongDataPools |> 
        filter(Species %in% species05)
      
  #add id column and move it and species to the front
    LongDataPools <- LongDataPools |> 
      mutate(id = paste0("id", row_number())) |> 
      select(id, Species, everything())
    
  
  #convert Depth_cm NAs into 0
    LongDataPools$Depth_cm[is.na(LongDataPools$Depth_cm)] <- 0

  #convert columns into number 
    LongDataPools$Depth_cm <- as.numeric(LongDataPools$Depth_cm)
    LongDataPools$Temp_Reg <- as.numeric(LongDataPools$Temp_Reg)
    LongDataPools$O2_Reg <- as.numeric(LongDataPools$O2_Reg)

#SECTION 2: (ALL TOGETHER) GEANGE ANALYSIS WITHOUT TEMP AND O2------------------
  #this section performs Genage analysis on all the data, excluding temp and O2 
  #data and not separating by year
    
##2.1: DATA PREP AND SETUP -----------------------------------------------------
  
    # read in the individual data file
    LD.df <- select(LongDataPools, id, Species, Julian_Day, Pool_Number, Recruitment_Amount, Depth_cm)
    
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
    vartypes <- c("cts", "cat", "cts", "cts") 
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
    replic <- 1000
    
    pseudo.no.array           <- array(1,c(no.spp,no.spp,no.vars,replic))
    dimnames(pseudo.no.array) <- list(spnames,spnames,varnames,NULL)
    
    # Set a temporary data frame, which will change each time
    # through the cycle by having its species column permuted.
    temp.df <- LD.df
    
    
    # For each replication, permute the species labels, run the
    # niche overlap calculations, and store the results in the
    # pseudo NO array
    
    ##WARNING: 1000 REPLICATES TAKES ABOUT 11 MINUTES. REDUCE WHEN TESTING
    
    tic()
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
    toc()
    
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
  
##2.3: SAVE RESULTS ------------------------------------------------------------
    #save results! 
    NO_all_4axes.results <- list(
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
    

#SECTION 3: (ALL TOGETHER) GEANGE ANALYSIS WITH TEMP AND O2 --------------------
    #this section performs Genage analysis on all the data, including temp and O2 
    #data and not separating by year
    
    ##3.1: DATA PREP AND SETUP -------------------------------------------------
    
    # read in the individual data file
    LD.df <- select(LongDataPools, id, Species, Julian_Day, Pool_Number, 
                    Recruitment_Amount, Depth_cm, Temp_Reg, O2_Reg)
    
    #remove all rows where temp and 02 data are missing
    LD.df <- LD.df[complete.cases(LD.df), ]
    
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
    vartypes <- c("cts", "cat", "cts", "cts", "cts", "cts") 
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
    
    ##3.2: PERMUTATION TESTING -----------------------------------------------------
    
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
    replic <- 1000
    
    pseudo.no.array           <- array(1,c(no.spp,no.spp,no.vars,replic))
    dimnames(pseudo.no.array) <- list(spnames,spnames,varnames,NULL)
    
    # Set a temporary data frame, which will change each time
    # through the cycle by having its species column permuted.
    temp.df <- LD.df
    
    
    # For each replication, permute the species labels, run the
    # niche overlap calculations, and store the results in the
    # pseudo NO array
    tic()
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
    toc()
    
    ### null model analysis ----------------------------------------------------
    
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
    
##3.3: SAVE RESULTS ------------------------------------------------------------
    #save results! 
    NO_all_6axes.results <- list(
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

#SECTION 4: (YEARS SEPARATE) GEANGE ANALYSIS EXCLUDING TEMP AND O2--------------
  #this section does identical analysis to section 2 but separates data by year
  
##4.1: DATA PREP AND GEANGE SETUP ----------------------------------------------
    
    # read in the individual data file
    LD.df <- select(LongDataPools, id, Species, Julian_Day, Pool_Number, 
                    Recruitment_Amount, Depth_cm) 
    
    
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
    vartypes <- c("cts", "cat", "cts", "cts") 
    #check they are correctly labeled: 
    cbind(varnames,vartypes)
    
    # Set up a list of objects which are NULL if this is not
    # a resource selection variable, and with the availability
    # vector if it is resource selection.
    avail.list        <- vector("list",no.vars)
    names(avail.list) <- varnames
    avail.list
    
    #separate LD.df into 2 files, one for year 0 and one for year 1
      LD0.df <- LD.df[LD.df$Julian_Day < 365, ] 
      LD1.df <- LD.df[LD.df$Julian_Day > 366, ]
    
      
##4.2: YEAR 0 GEANGE ANALYSIS --------------------------------------------------  
    ### Set up R objects to store results --------------------------------------
      
    # alpha.list
    
    # The object alpha.list has one component per variable.
    # The components are NULL for ordinary variables.
    
    alpha.list        <- vector("list",no.vars)
    names(alpha.list) <- varnames
    
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
      y <- LD0.df[[varnames[vv]]]
      #this ensures that y is a vector, and not a 1 column tibble 
      #the latter happened with my data and not the example dataset 
      #no idea why but this seems to work 
      
      #adding prints in here for trouble shooting:
      print(paste("vv =", vv))
      print(str(y))
      print(paste("vartype =", vartypes[vv])) 
      if (vartypes[vv] == "bin")
        no.array[,,vv] <- no.bin.fn(LD0.df$species,y)
      if (vartypes[vv] == "cat")
        no.array[,,vv] <- no.cat.fn(LD0.df$species,y)
      if (vartypes[vv] == "count")
        no.array[,,vv] <- no.count.fn(LD0.df$species,y)
      if (vartypes[vv] == "cts")
        no.array[,,vv] <- no.cts.fn(LD0.df$species,y)
      if (vartypes[vv] == "meas")
        no.array[,,vv] <- no.cts.fn(LD0.df$species,log(y))
      if (vartypes[vv] == "pcent")
        no.array[,,vv] <- no.cts.fn(LD0.df$species,
                                    log(y/(100 - y)))
      if (vartypes[vv] == "propn")
        no.array[,,vv] <- no.cts.fn(LD0.df$species,
                                    log(y/(1 - y)))
      if (vartypes[vv] == "rsel")
      {
        
        # Do Manly's alpha calculations, store.
        no.choices       <- length(avail.list[[vv]])
        choicenames      <- names(avail.list[[vv]])
        avail.vect       <- avail.list[[vv]]
        alpha.mat        <- alpha.fn(LD0.df$species,y,avail.vect)
        alpha.list[[vv]] <- alpha.mat         
        
        # Do niche overlaps, as proportions in categories:
        no.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
      }
    }
    
    #also calculate overall NO measures, averaged over dimensions
    no.overall.mat    <- apply(no.array,c(1,2),mean)
    no.overall.mat.sd <- apply(no.array,c(1,2),sd)
    
    ### Permutation Testing ------------------------------------------------
    
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
    replic <- 1000
    
    pseudo.no.array           <- array(1,c(no.spp,no.spp,no.vars,replic))
    dimnames(pseudo.no.array) <- list(spnames,spnames,varnames,NULL)
    
    # Set a temporary data frame, which will change each time
    # through the cycle by having its species column permuted.
    temp.df <- LD0.df
    
    
    # For each replication, permute the species labels, run the
    # niche overlap calculations, and store the results in the
    # pseudo NO array
    
    ##WARNING: 1000 REPLICATES TAKES ABOUT 11 MINUTES. REDUCE WHEN TESTING
    
    tic()
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
    toc()
    
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
    VV <- ncol(LD0.df[,-c(1:2)])
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
    
    ###save results ------------------------------------------------------------
    NO_0.results <- list(
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
    
##4.3: YEAR 1 GEANGE ANALYSIS ---------------------------------------------------
    ### Set up R objects to store results --------------------------------------
    
    # alpha.list
    
    # The object alpha.list has one component per variable.
    # The components are NULL for ordinary variables.
    
    alpha.list        <- vector("list",no.vars)
    names(alpha.list) <- varnames
    
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
      y <- LD1.df[[varnames[vv]]]
      #this ensures that y is a vector, and not a 1 column tibble 
      #the latter happened with my data and not the example dataset 
      #no idea why but this seems to work 
      
      #adding prints in here for trouble shooting:
      print(paste("vv =", vv))
      print(str(y))
      print(paste("vartype =", vartypes[vv])) 
      if (vartypes[vv] == "bin")
        no.array[,,vv] <- no.bin.fn(LD1.df$species,y)
      if (vartypes[vv] == "cat")
        no.array[,,vv] <- no.cat.fn(LD1.df$species,y)
      if (vartypes[vv] == "count")
        no.array[,,vv] <- no.count.fn(LD1.df$species,y)
      if (vartypes[vv] == "cts")
        no.array[,,vv] <- no.cts.fn(LD1.df$species,y)
      if (vartypes[vv] == "meas")
        no.array[,,vv] <- no.cts.fn(LD1.df$species,log(y))
      if (vartypes[vv] == "pcent")
        no.array[,,vv] <- no.cts.fn(LD1.df$species,
                                    log(y/(100 - y)))
      if (vartypes[vv] == "propn")
        no.array[,,vv] <- no.cts.fn(LD1.df$species,
                                    log(y/(1 - y)))
      if (vartypes[vv] == "rsel")
      {
        
        # Do Manly's alpha calculations, store.
        no.choices       <- length(avail.list[[vv]])
        choicenames      <- names(avail.list[[vv]])
        avail.vect       <- avail.list[[vv]]
        alpha.mat        <- alpha.fn(LD1.df$species,y,avail.vect)
        alpha.list[[vv]] <- alpha.mat         
        
        # Do niche overlaps, as proportions in categories:
        no.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
      }
    }
    
    #also calculate overall NO measures, averaged over dimensions
    no.overall.mat    <- apply(no.array,c(1,2),mean)
    no.overall.mat.sd <- apply(no.array,c(1,2),sd)
    
    ### Permutation Testing ------------------------------------------------
    
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
    replic <- 1000
    
    pseudo.no.array           <- array(1,c(no.spp,no.spp,no.vars,replic))
    dimnames(pseudo.no.array) <- list(spnames,spnames,varnames,NULL)
    
    # Set a temporary data frame, which will change each time
    # through the cycle by having its species column permuted.
    temp.df <- LD1.df
    
    
    # For each replication, permute the species labels, run the
    # niche overlap calculations, and store the results in the
    # pseudo NO array
    
    ##WARNING: 1000 REPLICATES TAKES ABOUT 11 MINUTES. REDUCE WHEN TESTING
    
    tic()
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
    toc()
    
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
    VV <- ncol(LD1.df[,-c(1:2)])
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
    
    ###save results ------------------------------------------------------------
    NO_1.results <- list(
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
    
      #export everything into 1 spreadsheet
    
#SECTION 5: NICHE OVERLAP COMPARISON (YEAR 0/1) --------------------------------
  
  #comparing all the data
  NO_0 <- NO_0.results$ests.overall[upper.tri(NO_0.results$ests.overall)] 
  NO_1 <- NO_1.results$ests.overall[upper.tri(NO_1.results$ests.overall)]
  
  NO_Long.df <- data.frame(
    value = c(NO_0, NO_1),
    group = factor(c(rep("Year 0", length(NO_0)), rep("Year 1", length(NO_1))))
  )
    
  NO_Comp_boxplot <- ggplot(NO_Long.df, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(title = "Niche Overlap by Year",
         x = "Year",
         y = "Niche Overlap Proportion") +
    theme_minimal() 
   
  tNO_Comp <- t.test(NO_0, NO_1)
    
  #not the result we expected, let's try it individually by niche axis 
  grabmat <- function(layer, samp) {
    temp <- samp$NOestimates[, , layer]
    return(temp[upper.tri(temp)])
  }
  
  NO_0.Jul <- grabmat(1, NO_0.results)
  NO_0.Pool <- grabmat(2, NO_0.results)
  NO_0.Rec <- grabmat(3, NO_0.results)
  NO_0.Depth <- grabmat(4, NO_0.results)
  
  NO_1.Jul <- grabmat(1, NO_1.results)
  NO_1.Pool <- grabmat(2, NO_1.results)
  NO_1.Rec <- grabmat(3, NO_1.results)
  NO_1.Depth <- grabmat(4, NO_1.results)
  
  NO_Long_Sep.df <- data.frame( 
    value = c(NO_0.Jul, NO_0.Pool, NO_0.Rec, NO_0.Depth, 
              NO_1.Jul, NO_1.Pool, NO_1.Rec, NO_1.Depth),
    
    measurement_type = factor(rep(c("Julian_Day", "Pool_Number", 
                                    "Recruitment", "Depth"), each = 36, times = 2)),
    
    year = factor(rep(c(0,1), each = 144))
    )
  
  NO_Sep_Comp_boxplot <- ggplot(NO_Long_Sep.df, aes(x = measurement_type, 
                                                    y = value, fill = year)) +
    geom_boxplot() +
    labs(title = "Niche Overlaps by Niche Axis",
         x = "Niche Axis",
         y = "Niche Overlap") +
    theme_minimal() 
  
  ##Bootstrap it 
  
  set.seed(456)
  n_boot <- 36 
  boot_diff <- numeric(n_boot)
  
  for(i in 1:n_boot) {
    year0 <- sample(NO_0, replace = TRUE) 
    year1 <- sample(NO_1, replace = TRUE) 
    boot_diff[i] <- mean(year0) - mean(year1)
  }
  
  hist(boot_diff, main = "Bootstrap dist of mean diff", col = "lightgreen")

  quantile(boot_diff, c(0.025, 0.975))
  #0 is in distribution -> difference is NOT significant  
  
  ex.df <- data.frame( 
    year_0 = year0, 
    year_1 = year1
    )
  
  boxplot(ex.df$year_0, ex.df$year_1)
  
  #ggridges with dataset separate (y0 vs y1) 
  #export all data in into one excel sheet
  #geange but remove winter (run march - september for both)
  
#SECTION 6: GGRIDGES WITH DATASETS SEPARATE ----------------------------------
  Y0_Ridges <- LD0.df %>%
    ggplot(aes(Julian_Day, reorder(species, desc(species)), fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = .5, scale = 8, show.legend = FALSE, color = FALSE) +
    #scale_fill_manual(values=c('#117733','#7EAA44','#88CCEE','#E4BF04','#CC6677','#882255'))+
    scale_x_continuous(limits = c(0, 365), breaks = seq(0, 400, by = 50), expand = c(0, 0)) +
    #coord_fixed(ratio = 10) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    theme(text = element_text(size = 10),
          #axis.text.y = element_text(size = 10, hjust = 2),
          axis.title = element_blank(),
          axis.line.x = element_line(linetype = "solid"))  
  
  Y1_Ridges <- LD1.df %>%
    ggplot(aes(Julian_Day, reorder(species, desc(species)), fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = .5, scale = 8, show.legend = FALSE, color = FALSE) +
    #scale_fill_manual(values=c('#117733','#7EAA44','#88CCEE','#E4BF04','#CC6677','#882255'))+
    scale_x_continuous(limits = c(366, 730), breaks = seq(350, 700, by = 50), expand = c(0, 0)) +
    #coord_fixed(ratio = 10) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    theme(text = element_text(size = 10),
          #axis.text.y = element_text(size = 10, hjust = 2),
          axis.title = element_blank(),
          axis.line.x = element_line(linetype = "solid"))  
  
#SAVE AND COMPILE ALL DATA -----------------------------------------------------
  View(NO_1.results)
  NO_1.results
  variables <- c("Julian_Day", "Pool_Number", )
  
  exportGeange <- function(results, filename) {
    
    #takes the results list from geange analysis and exports it as an xlsx file
    #results = the name of the results list 
    #filename = what you want the xlsx file to be called 
      #couldn't just use the same name for both because of the dot in the results
    
    for (i in 1:dim(results$NOestimates)[3]) {
      df         <- as.data.frame(results$NOestimates[, , i])
      sheet_name <- paste0("NO_Variable_", results[["info"]][["variables"]][i])
      
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet = sheet_name, df, rowNames = TRUE)
    }
    
    for (i in 1:dim(results$separate.pvalues)[3]) {
      df         <- as.data.frame(results$separate.pvalues[, , i])
      sheet_name <- paste0("pNO_Variable_", results[["info"]][["variables"]][i])
      
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet = sheet_name, df, rowNames = TRUE)
    } 
    
    addWorksheet(wb, "Ests_Overall")
    writeData(wb, sheet = "Ests_Overall", as.data.frame(results$ests.overall), 
              rowNames = TRUE)
    
    addWorksheet(wb, "Ests_Overall_sd")
    writeData(wb, sheet = "Ests_Overall_sd", as.data.frame(results$ests.overall.sd), 
              rowNames = TRUE)
    
    addWorksheet(wb, "overall_pvalues")
    writeData(wb, sheet = "overall_pvalues", as.data.frame(results$overall.pvalues), 
              rowNames = TRUE)
    
    sep.clus.p <-  as.data.frame(results[["separate.cluster.pvalues"]])
    sep.diff.p <-  as.data.frame(results[["separate.differentiated.pvalues"]])
    sep.clus.p <- sep.clus.p$`results[["separate.cluster.pvalues"]]`
    sep.diff.p <- sep.diff.p$`results[["separate.differentiated.pvalues"]]`
    
    sep.p <- data.frame(
      separate.cluster.pvalues = sep.clus.p,
      separate.differentiated.pvalues = sep.diff.p
    )
    
    addWorksheet(wb, "separate_values")
    writeData(wb, sheet = "separate_values", sep.p, 
              rowNames = TRUE)
    
    overall.p <- data.frame(
      overall.cluster.pvalues = results[["overall.cluster.pvalues"]],
      overall.differentiated.pvalues = results[["overall.differentiated.pvalues"]]
    )
    
    addWorksheet(wb, "overall_clusdiff_pvalues")
    writeData(wb, sheet = "overall_clusdiff_pvalues", overall.p, 
              rowNames = TRUE) 
    
    name = paste0(filename, ".xlsx")
    saveWorkbook(wb, file = name, overwrite = TRUE)
    
  }
  
  wb <- createWorkbook()
  exportGeange(NO_1.results, "Year_1")
  wb <- createWorkbook()
  exportGeange(NO_0.results, "Year_0")
  wb <- createWorkbook()
  exportGeange(NO_all_4axes.results, "Overall_4axes")
  wb <- createWorkbook() 
  exportGeange(NO_all_6axes.results, "Overall_6axes")
  
  
  
  