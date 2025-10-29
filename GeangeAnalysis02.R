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
    install.packages("lme4")
    install.packages("lmerTest") 
    install.packages("emmeans") 
    install.packages("pbkrtest")
  
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
    library(lmer)
    library(lmerTest) 
    library(emmeans) 
    library(pbkrtest)
  
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
      LD0.df <- LD.df[LD.df$Julian_Day < 250, ] 
      LD1.df <- LD.df[LD.df$Julian_Day > 250, ]
    
      
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
 
  ##5.1: OVERALL COMPARISON ---------------------------------------------------- 
  #extract the upper triangle of each ests.overall matrix for year 0 and 1
  NO_0 <- NO_0.results$ests.overall[upper.tri(NO_0.results$ests.overall)] 
  NO_1 <- NO_1.results$ests.overall[upper.tri(NO_1.results$ests.overall)]
  
  #combine the two into one dataframe that can be used for the ggplot boxplots
  NO_Long.df <- data.frame(
    
    value = c(NO_0, NO_1),
    group = factor(c(rep("Year 0", length(NO_0)), 
                     rep("Year 1", length(NO_1))))
    
  )
  
  #create boxplot comparing raw NO data from 0 and 1  
  NO_Comp_boxplot <- ggplot(NO_Long.df, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(title = "Niche Overlap by Year",
         x     = "Year",
         y     = "Niche Overlap Proportion") +
    theme_minimal() 
  
  #ttest comparing year 0 and 1 
  tNO_Comp <- t.test(NO_0, NO_1)
    
  #not the result we expected, let's try it individually by niche axis 
  
  ##5.2: BY NICHE AXIS COMPARISON ----------------------------------------------
  
  #create a function that will automatically grab the matrix we want from 
  #the NOestimates table 
  grabtri <- function(layer, samp) {
    
    temp <- samp$NOestimates[, , layer]
    return(temp[upper.tri(temp)])
    
  }
  
  #Create vectors for the upper triangle of the matrix for each niche axis
  #and each year
  NO_0.Jul   <- grabtri(1, NO_0.results)
  NO_0.Pool  <- grabtri(2, NO_0.results)
  NO_0.Rec   <- grabtri(3, NO_0.results)
  NO_0.Depth <- grabtri(4, NO_0.results)
  
  NO_1.Jul   <- grabtri(1, NO_1.results)
  NO_1.Pool  <- grabtri(2, NO_1.results)
  NO_1.Rec   <- grabtri(3, NO_1.results)
  NO_1.Depth <- grabtri(4, NO_1.results)
  
  #combine all of that into one dataframe that can be used to make a boxplot
  NO_Long_Sep.df <- data.frame( 
    
    value = c(NO_0.Jul, NO_0.Pool, NO_0.Rec, NO_0.Depth, 
              NO_1.Jul, NO_1.Pool, NO_1.Rec, NO_1.Depth),
    
    measurement_type = factor(rep(c("Julian_Day", "Pool_Number", 
                                    "Recruitment", "Depth"), each = 36, times = 2)),
    
    year = factor(rep(c(0,1), each = 144))
    
    )
  
  #make the boxplot
  NO_Sep_Comp_boxplot <- ggplot(NO_Long_Sep.df, aes(x = measurement_type, 
                                                    y = value, fill = year)) +
    geom_boxplot() +
    labs(title = "Niche Overlaps by Niche Axis",
         x     = "Niche Axis",
         y     = "Niche Overlap") +
    theme_minimal() 
  
  ##Bootstrap time!! 
  
  #bootstrap resampling the year 0 and 1 data
  set.seed(456)
  n_boot    <- 36 
  boot_diff <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    year0        <- sample(NO_0, replace = TRUE) 
    year1        <- sample(NO_1, replace = TRUE) 
    boot_diff[i] <- mean(year0) - mean(year1)
  }
  
  #I don't really remember what I was doing here
  hist(c(year0, year1), main = "Bootstrap dist of mean diff", col = "lightgreen")
  quantile(boot_diff, c(0.025, 0.975))
  #0 is in distribution -> difference is NOT significant  
  
  #quick and dirty boxplot in baseR comparing bootstrapped dataset
  NO_Sep_boot.df <- data.frame( 
    value = c(year0, year1), 
    year = factor(rep(c("Year 0", "Year 1"), each = 36))
    )
  
  NO_Sep_boot.boxplot <- ggplot(NO_Sep_boot.df, aes(x = year, 
                                                    y = value, fill = year)) +
    geom_boxplot() +
    labs(title = "Niche Overlap by Year",
         x     = "Year",
         y     = "Niche Overlap Proportion") +
    theme_minimal() 
  t.test(year0, year1)
  
  
  ##5.3: JULIAN DAY COMPARISON -------------------------------------------------
  
  #create a function that will grab the entire matrix
  grabmat <- function(layer, samp) {
    
    temp <- samp$NOestimates[, , layer]
    return(temp)
    
  }
  
  #grab the julain day matrix for each year
  Jul0 <- grabmat(1, NO_0.results)
  Jul1 <- grabmat(1, NO_1.results)
  
  #make them a real matrix so we can use diag() function
  Jul0 <- as.matrix(Jul0)
  Jul1 <- as.matrix(Jul1)
  
  #set the diagonals to NA so we can remove them 
  diag(Jul0) <- NA 
  diag(Jul1) <- NA
  
  #convert overlap matrices to a dataframe so we can work for it
  Jul0 <- as.data.frame(Jul0)
  Jul1 <- as.data.frame(Jul1)
  
  #convert all the species columns into separate vectors and remove NAs 
  CladJul0 <- Jul0$Cladoceran[!is.na(Jul0$Cladoceran)]
  CopeJul0 <- Jul0$Copepod[!is.na(Jul0$Copepod)] 
  DiptJul0 <- Jul0$Diptera_Pupae[!is.na(Jul0$Diptera_Pupae)]
  MidgJul0 <- Jul0$Midge_Larvae[!is.na(Jul0$Midge_Larvae)]
  MiteJul0 <- Jul0$Mites[!is.na(Jul0$Mites)]
  MosqJul0 <- Jul0$Mosquito_Larvae[!is.na(Jul0$Mosquito_Larvae)]
  OstrJul0 <- Jul0$Ostracod[!is.na(Jul0$Ostracod)]
  RounJul0 <- Jul0$Roundworm[!is.na(Jul0$Roundworm)]
  SpriJul0 <- Jul0$Springtail[!is.na(Jul0$Springtail)]
  
  CladJul1 <- Jul1$Cladoceran[!is.na(Jul1$Cladoceran)]
  CopeJul1 <- Jul1$Copepod[!is.na(Jul1$Copepod)] 
  DiptJul1 <- Jul1$Diptera_Pupae[!is.na(Jul1$Diptera_Pupae)]
  MidgJul1 <- Jul1$Midge_Larvae[!is.na(Jul1$Midge_Larvae)]
  MiteJul1 <- Jul1$Mites[!is.na(Jul1$Mites)]
  MosqJul1 <- Jul1$Mosquito_Larvae[!is.na(Jul1$Mosquito_Larvae)]
  OstrJul1 <- Jul1$Ostracod[!is.na(Jul1$Ostracod)]
  RounJul1 <- Jul1$Roundworm[!is.na(Jul1$Roundworm)]
  SpriJul1 <- Jul1$Springtail[!is.na(Jul1$Springtail)]

  #convert into one big df for boxplot
  NO_Jul_Sep.df <- data.frame(
    NO      = c(CladJul0, CopeJul0, DiptJul0, MidgJul0, MiteJul0, MosqJul0, OstrJul0, RounJul0, SpriJul0, 
                CladJul1, CopeJul1, DiptJul1, MidgJul1, MiteJul1, MosqJul1, OstrJul1, RounJul1, SpriJul1),
    Species = factor(rep(c("Cladoceran", "Copepod", "Diptera_Pupae", 
                           "Midge_Larvae", "Mites", "Mosquito_Larvae", 
                           "Ostracod", "Roundworm", "Springtail"), 
                         each = 8, times = 2)),
    Year    = factor(rep(c(0,1), each = 72))
  )
  
  #make a boxplot for julian day niche overlap 
  NO_Jul_Sep.boxplot <- ggplot(NO_Jul_Sep.df, 
                               aes(x = Species, y = NO, fill = Year)) + 
    geom_boxplot() +
    labs(title = "Niche Overlap by Species (Julian Day Only)",
         x     = "Species",
         y     = "Niche Overlap") +
    theme_minimal() 
  
  #ok now let's just do the overall Julian day comparison 
  
  #use grabtri() to get just the upper triangle for both years
  Jul0tri <- grabtri(1, NO_0.results)
  Jul1tri <- grabtri(1, NO_1.results)
  
  #make into a dataframe we can use for boxplot purposes
  NO_Jul.df <- data.frame(
    values = c(Jul0tri, Jul1tri), 
    year = factor(rep(c("Year 0", "Year 1"), , each = 36))
    
  )
  NO_Jul.df
  #make said boxplot
  NO_Jul_Only.boxplot <- ggplot(NO_Jul.df, aes(x = year, y = values, fill = year)) +
    geom_boxplot() + 
    labs(title = "Niche Overlap Year 0 vs Year 1, Julian Day", 
         x     = "Year", 
         y     = "Niche Overlap") + 
    theme_minimal()
  
  t.test(Jul0tri, Jul1tri) #not significant 
  
  #bootstrap it
  set.seed(789)
  n_boot = 36
  
  for (i in n_boot) { 
    
    Jul0tri.boot <- sample(Jul0tri, replace = TRUE) 
    Jul1tri.boot <- sample(Jul0tri, replace = TRUE)
    
    }
  
  #make into boxplotable dataframe
  NO_Jul.boot.df <- data.frame(
    values = c(Jul0tri.boot, Jul1tri.boot), 
    year   = factor(rep(c("Year 0", "Year 1"), , each = 36))
    
  )
  
  #boxplot it
  NO_Jul_Only.boot.boxplot <- ggplot(NO_Jul.boot.df, aes(x = year, y = values, fill = year)) +
    geom_boxplot() + 
    labs(title = "Niche Overlap Year 0 vs Year 1, Julian Day (Bootstrapped)", 
         x     = "Year", 
         y     = "Niche Overlap") + 
    theme_minimal()
  
  t.test(Jul0tri.boot, Jul1tri.boot) # not significant
  
#SECTION 6: GGRIDGES WITH DATASETS SEPARATE ----------------------------------
  #this section recreates the ggridges graph for each year separately
  
  #year 0
  Jul0_Ridges <- LD0.df |> 
    ggplot(aes(Julian_Day, reorder(species, desc(species)), fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = .5, scale = 8, 
                        show.legend = FALSE, color = FALSE) +
    #scale_fill_manual(values=c('#117733','#7EAA44','#88CCEE','#E4BF04','#CC6677','#882255'))+
    scale_x_continuous(limits = c(0, 365), 
                       breaks = seq(0, 400, by = 50), 
                       expand = c(0, 0)) +
    #coord_fixed(ratio = 10) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    theme(text = element_text(size = 10),
          #axis.text.y = element_text(size = 10, hjust = 2),
          axis.title  = element_blank(),
          axis.line.x = element_line(linetype = "solid"))  
  
  #year 1
  Jul1_Ridges <- LD1.df |> 
    ggplot(aes(Julian_Day, reorder(species, desc(species)), fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = .5, scale = 8, 
                        show.legend = FALSE, color = FALSE) +
    #scale_fill_manual(values=c('#117733','#7EAA44','#88CCEE','#E4BF04','#CC6677','#882255'))+
    scale_x_continuous(limits = c(366, 730), 
                       breaks = seq(350, 700, by = 50), 
                       expand = c(0, 0)) +
    #coord_fixed(ratio = 10) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    theme(text = element_text(size = 10),
          #axis.text.y = element_text(size = 10, hjust = 2),
          axis.title = element_blank(),
          axis.line.x = element_line(linetype = "solid")) 
  
  JulAll_Ridges <- LD.df |> 
    ggplot(aes(Julian_Day, reorder(species, desc(species)), fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = .5, scale = 8, 
                        show.legend = FALSE, color = FALSE) +
    #scale_fill_manual(values=c('#117733','#7EAA44','#88CCEE','#E4BF04','#CC6677','#882255'))+
    scale_x_continuous(limits = c(0, 730), 
                       breaks = seq(0, 700, by = 50), 
                       expand = c(0, 0)) +
    #coord_fixed(ratio = 10) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    theme(text = element_text(size = 10),
          #axis.text.y = element_text(size = 10, hjust = 2),
          axis.title = element_blank(),
          axis.line.x = element_line(linetype = "solid"))  
  
  DepAll_Ridges <- LD.df |> 
    ggplot(aes(Depth_cm, reorder(species, desc(species)), fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = .5, scale = 8, 
                        show.legend = FALSE, color = FALSE) +
    #scale_fill_manual(values=c('#117733','#7EAA44','#88CCEE','#E4BF04','#CC6677','#882255'))+
    scale_x_continuous(limits = c(0, 50), 
                       breaks = seq(0, 50, by = 10), 
                       expand = c(0, 0)) +
    #coord_fixed(ratio = 10) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    theme(text = element_text(size = 10),
          #axis.text.y = element_text(size = 10, hjust = 2),
          axis.title = element_blank(),
          axis.line.x = element_line(linetype = "solid"))  
  
  RecAll_Ridges <- LD.df |> 
    ggplot(aes(Recruitment_Amount, reorder(species, desc(species)), fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, alpha = .5, scale = 8, 
                        show.legend = FALSE, color = FALSE) +
    #scale_fill_manual(values=c('#117733','#7EAA44','#88CCEE','#E4BF04','#CC6677','#882255'))+
    scale_x_continuous(limits = c(0, 600), 
                       breaks = seq(0, 600, by = 50), 
                       expand = c(0, 0)) +
    #coord_fixed(ratio = 10) +
    theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
    theme(text = element_text(size = 10),
          #axis.text.y = element_text(size = 10, hjust = 2),
          axis.title = element_blank(),
          axis.line.x = element_line(linetype = "solid"))  
  
  WD.df <- WideData |> 
    select(Pool_Number, Cladoceran, Copepod, Diptera_Pupae, Midge_Larvae, Mites,
           Mosquito_Larvae, Ostracod, Roundworm, Springtail)    |> 
    mutate(Cladoceran      = na_if(Cladoceran, "NA"))           |>
    mutate(Copepod         = na_if(Copepod, "NA"))              |>
    mutate(Diptera_Pupae   = na_if(Diptera_Pupae, "NA"))        |>
    mutate(Midge_Larvae    = na_if(Midge_Larvae, "NA"))         |>
    mutate(Mites           = na_if(Mites, "NA"))                |>
    mutate(Mosquito_Larvae = na_if(Mosquito_Larvae, "NA"))      |> 
    mutate(Ostracod        = na_if(Ostracod, "NA"))             |>
    mutate(Roundworm       = na_if(Roundworm, "NA"))            |>
    mutate(Springtail      = na_if(Springtail, "NA"))           |> 
    filter(!(Pool_Number %in% "Pond"))                          |>  
    filter(!(Pool_Number %in% "Creek"))                         |> 
    mutate(across(everything(), ~ as.numeric(as.character(.)))) |> 
    group_by(Pool_Number)                                       |> 
    summarise(across(everything(), sum, na.rm = TRUE))  
  
  WD.df <- WD.df %>%
    pivot_longer(
      cols = -Pool_Number,             
      names_to = "Species",
      values_to = "Abundance"
    )  
  
  PoolAll_Bars <-  ggplot(WD.df, aes(x = Species, y = Abundance, fill = Species)) +
    geom_col(position = "identity") +
    facet_wrap(~ Pool_Number, ncol = 1, strip.position = "right") +  # vertical stacking
    theme_minimal() +
    theme(
      strip.text = element_text(angle = 0),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = "Species Abundance Across Pools",
         y = "Abundance",
         x = "Species")
  
  
  
  
  
#SECTION 7: NICHE OVERLAP BY SPECIES -------------------------------------------
  #this section examines the niche overlap values for each species individually
  
  
  ##7.1: NICHE OVERLAP BETWEEN YEARS COMPARISON BY SPECIES ---------------------
  #covert overlap estimates to a matrix
  NO0_ests <- as.matrix(NO_0.results[["ests.overall"]])
  NO1_ests <- as.matrix(NO_1.results[["ests.overall"]])
  
  #set the diagonals to NA so we can remove them 
  diag(NO0_ests) <- NA 
  diag(NO1_ests) <- NA
  
  #convert overlap matrices to a dataframe so we can work for it
  NO0_ests <- as.data.frame(NO0_ests)
  NO1_ests <- as.data.frame(NO1_ests)
  
  #convert all the species columns into separate vectors and remove NAs 
  Clad0 <- NO0_ests$Cladoceran[!is.na(NO0_ests$Cladoceran)]
  Cope0 <- NO0_ests$Copepod[!is.na(NO0_ests$Copepod)] 
  Dipt0 <- NO0_ests$Diptera_Pupae[!is.na(NO0_ests$Diptera_Pupae)]
  Midg0 <- NO0_ests$Midge_Larvae[!is.na(NO0_ests$Midge_Larvae)]
  Mite0 <- NO0_ests$Mites[!is.na(NO0_ests$Mites)]
  Mosq0 <- NO0_ests$Mosquito_Larvae[!is.na(NO0_ests$Mosquito_Larvae)]
  Ostr0 <- NO0_ests$Ostracod[!is.na(NO0_ests$Ostracod)]
  Roun0 <- NO0_ests$Roundworm[!is.na(NO0_ests$Roundworm)]
  Spri0 <- NO0_ests$Springtail[!is.na(NO0_ests$Springtail)]
  
  Clad1 <- NO1_ests$Cladoceran[!is.na(NO1_ests$Cladoceran)] 
  Cope1 <- NO1_ests$Copepod[!is.na(NO1_ests$Copepod)]
  Dipt1 <- NO1_ests$Diptera_Pupae[!is.na(NO1_ests$Diptera_Pupae)]
  Midg1 <- NO1_ests$Midge_Larvae[!is.na(NO1_ests$Midge_Larvae)]
  Mite1 <- NO1_ests$Mites[!is.na(NO1_ests$Mites)]
  Mosq1 <- NO1_ests$Mosquito_Larvae[!is.na(NO1_ests$Mosquito_Larvae)]
  Ostr1 <- NO1_ests$Ostracod[!is.na(NO1_ests$Ostracod)]
  Roun1 <- NO1_ests$Roundworm[!is.na(NO1_ests$Roundworm)]
  Spri1 <- NO1_ests$Springtail[!is.na(NO1_ests$Springtail)]
  
  #turn all the above vectors into a dataframe for a boxplot
  Ests_By_Sp.df <- data.frame(
    
    NO      = c(Clad0, Cope0, Dipt0, Midg0, Mite0, Mosq0, Ostr0, Roun0, Spri0, 
                Clad1, Cope1, Dipt1, Midg1, Mite1, Mosq1, Ostr1, Roun1, Spri1),
    Species = factor(rep(c("Cladoceran", "Copepod", "Diptera_Pupae", 
                           "Midge_Larvae", "Mites", "Mosquito_Larvae", 
                           "Ostracod", "Roundworm", "Springtail"), 
                         each = 8, times = 2)),
    Year    = factor(rep(c(0,1), each = 72, times = 2))
    
    )
  
  #create a box plot for niche overlap estimates by species by year
  Ests_By_Sp.boxplot <- ggplot(Ests_By_Sp.df, 
                               aes(x = Species, y = NO, fill = Year)) + 
    geom_boxplot() +
    labs(title = "Niche Overlap by Species",
         x     = "Species",
         y     = "Niche Overlap") +
    theme_minimal() 
  
  
  #bootstrap the values for niche overlap estimates by species by year
  n_boot <- 8
  for (i in 1:n_boot) { 
    
    Clad0_boot <- sample(Clad0, replace = TRUE)
    Cope0_boot <- sample(Cope0, replace = TRUE)
    Dipt0_boot <- sample(Dipt0, replace = TRUE)
    Midg0_boot <- sample(Midg0, replace = TRUE)
    Mite0_boot <- sample(Mite0, replace = TRUE)
    Mosq0_boot <- sample(Mosq0, replace = TRUE)
    Ostr0_boot <- sample(Ostr0, replace = TRUE)
    Roun0_boot <- sample(Roun0, replace = TRUE)
    Spri0_boot <- sample(Spri0, replace = TRUE)
    
    Clad1_boot <- sample(Clad1, replace = TRUE)
    Cope1_boot <- sample(Cope1, replace = TRUE)
    Dipt1_boot <- sample(Dipt1, replace = TRUE)
    Midg1_boot <- sample(Midg1, replace = TRUE)
    Mite1_boot <- sample(Mite1, replace = TRUE)
    Mosq1_boot <- sample(Mosq1, replace = TRUE)
    Ostr1_boot <- sample(Ostr1, replace = TRUE)
    Roun1_boot <- sample(Roun1, replace = TRUE)
    Spri1_boot <- sample(Spri1, replace = TRUE)
    
  }
  
  #recreate the same boxplot from above but with the bootstrapped data
  #new dataframe
  Ests_By_Sp.df.boot <- data.frame( 
    
    NO      = c(Clad0_boot, Cope0_boot, Dipt0_boot, Midg0_boot, Mite0_boot, 
                Mosq0_boot, Ostr0_boot, Roun0_boot, Spri0_boot, 
                Clad1_boot, Cope1_boot, Dipt1_boot, Midg1_boot, Mite1_boot, 
                Mosq1_boot, Ostr1_boot, Roun1_boot, Spri1_boot),
    Species = factor(rep(c("Cladoceran", "Copepod", "Diptera_Pupae", 
                           "Midge_Larvae", "Mites", "Mosquito_Larvae", 
                           "Ostracod", "Roundworm", "Springtail"), 
                          each = 8, times = 2)),
    Year    = factor(rep(c(0,1), each = 72, times = 2))
    
  )
  
  #bootstrapped boxplot
  Ests_By_Sp.boxplot.boot <- ggplot(Ests_By_Sp.df, 
                                    aes(x = Species, y = NO, fill = Year)) + 
    geom_boxplot() +
    labs(title = "Niche Overlap by Species by Year",
         x     = "Species",
         y     = "Niche Overlap") +
    theme_minimal() 
  
  #ttest every pair of niche overlaps for each species between years
  Cladt <- t.test(Clad0_boot, Clad1_boot)
  Copet <- t.test(Cope0_boot, Cope1_boot)
  Diptt <- t.test(Dipt0_boot, Dipt1_boot)
  Midgt <- t.test(Midg0_boot, Midg1_boot)
  Mitet <- t.test(Mite0_boot, Mite1_boot)
  Mosqt <- t.test(Mosq0_boot, Mosq1_boot)
  Ostrt <- t.test(Ostr0_boot, Ostr1_boot)
  Rount <- t.test(Roun0_boot, Roun1_boot)
  Sprit <- t.test(Spri0_boot, Spri1_boot) #significant
  
  #7.2: COMPARING EACH SPECIES OVERALL NO TO EACH OTHER (NOT SEPARATED BY YEAR)
  
  #run a t-test on every pair of species comparing their niche overlap values
  
    #combine data from bootstraps of both years into one vector for each species
      Clad_boot <- c(Clad0_boot, Clad1_boot)
      Cope_boot <- c(Cope0_boot, Cope1_boot)
      Dipt_boot <- c(Dipt0_boot, Dipt1_boot)
      Midg_boot <- c(Midg0_boot, Midg1_boot)
      Mite_boot <- c(Mite0_boot, Mite1_boot)
      Mosq_boot <- c(Mosq0_boot, Mosq1_boot)
      Ostr_boot <- c(Ostr0_boot, Ostr1_boot)
      Roun_boot <- c(Roun0_boot, Roun1_boot)
      Spri_boot <- c(Spri0_boot, Spri1_boot)
    
    #turn those vectors into a list
      NO_boot_list <- list(Clad_boot, Cope_boot, Dipt_boot, Midg_boot, Mite_boot, 
                           Mosq_boot, Ostr_boot, Roun_boot, Spri_boot)
    
    #give each item in the list the correct name
      names(NO_boot_list) <- c("Cladoceran", "Copepod", "Diptera_Pupae", 
                               "Midge_Larvae", "Mites", "Mosquito_Larvae", 
                               "Ostracod", "Roundworm", "Springtail")
   
   #create an empty 9x9 matrix to fill in with data from our t tests
     t_NOmatrix <- matrix(NA, nrow = 9, ncol = 9) 
     #give the rows and columns correct names
     rownames(t_NOmatrix) <- colnames(t_NOmatrix) <- names(NO_boot_list)
    
    #for loop to run a t-test for every pair of species (except for diagonals)
     for (i in 1:9) {
       for (j in 1:9) {
         if (i == j) {
           t_NOmatrix[i, j] <- NA
         } else {
           test             <- t.test(NO_boot_list[[i]], NO_boot_list[[j]])
           t_NOmatrix[i, j] <- test$p.value
         }
       }
      }
   
  
   #round above matrix to a reasonable number of digits 
   t_NOmatrix <- round(t_NOmatrix, 3)
   
   #remove all of the non-significant items
   only_sigs_tNOmat <-  t_NOmatrix
   only_sigs_tNOmat[only_sigs_tNOmat >= 0.05] <- NA 
   
   #create a dataframe for the overall, not year-separated data
   Overall_By_Sp.df.boot <- data.frame(
     NO = c(Clad_boot, Cope_boot, Dipt_boot, Midg_boot, Mite_boot, 
             Mosq_boot, Ostr_boot, Roun_boot, Spri_boot),
     Species = factor(rep(c("Cladoceran", "Copepod", "Diptera_Pupae", 
                            "Midge_Larvae", "Mites", "Mosquito_Larvae", 
                            "Ostracod", "Roundworm", "Springtail"), 
                          each = 16))
   )
    
   #make a boxplot for the overall, not year-separated data
   Overall_By_Sp.boxplot.boot <- ggplot(Overall_By_Sp.df.boot, 
                                     aes(x = Species, y = NO, fill = Species)) + 
     geom_boxplot() +
     labs(title = "Overall Niche Overlap by Species",
          x     = "Species",
          y     = "Niche Overlap") +
     theme_minimal() 

 ##7.2: COMPARISON TO OVERALL FOR EACH SPECIES ---------------------------------
 #grab the upper tri for the overall data
   NO_Overall <- NO_all_4axes.results$ests.overall[
   upper.tri(NO_all_4axes.results$ests.overall)] 
 
   #bootstrap overall data
   n_boot <- 36
   for (i in n_boot) {
    NO_Overall.boot <- sample(NO_Overall, replace = TRUE)
   }
  
  #run a t.test comparing each species to the overall niche overlap 
  #I did this with a for loop because I thought it would be faster and easier. 
  #It wasn't but I'm proud of it so it's staying
  
    #create a vector of the shortenings I made for each species
    spshort        <- c("Clad", "Cope", "Dipt", "Midg", "Mite", "Mosq", "Ostr",
                        "Roun", "Spri")
    tNO_By_Species <- 1:9
    
    #for loop running a ttest for each species
    for (i in 1:9) {
      
      sp                <- spshort[i]
      temp              <- get(paste0(sp, "_boot"))
      test              <- t.test(temp, NO_Overall.boot)
      tNO_By_Species[i] <- test$p.value 
      
    }
    
    #convert to a list so I can name the values
    tNO_By_Species        <- as.list(tNO_By_Species)
    names(tNO_By_Species) <- species05
  
  ##7.3: COMPARING TO OVERALL BY YEAR FOR EACH SPECIES -------------------------
    #bootstrap year 0 and 1 data
    n_boot <- 36
    for (i in n_boot) {
      NO_0.boot <- sample(NO_0, replace = TRUE)
      NO_1.boot <- sample(NO_1, replace = TRUE)
    }
    
    #same for loop from above for years 0 and 1
    tNO_0 <- 1:9
    tNO_1 <- 1:9
    for (i in 1:9) {
      
      sp                <- spshort[i]
      temp              <- get(paste0(sp, "0_boot"))
      test              <- t.test(temp, NO_0.boot)
      tNO_0[i]          <- test$p.value 
      
    }
    for (i in 1:9) {
      
      sp                <- spshort[i]
      temp              <- get(paste0(sp, "1_boot"))
      test              <- t.test(temp, NO_1.boot)
      print(test$p.value)
      tNO_1[i]          <- test$p.value 
    }
    
    tNO_0 #only springtails are significantly different
    tNO_1 #springtails, ostracods, and mosquito larvae significantly different
    
    tNO_0        <- as.list(tNO_0)
    names(tNO_0) <- species05
    
    tNO_1        <- as.list(tNO_1)
    names(tNO_1) <- species05

 
    
    
    
  Jul.linmodel <- lm(values ~ year, data = NO_Jul.boot.df)
  anova(Jul.linmodel)
  
  
  Jul.emmeans <- emmeans(Jul.linmodel, ~ year)
  Jul.emmeans.df <- as.data.frame(Jul.emmeans)
  
#SAVE AND COMPILE ALL DATA -------Jul.emmeans.df#SAVE AND COMPILE ALL DATA -----------------------------------------------------
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
  
  
  
  