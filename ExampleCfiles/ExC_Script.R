# NOc

# A general method for combining different data types into
# a unified multivariate analysis of niche overlap.
# This program calculates niche overlaps over multiple niche axes
# when each axis of niche space is measured using DIFFERENT
# individuals.

# Program reads in individual data,  calculates 
# electivity scores for any resource usage variables,
# calculates niche overlaps between pairs of species, and
# runs null model tests of (i) differential use of niche 
# space by multiple species; and (ii) even distribution 
# of species across niche space.

# This illustration examines mean niche overlap between five species
# of plants based on two functional traits (1) Surface Leaf Area 
# (measurement data); and (2) Elevational distribution (continuous 
# data).

# Data for each axis is read into the program from an external datafile:
# 1. "MEE3_070_sm_Examle_C_sla.txt"
# 2. "MEE3_070_sm_Examle_C_elevation.txt"

# Likewise, functions are read into the program from
# the file "MEE3_070_sm_NicheFunctions.txt"

# The analysis is divided into six distinct sections, 
# each identified by
#############################################################
# 1. Analysis A: Surface Leaf Area
# 2. Analysis B: Elevation
# 3. The combination of the three axes into a single array
#    for oberved data, and for pseudo data, and the calulation
#    of mean niche overlap across axes
# 4. Tests determining if two species occupy different niche space,
#    seperately for each each axis, and for niche overlap averaged
#    across axes.
# 5. Tests determining if the species are evenly distributed across,
#    or clumped within niche space. This is done seperately for each
#    each axis, and for niche overlap averaged across axes.
# 6. The compilation of the above results into a single R object
#    for ease of comparison/printing/saving

# Cut and paste the following commands into R.
# They may be pasted into R in groups of rows - all comments 
# are preceded by a hash, so will not be acted on in R.

# Sections of the program which require input from the user
# are bracketed by
# ??????????????????????????????????????????????????

# Required packages
library(abind) 

# -----------------------------------------------------------
#             Functions Input
# -----------------------------------------------------------


# Read in the niche overlap functions:
setwd("C:/Users/jackb/OneDrive/Desktop/Fall2025PoolsAnalysis/ExampleCfiles")
source("MEE3_070_sm_NicheFunctions.txt")


############################################################# 
# -----------------------------------------------------------
#          Analysis A: SLA
# -----------------------------------------------------------
#############################################################

# Analysis for calculating niche overlap on a single niche axis,
# in this instance a measurement variable.

# The input data set is a .txt file.
# It needs to have its first two columns
# labelled "id" and "species". Subsequent columns are
# individual-based variables, one per column.

# Input the individual data file:

# ??????????????????????????????????????????????????  
A.df <- read.table("MEE3_070_sm_ExampleC_sla.txt",T)
# ??????????????????????????????????????????????????  

# Ensure the first two column names are "id" and "species".
colnames(A.df)[1] <- "id"
colnames(A.df)[2] <- "species"

# Ensure that the first 2 cols are factors.
A.df$id      <- as.factor(A.df$id)
A.df$species <- as.factor(A.df$species)

# Store some vectors of names:
spnames   <- sort(unique(as.character(A.df$species)))
no.spp    <- length(spnames)

varnames <- colnames(A.df)[-(1:2)]    
no.vars  <- length(varnames) 

# Make a vector of variable types to match the variable names:

# "cat"   = categorical, but not resource selection
# "bin"   = binary
# "cts"   = continuous, use raw data (no transformation)
# "meas"  = measurement, continuous positive, take logs
# "pcent" = percentage data, bounds at 0 and 100, use logits
# "propn" = proportion data, bounds at 0 and 1, use logits
# "count" = count
# "rsel"  = resource selection, categorical

# ?????????????????????????????????????????????????? 
vartypes <- c("meas")
# ?????????????????????????????????????????????????? 

# Check they are correctly labelled:
cbind(varnames,vartypes)

# Any variables of resource selection type should have an
# associated availability vector for the site.

# Set up a list of objects which are NULL if this is not
# a resource selection variable, and with the availability
# vector if it is resource selection.
avail.list <- vector("list",no.vars)
names(avail.list) <- varnames

# Set the right matrices into the list.
# Run through the following routine for each resource selection
# type of variable (not done here because this is a measurement variable).

# ??????????????????????????????????????????????????
# ----------------------------------------------------------------
# Routine start:
# Read in the availability vector for this each resource from a file.
# These are percentages of the various choices.

# avail.vect <- read.csv("Filename.csv",T)[1,]

# Sort alphabetically:

# avail.vect <- avail.vect[sort.list(names(avail.vect))]

# Check that the resource types used are all in the availability
# vector, need names to match.

# used <- levels(A.df$habitat)
# used %in% names(avail.vect)

# If not all TRUE, go back to data files, rename to make them match,
# start analysis again.

# Ensure the availabilities are percentages:

# avail.vect <- avail.vect/sum(avail.vect)*100
# If they were already percentages, this makes no change.

# avail.list[[1]] <- avail.vect  # Stored in first component, as
# "habitat" was the first variable.

# Routine end.
# -----------------------------------------------------------------

# Read in more availability data if required for other variables.
# Do the routine above, between the dashed lines, for each "rsel"
# variable.
# Put them in the correct component of the list.
# ??????????????????????????????????????????????????

# Next look at avail.list, check it seems right - availability
# vectors should match "rsel" type variables, NULL elsewhere.
# Different options within each resource should be in alphabetical
# order.
avail.list




# Set up R objects to store results
# ---------------------------------

# alpha.list

# The object alpha.list has one component per variable.
# The components are NULL for ordinary variables.
# For resource selection variables, the component is
# the matrix of Manly's alpha values for that variable. 
# The matrix has:
# Rows = species,
# Cols = choices for that resource

Aalpha.list <- vector("list",no.vars)
names(Aalpha.list) <- varnames

for (vv in 1:no.vars) if (vartypes[vv]=="rsel")
{
  choices <- unique(A.df[,vv+2])
  no.ch   <- length(choices)
  Aalpha.list[[vv]] <- matrix(NA,no.spp,no.ch)
  dimnames(Aalpha.list[[vv]]) <- list(spnames,choices)
}


# no.array

# Set up an array of niche overlaps.
# The object no.array is an array of niche overlaps.
# It is a 3-D array, with rows and columns being species 
# (a square symmetric matrix for pairwise niche overlaps), 
# and the layers are the dimensions for the multivariate 
# niche overlap measure (one dimension per variable).
# Rows and columns are species, layers are variables.

Ano.array  <- array(1,c(no.spp,no.spp,no.vars))
dimnames(Ano.array) <- list(spnames,spnames,varnames)

# Run through each variable in turn, identify its type,
# calculate the appropriate NO matrix and store it in
# the right layer of the no.array.
for (vv in 1:no.vars)
{
  y <- A.df[,colnames(A.df)==varnames[vv]]
  if (vartypes[vv] == "bin")
    Ano.array[,,vv] <- no.bin.fn(A.df$species,y)
  if (vartypes[vv] == "cat")
    Ano.array[,,vv] <- no.cat.fn(A.df$species,y)
  if (vartypes[vv] == "count")
    Ano.array[,,vv] <- no.count.fn(A.df$species,y)
  if (vartypes[vv] == "cts")
    Ano.array[,,vv] <- no.cts.fn(A.df$species,y)
  if (vartypes[vv] == "meas")
    Ano.array[,,vv] <- no.cts.fn(A.df$species,log(y))
  if (vartypes[vv] == "pcent")
    Ano.array[,,vv] <- no.cts.fn(A.df$species,
                                 log(y/(100-y)))
  if (vartypes[vv] == "propn")
    Ano.array[,,vv] <- no.cts.fn(A.df$species,
                                 log(y/(1-y)))
  if (vartypes[vv] == "rsel")
  {
    
    # Do Manly's alpha calculations, store.
    avail.vect <- avail.list[[vv]]
    alpha.mat <- alpha.fn(A.df$species,y,avail.vect)
    alpha.list[[vv]] <- alpha.mat         
    
    # Do niche overlaps, as proportions in categories:
    Ano.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
  }
}


# Analysis A  -  Permutation testing.
# -----------------------------------

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

# ??????????????????????????????????????????????????
# Choose no. of replications.
# Start low, eg. with 10 reps, to check it is working.
# Then do more reps, e.g. 1000 reps for 3 decimal places in p-values.
replic <- 10
# ??????????????????????????????????????????????????

# Set up array to store pseudo niche overlaps:
Apseudo.no.array  <- array(1,c(no.spp,no.spp,no.vars,replic))
dimnames(Apseudo.no.array) <- list(spnames,spnames,varnames,NULL)

# Set a temporary data frame, which will change each time
# through the cycle by having its species column permuted.
temp.df <- A.df

# For each replication, permute the species labels, run the
# niche overlap calculations, and store the results in the
# pseudo NO array.
for (rr in 1:replic)
{
  
  # Permute the species labels in the temporary dataframe:
  temp.df$species <- sample(temp.df$species)
  for (vv in 1:no.vars)
  {
    
    # Read out the column from this variable:
    y <- temp.df[,colnames(A.df)==varnames[vv]]
    
    # Run through the variable types, do appropriate analyses:
    if (vartypes[vv] == "bin")
      Apseudo.no.array[,,vv,rr] <- no.bin.fn(temp.df$species,y)
    if (vartypes[vv] == "cat")
      Apseudo.no.array[,,vv,rr] <- no.cat.fn(temp.df$species,y)
    if (vartypes[vv] == "count")
      Apseudo.no.array[,,vv,rr] <- no.count.fn(temp.df$species,y)
    if (vartypes[vv] == "cts")
      Apseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,y)
    if (vartypes[vv] == "meas")
      Apseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,log(y))
    if (vartypes[vv] == "pcent")
      Apseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(100-y)))
    if (vartypes[vv] == "propn")
      Apseudo.no.array[,,vv,rr] <- no.cts.fn(temp.df$species,
                                             log(y/(1-y)))
    if (vartypes[vv] == "rsel")
    {
      
      # Do Manly's alpha calculations, store.
      avail.vect <- avail.list[[vv]]
      alpha.mat  <- alpha.fn(temp.df$species,y,avail.vect)
      
      # Do niche overlaps, as proportions in categories:
      Apseudo.no.array[,,vv,rr] <- no.rsel.cat.fn(alpha.mat)
    }
  }
  print(paste("Rep",rr,"done"))
}

# Calculate p values for each pair of species 
# separately over all variables.


############################################################# 

# -----------------------------------------------------------
#          Analysis B: elevation
# -----------------------------------------------------------
############################################################# 

# Analysis for calculating niche overlap on a single niche axis,
# in this instance a continuous variable.

# The input data set is a .txt file.
# It needs to have its first two columns
# labelled "id" and "species". Subsequent columns are
# individual-based variables, one per column.

# Input the individual data file:

# ?????????????????????????????????????????????????????????? 
# read in the individual data file
B.df <- read.table("MEE3_070_sm_ExampleC_elevation.txt",T)
# ??????????????????????????????????????????????????????????

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

# Make a vector of variable types to match the variable names:
# "cat"   = categorical, but not resource selection
# "bin"   = binary
# "cts"   = continuous, use raw data (no transformation)
# "meas"  = measurement, continuous positive, take logs
# "pcent" = percentage data, bounds at 0 and 100, use logits
# "propn" = proportion data, bounds at 0 and 1, use logits
# "count" = count
# "rsel"  = resource selection, categorical


# ??????????????????????????????????????????????????????????
vartypes <- c("cts")
# ??????????????????????????????????????????????????????????


# Check they are correctly labelled:
cbind(varnames,vartypes)

# Any variables of resource selection type should have an
# associated availability vector for the site.

# Set up a list of objects which are NULL if this is not
# a resource selection variable, and with the availability
# vector if it is resource selection.
avail.list <- vector("list",no.vars)
names(avail.list) <- varnames

# Set the right matrices into the list.
# Run through the following routine for each resource selection
# type of variable (not done here because this is a continuous variable).

# ??????????????????????????????????????????????????
# ----------------------------------------------------------------
# Routine start:
# Read in the availability vector for this each resource from a file.
# These are percentages of the various choices.

# avail.vect <- read.csv("Filename.csv",T)[1,]

# Sort alphabetically:

# avail.vect <- avail.vect[sort.list(names(avail.vect))]

# Check that the resource types used are all in the availability
# vector, need names to match.

# used <- levels(A.df$habitat)
# used %in% names(avail.vect)

# If not all TRUE, go back to data files, rename to make them match,
# start analysis again.

# Ensure the availabilities are percentages:

# avail.vect <- avail.vect/sum(avail.vect)*100
# If they were already percentages, this makes no change.

# avail.list[[1]] <- avail.vect  # Stored in first component, as
# "habitat" was the first variable.

# Routine end.
# -----------------------------------------------------------------

# Read in more availability data if required for other variables.
# Do the routine above, between the dashed lines, for each "rsel"
# variable.
# Put them in the correct component of the list.
# ??????????????????????????????????????????????????

# Next look at avail.list, check it seems right - availability
# vectors should match "rsel" type variables, NULL elsewhere.
# Different options within each resource should be in alphabetical
# order.
avail.list


# Set up R objects to store results
# ---------------------------------

# alpha.list

# The object alpha.list has one component per variable.
# The components are NULL for ordinary variables.
# For resource selection variables, the component is
# the matrix of Manly's alpha values for that variable. 
# The matrix has:
# Rows = species,
# Cols = choices for that resource (e.g. the resource
# "habitat" may have choices "grass", "rock", "forest"),

Balpha.list <- vector("list",no.vars)
names(Balpha.list) <- varnames

for (vv in 1:no.vars) if (vartypes[vv]=="rsel")
{
  choices <- unique(B.df[,vv+2])
  no.ch   <- length(choices)
  Balpha.list[[vv]] <- matrix(NA,no.spp,no.ch)
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

Bno.array  <- array(1,c(no.spp,no.spp,no.vars))
dimnames(Bno.array) <- list(spnames,spnames,varnames)

# Run through each variable in turn, identify its type,
# calculate the appropriate NO matrix and store it in
# the right layer of the no.array.
for (vv in 1:no.vars)
{
  y <- B.df[,colnames(B.df)==varnames[vv]]
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
    avail.vect <- avail.list[[vv]]
    alpha.mat <- alpha.fn(B.df$species,y,avail.vect)
    alpha.list[[vv]] <- alpha.mat         
    
    # Do niche overlaps, as proportions in categories:
    Bno.array[,,vv] <- no.rsel.cat.fn(alpha.mat)
  }
}


# Analysis B  -  Permutation testing.
# -----------------------------------

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

# Set up array to store pseudo niche overlaps:
Bpseudo.no.array  <- array(1,c(no.spp,no.spp,no.vars,replic))
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
    y <- temp.df[,colnames(B.df)==varnames[vv]]
    
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



#############################################################   
# -----------------------------------------------------------
# Combine the two axes (sla and elevation) into a
# single dataframe and calculate mean niche overlap across axes
# ----------------------------------------------------------- 
############################################################# 

# Combine the two axes (sla and elevation) into a
# single dataframe   
no.all.mat <- abind(Ano.array,Bno.array,along=3) 

# calculate mean niche overlap across axes  
mean_overlap <- apply(no.all.mat,1:2,mean)

# calculate the associated standard deviation
sd_overlap <- apply(no.all.mat,1:2,sd)

# Combine the two pseudo datasets (sla and elevation) 
# into a single array   
pseudo.no.all.mat <- abind(Apseudo.no.array,Bpseudo.no.array,along=3) 

# For each replicate, calculate mean niche overlap across axes 
pseudo.mean_overlap <- apply(pseudo.no.all.mat,c(1:2,4),mean)



#############################################################
#--------------------------------------------------------
# Null model analysis determining if the niches
# of two species in niche space differ
#--------------------------------------------------------
#############################################################

# Calculate p values for each pair of species 
# separately for each variable.

# calculate number of axes
no.axes <- dim(no.all.mat)[3]

#????????????????????????????????????????????????????????
# assign name to each axis
axis.names <- c("sla", "elevation")
#????????????????????????????????????????????????????????

# calculate seperate p-values for each axis
sep.pvals     <- array(1,c(no.spp,no.spp,no.axes))
dimnames(sep.pvals) <- list(spnames,spnames,axis.names)

for (spa in 1:(no.spp-1)) for (spb in (spa+1):no.spp)
  for (vv in 1: no.axes)   
  {
    pseudo.nos <- pseudo.no.all.mat[spa,spb,vv,]
    data.no    <- no.all.mat[spa,spb,vv]
    sep.pvals[spa,spb,vv] <- mean(pseudo.nos<data.no) 
    length(pseudo.nos[data.no<pseudo.nos])
    sep.pvals[spb,spa,vv] <- sep.pvals[spa,spb,vv] 
  }


# Also calculate a p-value for overall NO measure averaged across axes
overall.pvals <- matrix(1,no.spp,no.spp)
dimnames(overall.pvals) <- list(spnames,spnames)

for (spa in 1:(no.spp-1)) for (spb in (spa+1):no.spp)
{
  temp.mat  <- pseudo.no.all.mat[spa,spb,,]
  pseudo.nos <- apply(temp.mat,2,mean)
  data.no    <- mean_overlap[spa,spb]
  overall.pvals[spa,spb] <- mean(pseudo.nos<data.no) 
  length(pseudo.nos[data.no<pseudo.nos])
  overall.pvals[spb,spa] <- overall.pvals[spa,spb] 
}



#############################################################
#--------------------------------------------------------
# Null model analysis determining if the distribution of
# species across niche space are more differentiated
# or more clustered than expected
#--------------------------------------------------------
#############################################################

# First, reformat the observed data to derive a matrix of niche overlaps
# with one row per species combination, and one column for each niche dimension

VV <- no.axes  # Number of axes
RR <- replic   # Number of replications.

no.mat <- matrix(NA,(no.spp*(no.spp-1)/2),VV)
for (vv in 1:VV)
  no.mat[,vv] <- as.vector(as.dist(no.all.mat[,,vv]))

# Next, reformat the pseudo data to derive a matrix of niche overlaps
# with one row per species, and one column for each niche dimension,
# with one extra dimension, one layer for each replication

pseudo.mat <- 	array(NA,c((no.spp*(no.spp-1)/2),VV,replic))
for (vv in 1:VV) for (rr in 1:RR)
  pseudo.mat[,vv,rr] <- as.vector(as.dist(pseudo.no.all.mat[,,vv,rr]))

# --------------------------------------------------------
# For each niche dimension, calculate mean and variance over the species
# pairs, and hence the test statistic ch = coefficient of heterogeneity.
# Note: Need to use variance formula based on n, not n-1.

KK <- ncol(no.mat)      # Number of niche dimensions
SS <- nrow(no.mat)      # Number of species pairs
RR <- replic            # Number of replications.


data.ch <- rep(NA,KK)
pseudo.ch <- matrix(NA,RR,KK)

for (kk in 1:KK)
{
  # Calculate data test statistic:
  x <- mean(no.mat[,kk])
  v <- var(no.mat[,kk])*(SS-1)/SS # Adjust for denom n, not n-1
  data.ch[kk] <- v/x/(1-x)
  
  # Calculate test stats for all pseudo-data:
  for (rr in 1:RR)
  {
    x <- mean(pseudo.mat[,kk,rr])
    v <- var(pseudo.mat[,kk,rr])*(SS-1)/SS
    pseudo.ch[rr,kk] <- v/x/(1-x)
  }
}

# For each niche dimension, see if data more differentiated than random.
p.dims.diff <- rep(NA,KK)
for (kk in 1:KK)
  p.dims.diff[kk] <- mean(data.ch[kk] > pseudo.ch[,kk])
names(p.dims.diff) <- paste("diff.dim",sort(axis.names))

# For each niche dimension, see if data more clustered than random.
p.dims.clus <- rep(NA,KK)
for (kk in 1:KK)
  p.dims.clus[kk] <- mean(data.ch[kk] < pseudo.ch[,kk])
names(p.dims.clus) <- paste("clus.dim",sort(axis.names))


# --------------------------------------------------------
# For average niche overlap, calculate mean and variance over the species
# pairs, and hence the test statistic ch = coefficient of heterogeneity.
# Note: Need to use variance formula based on n, not n-1.

overall.data.ch   <- mean(data.ch)
overall.pseudo.ch <- apply(pseudo.ch,1,mean)

# Test if this community is more differentiated than random:
p.all.diff <- mean(overall.data.ch > overall.pseudo.ch)

# Test if this community is more clustered than random:
p.all.clus <- mean(overall.data.ch < overall.pseudo.ch)


#############################################################
#--------------------------------------------------------
# Save all results of the analysis:
#--------------------------------------------------------
#############################################################

NOc.results <- list(
  info = list(variables = cbind(axis.names,c("meas","ctns")),
              perm.reps = replic),
  NOestimates = no.all.mat,
  separate.pvalues = sep.pvals,
  separate.cluster.pvalues = p.dims.clus,
  separate.differentiated.pvalues = p.dims.diff,
  ests.overall = mean_overlap,
  ests.overall.sd = sd_overlap,
  overall.pvalues = overall.pvals,
  overall.cluster.pvalues = p.all.clus,
  overall.differentiated.pvalues = p.all.diff)

# ??????????????????????????????????????????????????
# To inspect results later, type in
#    names(NOb.results)
# to decide what to look at. Then type (e.g.)
#    NOb.results$NOestimates
# to see that component of the list.

# Save the NOb.results with a more informative names for 
# your own data set.
NOc.results

# NOc.results$ests.overall may be viewed as a measure of association
# between species pairs; therefore, 1 - NOb.results$ests.overall
# may be viwed as a measure of distance, which can be used in non-Metric
# Multidimensional Scaling (nMDS) or Principle Coordinates Analysis (PCoA)
# to display relationships between the species. We recommend using a 
# package such as VEGAN to conduct nMDS or ECODIST to conduct PCoA
# and construct visual representations of species associations.
