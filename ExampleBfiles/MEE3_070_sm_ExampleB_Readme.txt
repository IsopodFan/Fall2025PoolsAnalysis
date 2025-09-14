# NOb.R

# Program to calculate niche overlaps over multiple niche
# axes per individual.

# Program reads in individual data,  calculates electivity scores for any resource usage variables,
# calculates niche overlaps between pairs of species, and runs null model tests of (i) differential use of niche 
# space by multiple species; and (ii) even distribution of species across niche space.

# Cut and paste the commands into R. They may be pasted into R in groups of rows - all comments 
# are preceded by a hash, so will not be acted on in R.

# Sections of the program which require input from the user are bracketed by

# ??????????????????????????????????????????????????

# The user should have the following files in the same directory:
# one program file NOb.R,
# one data file association.txt
# one availability file habitat.txt
# one data file pfar.txt
# one data file position.txt
# and one file of R functions, NicheFunctions.txt

# This illustration examines mean niche overlap between five species
# of reef fishes based on three functional traits (1) Habitat association 
# (resource selection); (2) Position within a lagoon (categorical 
# data); and (3) Pectoral fin aspect ratio (measurement data).

# In this example, different species are compared for niche overlaps.
# The basic data file accosiation.txt has its first three columns labelled site, id (for individual)
# and species. Subsequent columns are the variables of different types which were measured (in this case 
# only habitat association). The data files pfar.txt and position.txt have their first two columns 
# labelled id (for individual) and species. Subsequent columns are the variables of different types 
# which were measured.

# Since one of the variables, habitat association, is of the Òresource selectionÓ type, we must also
# supply a matching availability file habitat.txt. There must be
# one availability file for each resource-selection type of variable.

# The variables may be habitat where found, food eaten, morphological columns, etc.
# They may refer to resource usage, or be categories, counts, or continuous data. If continuous,
# they may be measurements or ratios of measurement. The possible types of
# variables are:
# ÓcatÓ = categorical, but not resource selection
# ÓbinÓ = binary
# ÓctsÓ = continuous, use raw data (no transformation)
# ÓmesaÓ = measurement, continuous positive, take logs
# ÓpentÓ = percentage data, bounds at 0 and 100, use logits
# ÓpropaneÓ = proportion data, bounds at 0 and 1, use logits
# ÓcountÓ = count
# ÓrselÓ = resource selection, categorical

# The variable types need to be specified while running the program.
# The commands in the program file NOb.txt are to be copied and pasted into
# R, using the directory which contains the program and data files. At some stages in
# the program file, the user needs to type in appropriate details, e.g. specifying the types of
# variables. 

# The analysis is divided into seven distinct sections, 
# each identified by
#############################################################
# 1. Analysis A: habitat association
# 2. Analysis B: position within a lagoon
# 3. Analysis C: pectoral fin aspect ratio
# 4. The combination of the three axes into a single array
#    for oberved data, and for pseudo data, and the calulation
#    of mean niche overlap across axes
# 5. Tests determining if two species occupy different niche space,
#    seperately for each each axis, and for niche overlap averaged
#    across axes.
# 6. Tests determining if the species are evenly distributed across,
#    or clumped within niche space. This is done seperately for each
#    each axis, and for niche overlap averaged across axes.
# 7. The compilation of the above results into a single R object
#    for ease of comparison/printing/saving

Sections of the program which require input from the user
# are bracketed by

# ??????????????????????????????????????????????????


# When supplying your own data, make the first two columns id = individual and
# species = species or some other taxonomic group. The data should have the column
# names across the top, and thereafter one row per individual.
# For any resource type of variable, there must be an associated file giving the availabilities.
# The first row should be the names of the choices, and the second row is either
# the percentages or the proportions of the different choices. 
# There are comments throughout the program to explain what is being calculated at
# each stage.
