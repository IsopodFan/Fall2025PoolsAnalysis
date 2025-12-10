#this Rscript creates panel figures for abundance, richness, and diversity 
#as well as other miscellaneous tables that don't require Geange analysis data 

#SECTION 1: DATA IMPORT, PREP, AND CLEANING ------------------------------------

#import and install packages
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
    install.packages("viridis") 
    install.packages("lubridate")
  
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
    library(lme4)
    library(lmerTest) 
    library(emmeans) 
    library(pbkrtest) 
    library(viridis) 
    library(lubridate)
    library(forcats)

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
        mutate(Count            = na_if(Count, "NA"))    |> 
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
  
  #make new LongData with just the pools (removing the creek) 
    LongDataPools <- LongData |> 
      filter(!(Pool_Number %in% "Creek"))
  
  
  #remove uncommon groups from analysis (uncommon = prevelance under 0.5%)
    #summarise data and sort by numerical prevelance
      LD_counts <- LongDataPools |> 
        group_by(Species)        |> 
        summarise(count   = n()) |> 
        mutate(percentage = (count/sum(count))*100)
  
  #7 groups with percentage > 1%
  #10 groups with percentage > 0.5%
  #18 groups with percentage > 0.1%
  
  species05 <- c("Cladoceran", "Copepod", "Diptera_Pupae", 
                 "Midge_Larvae", "Mites", "Mosquito_Larvae", 
                 "Ostracod", "Roundworm", "Springtail") 
  
  LongDataPools <- LongDataPools |> 
    filter(Species %in% species05)
  
  #add id column and move it and species to the front
    LongDataPools <- LongDataPools            |> 
      mutate(id = paste0("id", row_number())) |> 
      select(id, Species, everything())
  
  
  #convert Depth_cm NAs into 0
    LongDataPools$Depth_cm[is.na(LongDataPools$Depth_cm)] <- 0
  
  #convert columns into number 
    LongDataPools$Depth_cm <- as.numeric(LongDataPools$Depth_cm)
    LongDataPools$Temp_Reg <- as.numeric(LongDataPools$Temp_Reg)
    LongDataPools$O2_Reg   <- as.numeric(LongDataPools$O2_Reg)
  
#SECTION 2:Abundance and Richness 
  ## Pie chart of relative abundances ------------------------------------------ 
     species_counts <- LongDataPools                                |> 
      count(Species, name = "abundance")                            |> 
      mutate(
        species_grouped   = ifelse(Species %in% species05, 
                                   Species, "Other")
      )                                                             |> 
      group_by(species_grouped)                                     |> 
      summarise(abundance = sum(abundance))                         |> 
      ungroup()                                                     |> 
      mutate(
        species_grouped   = fct_reorder(species_grouped, abundance)
      )                                                             |> 
      mutate(
        species_grouped   = fct_relevel(species_grouped, "Other", 
                                        after = 0)
      )                                                             |> 
      mutate(
        species_grouped = recode(species_grouped, 
                                "Roundworm"       = "Nematodes", 
                                "Mosquito_Larvae" = "Mosquito Larvae", 
                                "Ostracod"        = "Ostracods", 
                                "Springtail"      = "Springtails", 
                                "Midge_Larvae"    = "Midge Larvae", 
                                "Diptera_Pupae"   = "Diptera Pupae", 
                                "Copepod"         = "Copepods", 
                                "Cladoceran"      = "Cladocerans")) |> 
      mutate(
        pct = abundance / sum(abundance) * 100
      )                                                             |> 
      mutate(
        cumulative = cumsum(pct),
        midpoint = cumulative - pct / 2 
      )
    
    Abun.pie <- ggplot(species_counts, aes(x    = "", 
                                           y    = abundance, 
                                           fill = species_grouped)) +
      geom_col(width                            = 1, color = "white") +
      coord_polar(theta                         = "y") +
      labs(title                                = "Relative Abundance of Prominent Taxa") +
      theme_void() +
      scale_fill_manual(values = c("gray40", "#7e2954", "#e69f00", "#56b3e9", 
                                   "#f0e442", "#009e74", "#cc79a7", "#1b1557", 
                                   "#0071b2", "#d55c00"))
    

  ## find total species richness 
    length(unique(LongData$Species))
    length(LongData$Julian_Day)
    
    
