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
    #install.packages("hrbrthemes")
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
    install.packages("stringr")
    install.packages("patchwork")
  
  #call packages
    library(tidyverse)
    library(vegan)
    library(ggplot2)
    library(dplyr)
    library(ggridges) 
    library(readxl) 
    #library(systemfonts)
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
    library(stringr)
    library(patchwork)

  #import and prep data 
    #ensure wd is correct
      setwd(here())
    #import data
      WideData <- read_excel(here("FINAL_POOLS_DATA_(5-26).xlsx"))
    #flip data
      LongData <- gather(data  = WideData, 
                   key   = "Species", 
                   value = "Count", "Ostracod":"Springtail") 

    #turn NA strings into real R na values
      LongData <- LongData |> 
        mutate(Count            = na_if(Count, "NA"))    |> 
        mutate(Depth_cm         = na_if(Depth_cm, "NA")) |> 
        mutate(Temp_C    = na_if(Temp_C, "NA")) |> 
        mutate(`O2_%` = na_if(`O2_%`, "NA"))
      
    LongData <-  rename(LongData, O2 = `O2_%`)

  #check to make sure there aren't any NAs remaining in other columns
    sum(is.na(LongData[, !(names(LongData) %in% c("Count", "Depth_cm", 
                                                  "Temp_C", 
                                                  "O2",
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
                 "Ostracod", "Nematode", "Springtail") 
  
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
    LongDataPools$Temp_C <- as.numeric(LongDataPools$Temp_C)
    LongDataPools$O2   <- as.numeric(LongDataPools$O2)
  
#SECTION 2: ABUNDANCE ----------------------------------------------------------
  ## 2.1: Pie chart of relative abundances -------------------------------------
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
    
  ## 2.2: Total Abundance Over Time (Focal Groups) -----------------------------
    #Find total abundance for each month (doing monthly abundance because 
    #theoretically sampling should be representative of each month)
  
    TotalAbun <- LongDataPools |> 
      mutate(
        year               = year(ymd(Date_Sampled)),
        month              = str_remove(Sampling_Sequence, "^1st "),
        month_num          = match(month, month.abb)
      )                             |> 
      group_by(year, month_num)     |> 
      summarise(
        total_observations = n(),
        avg_julian_day = mean(Julian_Day, na.rm = TRUE),
        .groups            = "drop"
      )                             |> 
      arrange(year, month_num)      |> 
      mutate(
        year_month = make_date(year, month_num, 1)
      ) |> 
      add_row(year = 2022, month_num = 8, total_observations = 0, avg_julian_day = 208) |> 
      add_row(year = 2023, month_num = 2, total_observations = 0, avg_julian_day = 404) |> 
      add_row(year = 2023, month_num = 6, total_observations = 0, avg_julian_day = 521) |> 
      arrange(year, month_num)      |> 
      # complete(
      #   year_month = seq(
      #     as.Date("2022-03-01"),  # explicitly start at March 2022
      #     as.Date("2023-09-01"),  # explicitly end at September 2023
      #     by = "month"
      #   ),
      #   fill = list(total_observations = 0)
      # ) |> 
      head(18)
      
      sample_jday <- TotalAbun$avg_julian_day
      
    TotalAbun.sd <- LongDataPools |> 
      mutate(
        year       = year(ymd(Date_Sampled)),
        month      = str_remove(Sampling_Sequence, "^1st "),
        month_num  = match(month, month.abb),
        year_month = make_date(year, month_num, 1)
      ) |> 
      group_by(Pool_Number, year_month) |> 
      summarise(
        monthly_total = n(),
        .groups = "drop"
      ) |> 
      group_by(year_month) |> 
      summarise(
        abun.sd = sd(monthly_total, na.rm = TRUE), 
        .groups = "drop"
      ) |> 
      complete(
        year_month = seq(
          as.Date("2022-03-01"),  # explicitly start at March 2022
          as.Date("2023-09-01"),  # explicitly end at September 2023
          by = "month"
        ),
        fill = list(abun.sd = 0)
      ) |> 
      head(18)
    
    TotalAbun <- TotalAbun |> 
      mutate(
        abun.sd = TotalAbun.sd$abun.sd
      )
    
  #plot those abundances 
   TotalAbun.plot <- ggplot(TotalAbun, aes(x = avg_julian_day, 
                                                y = total_observations)) +
      geom_line()  +
      geom_point() +
      geom_errorbar(
       aes(
         ymin = total_observations - abun.sd,
         ymax = total_observations + abun.sd
       ),
       width = 10
      ) +
     scale_x_continuous(
       limits = c(0,625),
       breaks = seq(0,600, by = 100)
     ) +
      labs(
        x     = "Days Since January 1st, 2022",
        y     = "Total Organismal Abundance",
        title = NULL
      ) +
      theme_minimal()
  
  ##2.3: Abundance by Pool Over Time (Focal Groups) ----------------------------
  PoolNum <- 1
   PoolAbun <- LongDataPools |> 
     filter(Pool_Number == PoolNum) %>%
     mutate(
       year       = year(ymd(Date_Sampled)),
       month      = str_remove(Sampling_Sequence, "^1st "),
       month_num  = match(month, month.abb),
       year_month = make_date(year, month_num, 1)
     ) |> 
     group_by(year_month) %>%
     summarise(total_observations = n(), .groups = "drop") |> 
     complete(
       year_month = seq(
         as.Date("2022-03-01"),  # explicitly start at March 2022
         as.Date("2023-09-01"),  # explicitly end at September 2023
         by = "month"
       ),
       fill = list(total_observations = 0)
     ) |> 
     head(18) |> 
     mutate(
       year      = year(year_month),
       month_num = month(year_month),
       month     = month(year_month, label = TRUE, abbr = TRUE),
      julian_day = sample_jday
     )
   
   
  #Function to do this for each pool individually
   AbunPlot <- function(PoolNum) {
    
     PoolAbun <- LongDataPools |> 
       filter(Pool_Number == PoolNum) %>%
       mutate(
         year       = year(ymd(Date_Sampled)),
         month      = str_remove(Sampling_Sequence, "^1st "),
         month_num  = match(month, month.abb),
         year_month = make_date(year, month_num, 1)
       ) |> 
       group_by(year_month) %>%
       summarise(total_observations = n(), .groups = "drop") |> 
       complete(
         year_month = seq(
           as.Date("2022-03-01"),  # explicitly start at March 2022
           as.Date("2023-09-01"),  # explicitly end at September 2023
           by = "month"
         ),
         fill = list(total_observations = 0)
       ) |> 
       head(18) |> 
       mutate(
         year      = year(year_month),
         month_num = month(year_month),
         month     = month(year_month, label = TRUE, abbr = TRUE),
         julian_day = sample_jday
       )
  
  abun.sd <- sd(PoolAbun$total_observations)
   
  PoolAbun.plot <- ggplot(PoolAbun, aes(x = julian_day, 
                                               y = total_observations)) +
    geom_line()  +
    geom_point() +
    geom_errorbar(
      aes(
        ymin = pmax(0, total_observations - abun.sd),
        ymax = total_observations + abun.sd
      ),
      width = 10
    ) +
    scale_y_continuous(
      limits = c(0,700), 
      breaks = seq(0,600, by = 200)
    ) +
    labs(
      title = paste("Pool", PoolNum),
      x = NULL, 
      y = NULL
    ) +
    theme_minimal()
  
    return(PoolAbun.plot)
   }
   
   
   #use this function on each pool individually 
   AbunPlot.list <- list(AbunPlot(1), AbunPlot(2), AbunPlot(3), AbunPlot(4), 
                     AbunPlot(5), AbunPlot(6), AbunPlot(7), AbunPlot(8),
                     AbunPlot(9), AbunPlot(10))
   
   AbunPanel <- wrap_plots(AbunPlot.list, nrow = 5, ncol = 2) 
   
   
  ##2.4: Total Abundance (All Groups) ------------------------------------------
   #create new version of LongData, but without filtering specific species
   WideData <- rename(WideData, O2 = `O2_%`)
   
   LDAll <- gather(data  = WideData, 
                   key   = "Species", 
                   value = "Count", "Ostracod":"Springtail") 
   
   #turn NA strings into real R na values
   LDAll <- LDAll |> 
     mutate(Count            = na_if(Count, "NA"))    |> 
     mutate(Depth_cm         = na_if(Depth_cm, "NA")) |> 
     mutate(Temperature_C    = na_if(Temp_C, "NA")) |> 
     mutate(Dissolved_Oxygen = na_if(O2, "NA"))
   
   #check to make sure there aren't any NAs remaining in other columns
   sum(is.na(LDAll[, !(names(LDAll) %in% c("Count", "Depth_cm", 
                                           "Temperature_C", 
                                           "Dissolved_Oxygen",
                                           "NOTES"))]))
   #we good if you see 0
   
   #turn Count column into numeric
   LDAll$Count <- as.numeric(LDAll$Count)
   #note there are now going to be some NAs where data was missing 
   
   #remove any sampling sequence other than the first of each month 
   #define sequences to remove
   other_sequences = c("2nd Feb", "2nd Mar", "3rd Mar", "4th Mar", "2nd Apr", 
                       "2nd May")
   #use filter function to remove them
   LDAll <- LDAll |>  
     filter(!(Sampling_Sequence %in% other_sequences))
   #same process to remove pond data 
   LDAll <- LDAll |> 
     filter(!(Pool_Number %in% "Pond"))
   
   #Expand counts, so each observation is its own line
   #convert count to integer 
   LDAll$Count <- as.integer(LDAll$Count)
   #convert NAs in count  
   LDAll$Count[is.na(LDAll$Count)] <- 0
   #expand counts 
   LDAll <- uncount(LDAll, weights = Count)
   
   #make new LDAll with just the pools (removing the creek) 
   LDAllPools <- LDAll |> 
     filter(!(Pool_Number %in% "Creek"))
   
   #add id column and move it and species to the front
   LDAllPools <- LDAllPools            |> 
     mutate(id = paste0("id", row_number())) |> 
     select(id, Species, everything())
   
   
   #convert Depth_cm NAs into 0
   LDAllPools$Depth_cm[is.na(LDAllPools$Depth_cm)] <- 0
   
   #convert columns into number 
   LDAllPools$Depth_cm <- as.numeric(LDAllPools$Depth_cm)
   LDAllPools$Temp_C <- as.numeric(LDAllPools$Temp_C)
   LDAllPools$O2   <- as.numeric(LDAllPools$O2)
   
   
  #Repeat 2.2 with LDAll
   TotalAbunAll <- LDAllPools |> 
     mutate(
       year               = year(ymd(Date_Sampled)),
       month              = str_remove(Sampling_Sequence, "^1st "),
       month_num          = match(month, month.abb)
     )                             |> 
     group_by(year, month_num)     |> 
     summarise(
       total_observations = n(),
       .groups            = "drop"
     )                             |> 
     arrange(year, month_num)      |> 
     mutate(
       year_month = make_date(year, month_num, 1)
     ) |> 
     complete(
       year_month = seq(
         as.Date("2022-03-01"),  # explicitly start at March 2022
         as.Date("2023-09-01"),  # explicitly end at September 2023
         by = "month"
       ),
       fill = list(total_observations = 0)
     ) |> 
     head(18)
   
   TotalAbunAll.sd <- LDAllPools |> 
     mutate(
       year       = year(ymd(Date_Sampled)),
       month      = str_remove(Sampling_Sequence, "^1st "),
       month_num  = match(month, month.abb),
       year_month = make_date(year, month_num, 1)
     ) |> 
     group_by(Pool_Number, year_month) |> 
     summarise(
       monthly_total = n(),
       .groups = "drop"
     ) |> 
     group_by(year_month) |> 
     summarise(
       abun.sd = sd(monthly_total, na.rm = TRUE), 
       .groups = "drop"
     ) |> 
     complete(
       year_month = seq(
         as.Date("2022-03-01"),  # explicitly start at March 2022
         as.Date("2023-09-01"),  # explicitly end at September 2023
         by = "month"
       ),
       fill = list(abun.sd = 0)
     ) |> 
     head(18)
   
   TotalAbunAll <- TotalAbunAll |> 
     mutate(
       abun.sd = TotalAbunAll.sd$abun.sd
     )
   
   TotalAbunAll.plot <- ggplot(TotalAbunAll, aes(x = year_month, 
                                                y = total_observations)) +
     geom_line()  +
     geom_point() +
     geom_errorbar(
       aes(
         ymin = total_observations - abun.sd,
         ymax = total_observations + abun.sd
       ),
       width = 10
     ) +
     scale_x_date(
       breaks = seq(
         from = as.Date("2022-03-01"),
         to   = max(TotalAbunAll$year_month),
         by   = "3 months"
       ),
       date_labels = "%Y-%b") +
     labs(
       x     = "Month",
       y     = "Total Number of Organisms",
       title = "Total Pool Abundance Over Time (All Groups)"
     ) +
     theme_minimal()
  
  ##2.5: Abundance by Pool (All Groups) ---------------------------------------- 
   AbunPlotAll <- function(PoolNum) {
     
     PoolData <- LDAllPools |> 
       filter(Pool_Number == PoolNum) %>%
       mutate(
         year       = year(ymd(Date_Sampled)),
         month      = str_remove(Sampling_Sequence, "^1st "),
         month_num  = match(month, month.abb),
         year_month = make_date(year, month_num, 1)
       ) |> 
       group_by(year_month) %>%
       summarise(total_observations = n(), .groups = "drop") |> 
       complete(
         year_month = seq(
           as.Date("2022-03-01"),  # explicitly start at March 2022
           as.Date("2023-09-01"),  # explicitly end at September 2023
           by = "month"
         ),
         fill = list(total_observations = 0)
       ) |> 
       mutate(
         year      = year(year_month),
         month_num = month(year_month),
         month     = month(year_month, label = TRUE, abbr = TRUE),
         julian_day = sample_jday
       ) 
     
     abun.sd <- sd(PoolData$total_observations)
     
     PoolAbun.plot <- ggplot(PoolData, aes(x = year_month, 
                                           y = total_observations)) +
       geom_line()  +
       geom_point() +
       geom_errorbar(
         aes(
           ymin = pmax(0, total_observations - abun.sd),
           ymax = total_observations + abun.sd
         ),
         width = 10
       ) +
       scale_x_date(
         breaks = seq(
           from = as.Date("2022-03-01"),
           to   = max(PoolData$year_month),
           by   = "3 months"
         ),
         date_labels = "%Y-%b") +
       scale_y_continuous(
         limits = c(0,700), 
         breaks = seq(0,600, by = 100)
       ) +
       labs(
         x     = "Month",
         y     = "Total Number of Organisms",
         title = paste("Pool", PoolNum, "Abundance Over Time (All Groups)")
       ) +
       theme_minimal()
     
     return(PoolAbun.plot)
   }
   
   AbunPlotAll(1)
   AbunPlotAll(2)
   AbunPlotAll(3)
   AbunPlotAll(4)
   AbunPlotAll(5)
   AbunPlotAll(6)
   AbunPlotAll(7)
   AbunPlotAll(8)
   AbunPlotAll(9) 
   AbunPlotAll(10)
   

  #SECTION 3: RICHNESS --------------------------------------------------------- 
   ##3.1: Overall Richness (Focal Groups Only) ---------------------------------
   
   TotalRich <- LongDataPools |> 
     mutate(
       year               = year(ymd(Date_Sampled)),
       month              = str_remove(Sampling_Sequence, "^1st "),
       month_num          = match(month, month.abb)
     )                             |> 
     group_by(year, month_num)     |> 
     summarise(
       richness = n_distinct(Species),
       avg_julian_day = mean(Julian_Day, na.rm = TRUE),
       .groups            = "drop"
     )                             |> 
     arrange(year, month_num)      |> 
     mutate(
       year_month = make_date(year, month_num, 1)
     ) |> 
     add_row(year = 2022, month_num = 8, richness = 0, avg_julian_day = 208) |> 
     add_row(year = 2023, month_num = 2, richness = 0, avg_julian_day = 404) |> 
     add_row(year = 2023, month_num = 6, richness = 0, avg_julian_day = 521) |>
     arrange(year, month_num)      |> 
     head(18)
   
   TotalRich.sd <- LongDataPools |> 
     mutate(
       year       = year(ymd(Date_Sampled)),
       month      = str_remove(Sampling_Sequence, "^1st "),
       month_num  = match(month, month.abb),
       year_month = make_date(year, month_num, 1)
     ) |> 
     group_by(Pool_Number, year_month) |> 
     summarise(
       richness = n_distinct(Species),
       .groups            = "drop"
     ) |> 
     group_by(year_month) |> 
     summarise(
       rich.sd = sd(richness, na.rm = TRUE), 
       .groups = "drop"
     ) |> 
     complete(
       year_month = seq(
         as.Date("2022-03-01"),  # explicitly start at March 2022
         as.Date("2023-09-01"),  # explicitly end at September 2023
         by = "month"
       ),
       fill = list(rich.sd = 0)
     ) |> 
     head(18)
   
   TotalRich <- TotalRich |> 
     mutate(
       rich.sd = TotalRich.sd$rich.sd
     )
   
   TotalRich.plot <- ggplot(TotalRich, aes(x = avg_julian_day, 
                                                y = richness)) +
     geom_line()  +
     geom_point() +
     geom_errorbar(
       aes(
         ymin = richness - rich.sd,
         ymax = richness + rich.sd
       ),
       width = 10
     ) +
     scale_x_continuous(
       limits = c(0,625),
       breaks = seq(0,600, by = 100)
     ) +
     scale_y_continuous(
       limits = c(0,11),
       breaks = seq(0,10, by = 2)
     ) +
     labs(
       x     = "Days Since January 1st, 2022",
       y     = "Total Species Richness",
     ) +
     theme_minimal()
   
   ##3.2: Overall Richness by Pool --------------------------------------------- 
   #create a function to create the richness plot for each pool
   
   RichPlot <- function(PoolNum) {
     PoolData <- LongDataPools |> 
       filter(Pool_Number == PoolNum) |> 
       mutate(
         year       = year(ymd(Date_Sampled)),
         month      = str_remove(Sampling_Sequence, "^1st "),
         month_num  = match(month, month.abb),
         year_month = make_date(year, month_num, 1)
       ) |> 
       group_by(year_month) |> 
       summarise(
         richness = n_distinct(Species),
         .groups            = "drop"
       )  |> 
       complete(
         year_month = seq(
           as.Date("2022-03-01"),  # explicitly start at March 2022
           as.Date("2023-09-01"),  # explicitly end at September 2023
           by = "month"
         ),
         fill = list(richness = 0)
       ) |> 
       head(18) |> 
       mutate(
         year      = year(year_month),
         month_num = month(year_month),
         month     = month(year_month, label = TRUE, abbr = TRUE),
         julian_day = sample_jday
       )
   
     rich.sd <- sd(PoolData$richness)
     
     PoolRich.plot <- ggplot(PoolData, aes(x = julian_day, 
                                                   y = richness)) +
       geom_line()  +
       geom_point() +
       geom_errorbar(
         aes(
           ymin = pmax(0, richness - rich.sd),
           ymax = pmin(richness + rich.sd, 9)
         ),
         width = 10
       ) +
       scale_y_continuous(
         limits = c(0,12), 
         breaks = seq(0,15, by = 5)
       ) +
       labs(
         x     = NULL,
         y     = NULL,
         title = paste("Pool", PoolNum)
       ) +
       theme_minimal()
     
     return(PoolRich.plot)
   }
   
   RichPlot.list <- list(RichPlot(1), RichPlot(2), RichPlot(3), RichPlot(4), 
                         RichPlot(5), RichPlot(6), RichPlot(7), RichPlot(8),
                         RichPlot(9), RichPlot(10))
   
   RichPanel <- wrap_plots(RichPlot.list, nrow = 5, ncol = 2) 

   
   
   
  ##3.3: Total Richness (All Groups) -------------------------------------------
  
  #Repeat 3.1 with all species dataset 
     TotalRichAll <- LDAllPools |> 
       mutate(
         year               = year(ymd(Date_Sampled)),
         month              = str_remove(Sampling_Sequence, "^1st "),
         month_num          = match(month, month.abb)
       )                             |> 
       group_by(year, month_num)     |> 
       summarise(
         richness = n_distinct(Species),
         .groups            = "drop"
       )                             |> 
       arrange(year, month_num)      |> 
       mutate(
         year_month = make_date(year, month_num, 1)
       ) |> 
       complete(
         year_month = seq(
           as.Date("2022-03-01"),  # explicitly start at March 2022
           as.Date("2023-09-01"),  # explicitly end at September 2023
           by = "month"
         ),
         fill = list(richness = 0)
       ) |> 
        head(18)
   
   TotalRichAll.sd <- LDAllPools |> 
     mutate(
       year       = year(ymd(Date_Sampled)),
       month      = str_remove(Sampling_Sequence, "^1st "),
       month_num  = match(month, month.abb),
       year_month = make_date(year, month_num, 1)
     ) |> 
     group_by(Pool_Number, year_month) |> 
     summarise(
       richness = n_distinct(Species),
       .groups = "drop"
     ) |> 
     group_by(year_month) |> 
     summarise(
       rich.sd = sd(richness, na.rm = TRUE), 
       .groups = "drop"
     ) |> 
     complete(
       year_month = seq(
         as.Date("2022-03-01"),  # explicitly start at March 2022
         as.Date("2023-09-01"),  # explicitly end at September 2023
         by = "month"
       ),
       fill = list(abun.sd = 0)
     ) |> 
     head(18)
   
   TotalRichAll <- TotalRichAll |> 
     mutate(
       rich.sd = TotalRichAll.sd$rich.sd
     )
     
     TotalRichALL.plot <- ggplot(TotalRichAll, aes(x = year_month, 
                                                    y = richness)) +
       geom_line()  +
       geom_point() +
       geom_errorbar(
         aes(
           ymin = pmax(0, richness - rich.sd),
           ymax = richness + rich.sd
         ),
         width = 10
       ) +
       scale_x_date(
         breaks = seq(
           from = as.Date("2022-03-01"),
           to   = max(TotalRichAll$year_month),
           by   = "3 months"
         ),
         date_labels = "%Y-%b") + 
       scale_y_continuous(
         limits = c(0,22), 
         breaks = seq(0,22, by = 5)
       ) +
       labs(
         x     = "Month",
         y     = "Total Richness",
         title = "Total Pool Richness Over Time (All Groups)"
       ) +
       theme_minimal()
   
  ## 3.4: Richness by Pool (All groups) ----------------------------------------------
    #repeat 3.2 with LDAll dataset
     
     RichPlotAll <- function(PoolNum) {
       PoolData <- LDAllPools |> 
         filter(Pool_Number == PoolNum) |> 
         mutate(
           year       = year(ymd(Date_Sampled)),
           month      = str_remove(Sampling_Sequence, "^1st "),
           month_num  = match(month, month.abb),
           year_month = make_date(year, month_num, 1)
         ) |> 
         group_by(year_month) |> 
         summarise(
           richness = n_distinct(Species),
           .groups            = "drop"
         )  |> 
         complete(
           year_month = seq(
             as.Date("2022-03-01"),  # explicitly start at March 2022
             as.Date("2023-09-01"),  # explicitly end at September 2023
             by = "month"
           ),
           fill = list(richness = 0)
         ) |> 
         mutate(
           year      = year(year_month),
           month_num = month(year_month),
           month     = month(year_month, label = TRUE, abbr = TRUE)
         )
       
       rich.sd <- sd(PoolData$richness)
       
       PoolRich.plot <- ggplot(PoolData, aes(x = year_month, 
                                             y = richness)) +
         geom_line()  +
         geom_point() +
         geom_errorbar(
           aes(
             ymin = pmax(0, richness - rich.sd),
             ymax = richness + rich.sd
           ),
           width = 10
         ) +
         scale_x_date(
           breaks = seq(
             from = as.Date("2022-03-01"),
             to   = max(PoolData$year_month),
             by   = "3 months"
           ),
           date_labels = "%Y-%b") + 
         scale_y_continuous(
           limits = c(0,20), 
           breaks = seq(0,20, by = 5)
         ) +
         labs(
           x     = "Month",
           y     = "Total Richness",
           title = paste("Pool", PoolNum, "Richness Over Time (All Groups)")
         ) +
         theme_minimal()
       
       return(PoolRich.plot)
     }
     
     RichPlotAll(1)
     RichPlotAll(2)
     RichPlotAll(3)
     RichPlotAll(4)   
     RichPlotAll(5)
     RichPlotAll(6)
     RichPlotAll(7)
     RichPlotAll(8)
     RichPlotAll(9)
     RichPlotAll(10) 

     
     
     
# SECTION 4: DIVERSITY ---------------------------------------------------------
  ##4.1: Total Diversity Over Time (Focal Groups) ------------------------------
 
  #filter unneeded data out of WideData
  WDna <- WideData |> 
       #use mutate across everything to turn all instances of "NA" into a 0
       mutate(across(12:20, ~ ifelse(. == "NA", 0, .))) |> 
       #convert all data columns into numeric 
       mutate(across(12:20, as.numeric)) |> 
       #remove everything that's not considered
       filter(!(Sampling_Type %in% "Surface skim")) |> 
       filter(!(Pool_Number %in% "Pond")) |>
       filter(!(Pool_Number %in% "Creek")) |> 
       filter(!(Sampling_Sequence %in% other_sequences)) |> 
       select(
         1:11,
         any_of(species05)
       ) 
     
    #calculate diversity
     div <- diversity(WDna[,12:20], index = "shannon")
     even <- div/log(specnumber(WDna[,12:20]))
     
     #make it a new column in the dataset
     WDna$shannon_div <- div
     WDna$even <- even
       
    # make the sampling sequences work and then summarize by diversity
     WDna.div <- WDna |> 
       mutate(
         year               = year(ymd(Date_Sampled)),
         month              = str_remove(Sampling_Sequence, "^1st "),
         month_num          = match(month, month.abb)
       )                             |> 
       group_by(year, month_num)     |> 
       summarise(
         diversity = mean(shannon_div),
         avg_julian_day = mean(Julian_Day, na.rm = TRUE),
         .groups            = "drop"
       )                             |> 
       arrange(year, month_num)      |> 
       mutate(
         year_month = make_date(year, month_num, 1)
       ) |> 
       head(18)
     
    WDna.div.sd <- WDna |>  
      mutate(
        year       = year(ymd(Date_Sampled)),
        month      = str_remove(Sampling_Sequence, "^1st "),
        month_num  = match(month, month.abb),
        year_month = make_date(year, month_num, 1)
      ) |> 
      group_by(Pool_Number, year_month) |> 
      summarise(
        diversity = mean(shannon_div),
        .groups            = "drop"
      ) |> 
      group_by(year_month) |> 
      summarise(
        div.sd = sd(diversity, na.rm = TRUE), 
        .groups = "drop"
      ) |> 
      complete(
        year_month = seq(
          as.Date("2022-03-01"),  # explicitly start at March 2022
          as.Date("2023-09-01"),  # explicitly end at September 2023
          by = "month"
        ),
        fill = list(div.sd = 0)
      ) |> 
      head(18)
    
    WDna.div <- WDna.div |> 
      mutate(
        div.sd = WDna.div.sd$div.sd
      )
    
    WDna.div$diversity[is.na(WDna.div$diversity)] <- 0
    
    #plot it
     TotalDiv.plot <- ggplot(WDna.div, aes(x = avg_julian_day, y = diversity)) +
       geom_line()  +
       geom_point() +
       geom_errorbar(
         aes(
           ymin = pmax(0, diversity - div.sd),
           ymax = diversity + div.sd
         ),
         width = 10
       ) +
       scale_y_continuous(
         limits = c(0,1), 
         breaks = seq(0,1, by = 0.2)
       ) +
       scale_x_continuous(
         limits = c(0,625),
         breaks = seq(0,600, by = 100)
       ) +
       labs(
         x     = "Days Since January 1st, 2022",
         y     = "Shannon Diversity",
       ) +
       theme_minimal() 
     
  ## 4.2: Diversity by Pool (Focal Groups) -------------------------------------
     #create a function to create the above graph for each pool individually

    DivPlot <- function(PoolNum) {

    Pool.div <- WDna |>  
      filter(Pool_Number == PoolNum) |> 
      mutate(
        year               = year(ymd(Date_Sampled)),
        month              = str_remove(Sampling_Sequence, "^1st "),
        month_num          = match(month, month.abb)
      )                             |> 
      group_by(year, month_num)     |> 
      summarise(
        diversity = mean(shannon_div),
        .groups            = "drop"
      )                             |> 
      arrange(year, month_num)      |> 
      head(18)                      |> 
      mutate(
        year_month = make_date(year, month_num, 1),
        julian_day = sample_jday
      )
    
    Pool.div$diversity[is.na(Pool.div$diversity)] <- 0
    
    div.sd <- sd(Pool.div$diversity)
    
     PoolDiv.plot <- ggplot(Pool.div, aes(x = julian_day, y = diversity)) +
      geom_line()  +
      geom_point() +
       geom_errorbar(
         aes(
           ymin = pmax(0, diversity - div.sd),
           ymax = pmin(diversity + div.sd, 1)
         ),
         width = 10
       ) +
      scale_y_continuous(
        limits = c(0,1.1), 
        breaks = seq(0,1, by = 0.5)
      ) +
      labs(
        x     = NULL,
        y     = NULL,
        title = paste("Pool", PoolNum)
      ) +
      theme_minimal() 
     
     return(PoolDiv.plot)
     
    }
    
     DivPlot.list <- list(DivPlot(1), DivPlot(2), DivPlot(3), DivPlot(4), 
                           DivPlot(5), DivPlot(6), DivPlot(7), DivPlot(8),
                           DivPlot(9), DivPlot(10))
     
     DivPanel <- wrap_plots(DivPlot.list, nrow = 5, ncol = 2) 
    
  
  ##4.3: Total Diversity Over Time (All Groups) --------------------------------
    
    WDna.All <- WideData |> 
      #use mutate across everything to turn all instances of "NA" into a 0
      mutate(across(17:50, ~ ifelse(. == "NA", 0, .))) |> 
      #convert all data columns into numeric 
      mutate(across(17:50, as.numeric)) |> 
      #remove everything that's not considered
      filter(!(Sampling_Type %in% "Surface skim")) |> 
      filter(!(Pool_Number %in% "Pond")) |>
      filter(!(Pool_Number %in% "Creek")) |> 
      filter(!(Sampling_Sequence %in% other_sequences)) 
    
    #calculate diversity
    div <- diversity(WDna.All[,17:50], index = "shannon")
    even <- div/log(specnumber(WDna.All[,17:50]))
    
    #make it a new column in the dataset
    WDna.All$shannon_div <- div
    WDna.All$even <- even
    
    # make the sampling sequences work and then summarize by diversity
    WDna.All.div <- WDna.All |> 
      mutate(
        year               = year(ymd(Date_Sampled)),
        month              = str_remove(Sampling_Sequence, "^1st "),
        month_num          = match(month, month.abb)
      )                             |> 
      group_by(year, month_num)     |> 
      summarise(
        diversity = mean(shannon_div),
        .groups            = "drop"
      )                             |> 
      arrange(year, month_num)      |> 
      mutate(
        year_month = make_date(year, month_num, 1)
      ) |> 
      head(18)
    
    WDna.All.div.sd <- WDna.All |>  
      mutate(
        year       = year(ymd(Date_Sampled)),
        month      = str_remove(Sampling_Sequence, "^1st "),
        month_num  = match(month, month.abb),
        year_month = make_date(year, month_num, 1)
      ) |> 
      group_by(Pool_Number, year_month) |> 
      summarise(
        diversity = mean(shannon_div),
        .groups            = "drop"
      ) |> 
      group_by(year_month) |> 
      summarise(
        div.sd = sd(diversity, na.rm = TRUE), 
        .groups = "drop"
      ) |> 
      complete(
        year_month = seq(
          as.Date("2022-03-01"),  # explicitly start at March 2022
          as.Date("2023-09-01"),  # explicitly end at September 2023
          by = "month"
        ),
        fill = list(div.sd = 0)
      ) |> 
      head(18)
    
    WDna.All.div <- WDna.All.div |> 
      mutate(
        div.sd = WDna.All.div.sd$div.sd
      )
    
    WDna.All.div$diversity[is.na(WDna.All.div$diversity)] <- 0
    
    #plot
    DivAll.plot <- ggplot(WDna.All.div, aes(x = year_month, y = diversity)) +
      geom_line()  +
      geom_point() +
      geom_errorbar(
        aes(
          ymin = pmax(0, diversity - div.sd),
          ymax = diversity + div.sd
        ),
        width = 10
      ) +
      scale_x_date(
        breaks = seq(
          from = as.Date("2022-03-01"),
          to   = max(WDna.All.div$year_month),
          by   = "3 months"
        ),
        date_labels = "%Y-%b") + 
      scale_y_continuous(
        limits = c(0,1.1), 
        breaks = seq(0,1.1, by = 0.2)
      ) +
      labs(
        x     = "Month",
        y     = "Diversity",
        title = "Shannon Diversity Over Time \n Averaged Across All Pools (All Groups)"
      ) +
      theme_minimal() 
    
    ##4.4: Diversity by Pool (All Groups) --------------------------------------
    
    DivPlotAll <- function(PoolNum) {
  
      Pool.div <- WDna.All |>  
        filter(Pool_Number == PoolNum) |> 
        mutate(
          year               = year(ymd(Date_Sampled)),
          month              = str_remove(Sampling_Sequence, "^1st "),
          month_num          = match(month, month.abb)
        )                             |> 
        group_by(year, month_num)     |> 
        summarise(
          diversity = mean(shannon_div),
          .groups            = "drop"
        )                             |> 
        arrange(year, month_num)      |> 
        mutate(
          year_month = make_date(year, month_num, 1)
        ) |> 
        head(18)
      
      # convert NAs to 0
      Pool.div$diversity[is.na(Pool.div$diversity)] <- 0
      
      div.sd <- sd(Pool.div$diversity)
      
      PoolDiv.plot <- ggplot(Pool.div, aes(x = year_month, y = diversity)) +
        geom_line()  +
        geom_point() +
        geom_errorbar(
          aes(
            ymin = pmax(0, diversity - div.sd),
            ymax = diversity + div.sd
          ),
          width = 10
        ) +
        scale_x_date(
          breaks = seq(
            from = as.Date("2022-03-01"),
            to   = max(WDna.div$year_month),
            by   = "3 months"
          ),
          date_labels = "%Y-%b") + 
        scale_y_continuous(
          limits = c(0,1.65), 
          breaks = seq(0,1.6, by = 0.2)
        ) +
        labs(
          x     = "Month",
          y     = "Diversity",
          title = paste("Pool", PoolNum, "Shannon Diversity Over Time \n All Groups")
        ) +
        theme_minimal() 
      
      return(PoolDiv.plot)
      
    }

    
    #use function on each pool
    DivPlotAll(1)
    DivPlotAll(2)
    DivPlotAll(3)
    DivPlotAll(4)   
    DivPlotAll(5)
    DivPlotAll(6)
    DivPlotAll(7)
    DivPlotAll(8)
    DivPlotAll(9)
    DivPlotAll(10)
    
# SECTION 5: EVENNESS ----------------------------------------------------------
  ##5.1: Total Evenness Over Time ----------------------------------------------
    
    WDna <- WideData |> 
      #use mutate across everything to turn all instances of "NA" into a 0
      mutate(across(17:50, ~ ifelse(. == "NA", 0, .))) |> 
      #convert all data columns into numeric 
      mutate(across(17:50, as.numeric)) |> 
      #remove everything that's not considered
      filter(!(Sampling_Type %in% "Surface skim")) |> 
      filter(!(Pool_Number %in% "Pond")) |>
      filter(!(Pool_Number %in% "Creek")) |> 
      filter(!(Sampling_Sequence %in% other_sequences)) |> 
      select(
        1:16,
        any_of(species05)
      ) 
    
    #calculate evenness
    div <- diversity(WDna[,17:25], index = "shannon")
    even <- div/log(9)
    
    #make it a new column in the dataset
    WDna$shannon_div <- div
    WDna$even <- even
    
    # make the sampling sequences work and then summarize by evenness
    WDna.even <- WDna |> 
      mutate(
        year               = year(ymd(Date_Sampled)),
        month              = str_remove(Sampling_Sequence, "^1st "),
        month_num          = match(month, month.abb)
      )                             |> 
      group_by(year, month_num)     |> 
      summarise(
        evenness = mean(even),
        avg_julian_day = mean(Julian_Day, na.rm = TRUE),
        .groups            = "drop"
      )                             |> 
      arrange(year, month_num)      |> 
      mutate(
        year_month = make_date(year, month_num, 1)
      ) |> 
      head(18)
    
    WDna.even.sd <- WDna |>  
      mutate(
        year       = year(ymd(Date_Sampled)),
        month      = str_remove(Sampling_Sequence, "^1st "),
        month_num  = match(month, month.abb),
        year_month = make_date(year, month_num, 1)
      ) |> 
      group_by(Pool_Number, year_month) |> 
      summarise(
        evenness = mean(even),
        .groups            = "drop"
      ) |> 
      group_by(year_month) |> 
      summarise(
        even.sd = sd(even, na.rm = TRUE), 
        .groups = "drop"
      ) |> 
      complete(
        year_month = seq(
          as.Date("2022-03-01"),  # explicitly start at March 2022
          as.Date("2023-09-01"),  # explicitly end at September 2023
          by = "month"
        ),
        fill = list(even.sd = 0)
      ) |> 
      head(18)
    
    WDna.even <- WDna.even |> 
      mutate(
        even.sd = WDna.even.sd$even.sd
      )
    
    TotalEven.plot <- ggplot(WDna.even, aes(x = avg_julian_day, y = evenness)) +
      geom_line()  +
      geom_point() +
      geom_errorbar(
        aes(
          ymin = pmax(0, evenness - even.sd),
          ymax = evenness + even.sd
        ),
        width = 10
      ) +
      scale_y_continuous(
        limits = c(0,1), 
        breaks = seq(0,1, by = 0.2)
      ) +
      scale_x_continuous(
        limits = c(0,610),
        breaks = seq(0,600, by = 100)
      ) +
      labs(
        x     = "Days Since January 1st, 2022",
        y     = "Pielou Evenness",
      ) +
      theme_minimal() 
    
  ## 5.2: Evenness by Pool -----------------------------------------------------
    PoolNum <- 6
    EvenPlot <- function(PoolNum) {
      
      Pool.even <- WDna |>  
        filter(Pool_Number == PoolNum) |> 
        mutate(
          year               = year(ymd(Date_Sampled)),
          month              = str_remove(Sampling_Sequence, "^1st "),
          month_num          = match(month, month.abb)
        )                             |> 
        group_by(year, month_num)     |> 
        summarise(
          evenness = mean(even),
          .groups            = "drop"
        )                             |> 
        arrange(year, month_num)      |> 
        head(18)                      |> 
        mutate(
          year_month = make_date(year, month_num, 1),
          julian_day = sample_jday
        )
      
      print(Pool.even)
      
      Pool.even$evenness[is.na(Pool.even$evenness)] <- 0
      
      even.sd <- sd(Pool.even$evenness)
      
      PoolEven.plot <- ggplot(Pool.even, aes(x = julian_day, y = evenness)) +
        geom_line()  +
        geom_point() +
        geom_errorbar(
          aes(
            ymin = pmax(0, evenness - even.sd),
            ymax = evenness + even.sd
          ),
          width = 10
        ) +
        scale_y_continuous(
          limits = c(0,1.1), 
          breaks = seq(0,1, by = 0.5)
        ) +
        scale_x_continuous(
          limits = c(0,610),
          breaks = seq(0,600, by = 100)
        ) +
        labs(
          x     = NULL,
          y     = NULL,
          title = paste("Pool", PoolNum)
        ) +
        theme_minimal() 
      
      return(PoolEven.plot)
      
    }
    EvenPlot(1)
    
    EvenPlot.list <- list(EvenPlot(1), EvenPlot(2), EvenPlot(3), EvenPlot(4), 
                         EvenPlot(5), EvenPlot(6), EvenPlot(7), EvenPlot(8),
                         EvenPlot(9), EvenPlot(10))
    
    EvenPanel <- wrap_plots(EvenPlot.list, nrow = 5, ncol = 2)  
    
    
    
    
    
    
    
    
    