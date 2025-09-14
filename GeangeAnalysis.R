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
library(hrbrthemes)
library(gridExtra)
library(grid)

#import and prep data 
#import data
WideData <- read_excel("C:/Users/jackb/Downloads/Final_Cleaned_Pools_Dataset_11_10.xlsx", 
                       sheet = "FINAL_DATA_W_OTHERDATA")
#flip data
LongData <- gather(data = WideData, key = "Species", value = "Count", "Ostracod":"Springtail") 

#turn NA strings into real R na values
LongData <- LongData %>%
  mutate(Count = na_if(Count, "NA"))%>% 
  mutate(Depth_cm = na_if(Depth_cm, "NA")) %>%
  mutate(Temperature_C = na_if(Temperature_C, "NA")) %>% 
  mutate(Dissolved_Oxygen = na_if(Dissolved_Oxygen, "NA"))

#check to make sure there aren't any NAs remaining in other columns
LongDataCheck1 <- LongData[, !(names(LongData) %in% c("Count", "Depth_cm", "Temperature_C", "Dissolved_Oxygen", "NOTES"))]
sum(is.na(LongDataCheck1))
#we good

#turn Count column into numeric
LongData$Count <- as.numeric(LongData$Count)
#note there are now going to be some NAs where data was missing 

#remove any sampling sequence other than the first of each month 
#define sequences to remove
other_sequences = c("2nd Feb", "2nd Mar", "3rd Mar", "4th Mar", "2nd Apr", "2nd May") 
#use filter function to remove them
LongData <- LongData %>% 
  filter(!(Sampling_Sequence %in% other_sequences))
#same process to remove pond data 
LongData <- LongData %>% 
  filter(!(Pool_Number %in% "Pond"))
