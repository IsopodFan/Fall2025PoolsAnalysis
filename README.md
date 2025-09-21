# Fall2025PoolsAnalysis

# the goal of this project is to apply the analysis in Geange et al. (2011) to the TFO pools dataset, computing the niche overlap and partitioning dynamics of present species 

description of files: 
  - ExampleAfiles, ExampleBfiles, ExampleCfiles 
     - each of these are examples provided by the Geange paper on how to use their analysis technique 

  - GeangeAnalysis.R 
     - the R script where the first analysis took place. Can be run on a fresh session

     GeangeAnalysis02.R
      - R script where second analysis took place, basically replicating the first but with the addition of the recruitment, temp, and O2 columns, and removing the creek data from analysis

   - PoolsData_Final.xlsx 
      - final dataset cleaned and ready for analysis

  - PoolsDataWReg
      - final dataset but with extra columns added for recruitment, temperature, and O2 (the latter two using regression to fill the dataset)
        
  Output files: 
    the output files are currently unorganized and lack a good naming scheme. This will be fixed later but for now I'll simply describe them here: 
      - NOestimates.xlsx and Estsandpvalues.xlsx are from the first geange analysis 
      - NOestimates02 and Estsandpvalues02 are from the second geange analysis; they disclude the creek and include recruitment, but do NOT include temp and O2 
      - NOestimatesReg and EstsandpvaluesReg are from the second geange analysis; they disclude the creek and incldue recruitment, temp, and O2
