# Stochastically Simulating Waves

## Description

future wave simulation

* reads in future cp .csv file as dataframe
* loads cpstats dicts (stats and ARDict) can look to change these to dfs
* simulates N sequences of waves M timesteps long
* returns dataframe with following columns:
    
    time | CP | season | hs | tp | direction | seqNumber

## Test Data
Test pandas dataframes that are loaded into the programme can be found in the data folder. These are:

1. hadleyFuture_rcp%s.csv
2. cpWavObs.csv
3. simDf.csv

These dataframes contain relevant data and format to run the simulation code.
