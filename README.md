# Chinook-fecundity
Assimilation of data collected from WCVI hatchery broodstock programs enumerating fecundities of female Chinook.

Data from each year are stored in sub-folders. The principal R script (01_fecundity-data-analysis.R) draws data from the sub-folders and assimilates it into one multi-year data table, which is saved each time the script is run as "R-OUT_Fecundity-Data_all-years_{YYY-MM-DD}.xlsx"

In future years, we must strive to keep the raw data file format exactly the same. This will mitigate the need for annual updates to the script and ensure files are read properly. 
