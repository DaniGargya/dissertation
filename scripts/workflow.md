# Workflow

## Data selection/sources
- data inclusion criteria for BioTIME data set, merge data sets correctly
- standardize cell sizes of all datasets
  - partitioning BioTIME time-series? 
  - rarefaction BioTIME?
- extract/import mean accessiblity score from GEE matching the location of time-series with standardized cell size
- import population density data

## Data processing
- transformations: grid cell, coding scheme categorical value population density, bound accessiblity score between 0 and 1
- center other values of BioTIME?

## Data analysis
- calculate turnover trends using Jaccard
  - temporal: comparing last year to first year
    - 
  - spatial: using same value as TT, but relating it across space?
  
- create models
  - TT ~ A
  - ST ~ A
  - ST ~ A + A:T??
  - TT ~ A + A:HPD
  - ST ~ A + A:HPD

## Sensitivity analysis
- better temporally matched data
- sensitivity to cell sizes
- effect of latitude?
- number of studies that fall into protected areas

## Data visualisation
- graph RQ TT
- PCA RQ ST
- 

