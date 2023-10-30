# SeedCollectionSize
Tools for estimating the mass required to achieve an operational Seeds of Success collection.

Recent estimates indicate that the Federal land management agencies of the United States require roughly 1.9 trillion seeds for ecological restoration.
In order to supply this demand, seed collected from wild lands are 'increased' in agricultural settings, a process where wild-collected seed are raised to maturity and have there seed collected for up to several generations.
Large, 'Operational', collections of wild land seed are required to produce a standing crop which is economically feasible for farmers to raise to maturity and maintain.
Operational collections, requiring upwards of 80,000 seeds, can be difficult and timely to collect.
Further, many species are producing and dispersing, mature seed at similar times of year. 
This overlap in the timing of seed dispersal makes collecting operational collections from many species simultaneously challenging. 
The diligent allocation of human resources are required to obtain multiple operational collections per season.

This project aims to develop robust methods to estimate the mass of seeds required to achieve an operational collection. Two formats of developing these estimates exist:

1) Bootstrap re-sampling of raw, dried, collections which have subsequently had their total seeds counted, can produce summary statistics to guide in-field efforts. 
2) Generalized linear models, via poisson regression, can incorporate terms which reflect the environmental variability associated with the known collections. 

The latter models are used by contractors to more accurately determine the number of seeds collected for billing purposes.

## Data Sources

The Seeds of Success (SOS) database contains both known 1) raw (dried) collection weights, and 2) percent seed viability.   
These data were largely collected across Bureau of Land Management land in 13+ states, and for many species collections are representative of their ecological gradations in the Western United States. 
The SOS program contains seeds collected in two relevant modes of collection, 'conservation' collections target 10,000 seed collections in order to preserve the germplasm of a population, while the aforementioned 'operational' collections both preserve the germplasm and provide farmers adequate stocks to cultivate the lineage. 
Seeds of Success crews have been making conservation collections since 2002, given both their smaller size and the longer implementation time relative to operational collections - which began around 2015, the number of conservation collections exceeds the number of 'operational' collections, but effectively follow the same process. 
Both 'operational' and 'conservation' collections were considered simultaneously for our analyses in order to increase sample sizes. 
The dried collection weights served as the reponse variable, while seed viability served as a covariate. 

Spatial predictors include:   

| Layer |                       Description                       |              Source                   | 
| :---: | :-----------------------------------------------------: | :-----------------------------------: |
|   1   |                   soil % loam  (0-30 cm)                |              SoilGrids                |
|   2   |                     Aridity Index                       |               Chelsa                  |
|   3   |             Mean annual precipitation (BIO12)           |                PRISM                  |
|   4   |         Precipitation of Warmest Quarter (BIO18)        |              WorldClim                |
|   5   |           Number of growing degree days above 5*c       |               Chelsa                  |
|   6   |  Standardized Precip. Evapotranspiration Index   6 (mo) |                SPEI                   |
|   7   |  Standardized Precip. Evapotranspiration Index  12 (mo) |                SPEI                   |
|   8   |  Standardized Precip. Evapotranspiration Index  24 (mo) |                SPEI                   |
|   9   |                      Latitude                           |                 NA                    |
|  10   |                     Longitude                           |                 NA                    |

Wherein only one SPEI value may be retained in the model.   

## Models

### Non-parametric summary statistics
For both this, and linear models, the number of seeds per 10 grams will be used to avoid stochasticity in rounding to a whole number at smaller masses, e.g. 3.5 seeds per one gram could be rounded to 3 or 4, vastly increasing uncertainty in estimates, relative to 35 seeds per 10 grams. 
The number of seeds per unit mass will be non-parametrically re-sampled to develop estimates of the mean number of seeds per unit weight.

### Linear Models

#### Final Models

All initial stages of modelling will be performed using MuMin::pdredge. 
`pdredge`, is simply a form of dredge which will allow for parallel computations if supplied a cluster. 

Seeds/Unit Mass ~ with(% viability, spatial weights) * aridity index * MAP * MPWQ * !(SPEI6 && SPEI 12 && SPEI24),  m.lim = c(0, 4), family = poisson)

Percent viable seed are maintained in all models, this is an additional metric which is present both in the processed Bend data, and which is collected by crews on the date of collection. 
As I suspect the largest proportion of mass per, dried, collection is associated with seeds, whether viable or not. 
Only a single SPEI term will be retained in the final model. 
Models are then automatically selected via AIC tables.
Finally all models with <sub>d</sub>AIC < 2 are ensembled. 

#### Spatial auto-correlation

The effect of spatial auto-correlation is expected to be moderately strong for most species, which may result in improper estimates of the dispersion of residuals around our fit model. 
Subsequent to selecting or ensembling a top model, a GLS with the terms and interactions specified will be generated. 
This final model will be constructed with the package ... using the function ... which utilizes the eigenvectors across each site included in the training data. 

## Prediction

All final models will be saved as R Data Objects. 
Prediction are made at time of model generation using the generic 'predict' function, with 1000 rows of data which capture the feature space used to train models. 
A look up function is used to refer to these saved tables and find the predicted value most similar to the observed  parameters...  
The only variables which crews need to collect are the mass of collections, and the percent viability, metrics which they already use in the existing methods of sample size estimation. 
The spatial data associated with collections, and which are required to predict collection sizes, and gleaned from in-field electronic data collection. 

