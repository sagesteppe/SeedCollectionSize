# SeedCollectionSize
Tools for estimating the mass required to achieve an operational Seeds of Success collection.

Operational collections of wildland seed are required to generate a wild-type cultivar which may be reared in agricultural settings to produce seed for wildlands restoration.
Operational collections require either 80,000 or 90,000 seeds, and are time intensive to collect.
Further many species are producing, and dispersing, mature seed at similar times of year. 
Allocating human resources diligently is required to obtain multiple operational collections per season.

This project aims to develop robust methods to estimate the mass of seeds required to achieve an operational collection. Two formats of developing these estimates exist:

1) Bootstrap re-sampling of raw, dried, collections which have subsequently had their total seeds counted, can produce summary statistics to guide in-field efforts. 
2) Generalized linear models, via poisson regression, can incorporate terms which reflect the environmental variability associated with the known collections. 

## Data Sources

The Seeds of Success database contains known: Raw (dried) collection weights, percent seed viability, and number of seeds per unit volume.   

Spatial predictors include:   

| Layer |                       Description                       |              Source                   | 
| :---: | :-----------------------------------------------------: | :-----------------------------------: |
|   1   |                   soil % loam  (0-30 cm)                |              SoilGrids                |
|   2   |                       Aridity Index                     | Global Ariidty Index and PET database |
|   3   |             Mean annual precipitation (BIO12)           |        CLIMATENA/PRISM/DISMO          |
|   4   |         Precipitation of Warmest Quarter (BIO18)        |        PRISM / CLIMATENA / DISMO      |
|   5   |  Standardized Precip. Evapotranspiration Index   6 (mo) |                SPEI                   |
|   6   |  Standardized Precip. Evapotranspiration Index  12 (mo) |                SPEI                   |
|   7   |  Standardized Precip. Evapotranspiration Index  24 (mo) |                SPEI                   |

Wherein only one SPEI value may be retained in the model.   

## Models

### Non-parametric summary statistics
For both this, and linear models, the number of seeds per 10 grams will be used to avoid stochasticity in rounding to a whole number at smaller masses, e.g. 3.5 seeds per one gram could be rounded to 3 or 4, vastly increasing uncertainty in estimates. 
The number of seeds per unit mass will be non-parametrically re-sampled to develop estimates of the mean number of seeds per unit weight.


### Linear Models
All modelling will be performed using MuMin::dredge.

Seeds/Unit Mass ~ % viability * aridity index * MAP * MPWQ * SPEI6 * SPEI12 * SPEI24,  family = poisson)

Only a single SPEI term will be retained in the final model. 
Models are then automatically selected via AIC tables.
Finally all models with <sub>d</sub>AIC < 2 are ensembled. 

## Prediction

All final models will be saved as R Data Objects. Prediction will occur throughout the field season. 
The only variables which crews need to collect are the mass of collections, and the percent viability, metrics which they already use in the existing methods of sample size estimation. 
