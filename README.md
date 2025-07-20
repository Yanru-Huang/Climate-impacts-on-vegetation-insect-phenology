## Project Overview

This study evaluates how climate change has affected the phenological synchrony between vegetation and insects across Europe. Leveraging 34 years of remote‑sensing observations and citizen‑science records, we analysed observations of 1,495 phytophagous insect species (spanning four orders) together with the vegetation start‑of‑season (SOS) at the same spatio‑temporal points. Our goal was to quantify long-term trends in vegetation and insect phenology and analyse potential drivers.

---

## File Descriptions

### `01_Classification_of_Insect_Phenological_Patterns`

Classifies insect occurrence records by the number and timing of activity peaks.

* **Removing spatial & year effects**:Fit a linear mixed effects model for each species, retaining the residuals unaffected by latitude, longitude, and year for the next modeling step.
* **Peak detection**: Smooth residuals with a Gaussian kernel density estimate and Savitzky–Golay filter. Identify potential peaks (*P*) using `findpeaks()` from the pracma package.
* **Pattern classification**: For species with *P* > 1, run K‑means clustering on observation dates with `K = 2 … P`. Compute the Silhouette Score for each *K*; choose the *K* with the highest score as the optimal number of peaks (clusters).

---

### `02_Spatiotemporal_Differences_in_the_Analysis_of_Vegetation–Insect_Phenology`

Builds a linear mixed‑effects model to quantify 1982–2015 change rates in vegetation and insect phenology at both the continental scale and the 1° × 1° regional scale.
**Data pre‑filtering**: Retain datasets that meet all of the following:  ≥ 400 total records. ≥ 30 records and coverage of ≥ 10 years per grid cell.

---

### `03_Analysis_of_the_Impact_of_Environmental_Factors_on_Vegetation_and_Insect_Phenologies`

Assesses how environmental variables influence vegetation and insect phenology using partial correlation, ridge regression, and random forests.

#### Environmental variables considered

| Variable                              | Field name*                    |
|---------------------------------------|--------------------------------|
| Total evapotranspiration              | `total_evaporation_sum`        |
| Mean relative humidity                | `hu`                           |
| Precipitation                         | `rr`                           |
| Surface radiation                     | `qq`                           |
| Mean air temperature                  | `tg`                           |
| Mean wind speed                       | `fg`                           |
| Volumetric soil‑water content         | `volumetric_soil_water_layer_1`|
| Soil temperature (0 – 7 cm)           | `soil_temperature_level_1`     |

\*Field names match those in the raw data files.

**Data matching**: Extract environmental variables for the 0 – 6 months preceding each phenological date. Average values across that window.
 **Optimal pre‑season selection**: Detrend phenological dates for inter‑annual variation. Calculate correlations between phenology and each environmental variable across all candidate windows. Select the window with the highest absolute correlation as the optimal pre‑season.
**Statistical analyses**: Standardise variables from the optimal pre‑season. Run partial correlation, ridge regression, and random‑forest models to quantify each variable’s contribution.
