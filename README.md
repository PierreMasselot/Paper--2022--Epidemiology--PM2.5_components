# Differential mortality risks associated to PM2.5 components: a multi-country multi-city study.

R code and results attached to the publication:

Masselot P, et al. Differential mortality risks associated to PM2.5 components: a multi-country multi-city study. *Epidemiology*. In press.

### Data and Results

**Data are not available currently due to restricted data sharing agreement between the collaborators of this study. Therefore, the code is not fully reproducible.**

Data are normally included in a subfolder *Data*. It should contain:
- Mortality and pollution data, stored in a list of city-specific data.frames. Also contains a descriptive data.frame with one line for each city.
- PM2.5 components, stored in one csv files per year of data. Each csv file contains one line per city.
- City-specific characteristics, stored as a data.frame with one line per city.

### R code

The R code to reproduce the analysis is available. Scripts are meant to be executed in order:

- *0_PrepData.R* Loads and and links mortality, PM2.5, composition and city-specific characteristics. Performs city selection.
- *1a_DataSummary.R* Creates Table 1 providing descriptive statistics.
- *1b_CompositionSummary.R* Provides descriptive statistics of PM2.5 compositions. Produces Figure 2.
- *2_FirstStage.R* Runs the first-stage model on each selected city to produce city-level relative risks.
- *3_SecondStage_MetaRegComposition.R* Pools first-stage results in a meta analysis including city-specific compositions and socio-economic indicators. Produces Figures 1 and 3 as well as several eFigures.
- *4_SupplementaryResults.R* Produces several supplementary results and eFigures, including the PCA summary and residuals analysis.

### Results

Figures and Tables generated for the publication (including supplemental ones) are available in the Results folder of this repository. 

### Acknowledgements

This work was supported by the Medical Research Council of UK (Grant ID: MR/M022625/1), the Natural Environment Research Council of UK (Grant ID: NE/R009384/1), and the European Union’s Horizon 2020 Project Exhaustion (Grant ID: 820655).
