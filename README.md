# Differential mortality risks associated to PM2.5 components: a multi-country multi-city study.

Reproducible R code and results for the publication 

Masselot P, et al. Differential mortality risks associated to PM2.5 components: a multi-country multi-city study. *Epidemiology*. In press.

### Data and Results

Figures and Tables (including supplemental ones) are available in the Results folder of this repository. Data are not available currently due to restricted data sharing agreement between the collaborators of this study.

### R code

The R code to reproduce the analysis is available. Scripts are meant to be executed in order:

- *0_PrepData.R* Loads and and links mortality, PM2.5, composition and city-specific characteristics. Performs city selection.
- *1a_DataSummary.R* Creates Table 1 providing descriptive statistics.
- *1b_CompositionSummary.R* Provides descriptive statistics of PM2.5 compositions. Produces Figure 2.
- *2_FirstStage.R* Runs the first-stage model on each selected city to produce city-level relative risks.
- *3_SecondStage_MetaRegComposition.R* Pools first-stage results in a meta analysis including city-specific compositions and socio-economic indicators. Produces Figures 1 and 3 as well as several eFigures.
- *4_SupplementaryResults.R* Produces several supplementary results and eFigures, including the PCA summary and residuals analysis.
