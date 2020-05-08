# MCC-HetPoll

Here is a summary of what I have done so far.

## Data

For the first stage PM2.5 data, I roughly follow previous approach to clean the data. I exclude cities that have less than two years of pollution data. I also exclude the cities that have a too high proportion of missing data between the first and last valid data.

For the constituent data (SPEC), I hereby use the 10km buffer aggregation since it allows to have a lower proportion of zero values. But overall the results are really similar between the two datasets.

To harmonize both datasets (pollution and SPEC), I only selected common years between them. This corresponds to years between 2003 and 2017 for which there are PM2.5 values. It significantly reduces the time series length for some countries especially for North American ones (USA and Canada).

![TSlength](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/1_TSlength.png)

Here, I compare the sum of mean concentration of each constituent ofr each country to the mean of observed PM2.5 for each year. It shows that the constituent sum does not quite add up to the observed level of PM2.5.

![SpecCountries](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/1_SpecCountries.png)

## First stage model

For the first stage model, I compared the QAIC and RR of the city level models as published in the NEJM paper to alternative parametrizations. The only change the seems to dimish QAICs is more controlling for temperature, i.e. with a lag of 21 days and splines for the lag dimension (bottom middle panel). 

![QAIC](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/2bis_QAIC.png)

Also, when comparing with the models that account for the maximum data available (i.e. also prior to 2003), RRs are slightly more variables, as could be expected (bottom right panel).

![RR10](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/2bis_RR10.png)

Thus, I finally apply the same parametrization as in NEJM paper for PM2.5, but with a more aggressaive confounding for temperature (DLNM with lags up to 21) and no confounding of humidity since the variable is not always available.

## Second stage (in progress)

Below, I show the BLUP by applying the meta-regression with both city and country as random effects. As expected, fewer countries display significant RR than in the NEJM paper since less data are used (n.b. 2003 and after).  

![BLUP](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3a_forestplot_country.png)

In the following, I integrate the overall mean of each constituent as fixed-effect meta-predictors. According to this result, only ammonia (NH4) seems to have an aggravating effect.

![betaSPEC](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3b_forestplot_betaSPEC.png)
