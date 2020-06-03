# MCC-HetPoll

Here is a summary of what I have done so far.

## Data

For the first stage PM2.5 data, I roughly follow previous approach to clean the data. I exclude cities that have less than one year of pollution data. I also exclude the cities that have a too high proportion of missing data between the first and last valid data.

For the constituent data (SPEC), I hereby use the 10km buffer aggregation since it allows to have a lower proportion of zero values. But overall the results are really similar between the two datasets. 

As discussed, I make the hypothesis that pollution is time-invariant. This means that I use all PM2.5 data after 1999 for the first stage, and all constituent data (2003-2017).

Here, I compare the sum of mean concentration of each constituent ofr each country to the mean of observed PM2.5 for each common year. It shows that the constituent sum does not quite add up to the observed level of PM2.5. 

![SpecCountries](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/1a_SpecCountries.png)

Below is the same plot but showing proportions for all years. Overall, there is little variation of proportions within countries from year to years, athough there are few trends. For instance the proportion of DUST seems to increase in Greece while decreasing on Portugal. 

![PropCountries](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/1b_CompCountries.png)

The following plot shows the variation matrix of constituents proportion (averaged by city). This matrix computes the variance of the log-ratios between each constituent. When this variance is high, this means that the two constituents tend to vary in opposite directions (when one increases, it lowers down the proportion of the other). DUST and SS seem to vary a lot compared to other constituents. This may be explained by their high proportion of zero values, meaning that when one of them is present, all other constituents decrease. 

![VariationMat](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/1b_VariationMatrix.png)

It can also be interesting to do a PCA on proportions. Mainly the two "organic" components (BC and OC) have the same direction, so have NIT and NH4 which are two of the three secondary inorganix aerosols.

![PCA](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/1b_PCAbiplot.png)

## First stage model

For the first stage model, I compared the QAIC and RR of the city level models as published in the NEJM paper to alternative parametrizations. The only change the seems to dimish QAICs is more controlling for temperature, i.e. with a lag of 21 days and splines for the lag dimension (bottom middle panel). 

![QAIC](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/2bis_QAIC.png)

Also, when comparing with the models that account for the maximum data available (i.e. also prior to 2003), RRs are slightly more variables, as could be expected (bottom right panel).

![RR10](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/2bis_RR10.png)

I finally apply the same parametrization as in NEJM paper for PM2.5, but with no confounding of humidity since the variable is not always available. The confounding of temperature is slightly modified, replacing the 6 degrees of freedom of ns by knots at quantiles 10, 75 and 90 %.

## Second stage (in progress)

For the second stage, the constituents proportions cannot be used directly as predictors, since they sum to one. This results in a singular matrix and can also lead to spurious correlations as shown by Aitchison in his work.
I apply the procedure of Martín-Fernández et al. (2003), which goes through the isometric logratio transform. Basically, each constituent is replaced by the logarithm of its value divided by the geometric mean of other constituents. The resulting variables thus represent the variation of the constituent relatively to others. Mathematically, this maps the ratios to a classical euclidean geometry.

Below are the coefficients obtained using mixmeta with the logratios as meta-predictors. With this approach, no constituent seems to signifcantly modify the RR compared to all other.

![metaComp](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3c_forestplot_logratio.png)

The advantage of compositional data analysis, is that we can analyze subcompositions, i.e. to check the variation of several constituents compared to others. Below show effect modifications estimated when considering only secondary inorganic components (INOR, which are SO4, NH4 and NIT). This means that we estimate the impact of each compared only to other secondary inorganic components independently to the four other components. In this case, a higher proportion of SO4 seems to lead to higher RRs and conversely, a higher proportion of NH4 to lower RRs. 

The fourth value is the effect of high proportion of secondary inorganic components compared to other components (BC, OC, SS, DUST). And it seems that these components jointly lead to higher RRs.

![metaInor](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3c_forestplot_secondaryInorganic.png)

Next, we consider an aggregated composition in which components are gathered in three groups: secondary inorganic pollutants (INOR with SO4, NIT and NH4), organic pollutants (BC and OM) and natural particulate matter (SS and DUST). We want to estimate the effect modification of an increase of each of the two formers compared to the natural group as a baseline. Thus the meta-regression is performed here through the addtive log-ratio transform (ALR) of Aitchison, that divide each component by the last one (natural then) and compute the log. The forestplot is shown below and doesn't suggest any significant effect.

![metaAgg](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3c_forestplot_aggregated.png)