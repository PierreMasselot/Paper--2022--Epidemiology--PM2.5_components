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

The following plot shows the variation matrix of constituents proportion (averaged by city). This matrix computes the variance of the log-ratios between each constituent. When this variance is high, this means that the two constituents tend to vary in opposite directions (when one increases, it lowers down the proportion of the other). DUST and SS seem to vary a lot compared to other constituents, especially compared to NIT for DUST and to BC for SS. This may be explained by their high proportion of zero values, meaning that when one of them is present, all other constituents decrease. In addition, NIT and BC also seem to vary in opposite directions.

![VariationMat](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/1b_VariationMatrix.png)

It can also be interesting to do a PCA on proportions. The first component seem to mainly contrast SS and DUST to OC and BC. The second contrasts DUST and BC to NH4, NIT and SS.

![PCA](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/1b_PCAbiplot.png)

## First stage model

For the first stage model, I compared the QAIC and RR of the city level models as published in the NEJM paper to alternative parametrizations. The only change the seems to dimish QAICs is more controlling for temperature, i.e. with a lag of 21 days and splines for the lag dimension (bottom middle panel). 

![QAIC](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/2bis_QAIC.png)

Also, when comparing with the models that account for the maximum data available (i.e. also prior to 2003), RRs are slightly more variables, as could be expected (bottom right panel).

![RR10](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/2bis_RR10.png)

I finally apply the same parametrization as in NEJM paper for PM2.5, but with no confounding of humidity since the variable is not always available. The confounding of temperature is slightly modified, replacing the 6 degrees of freedom of ns by knots at quantiles 10, 75 and 90 %.

## Second stage

For the second stage, I use `mixmeta` with random effects for city and country. To account for confounders, I selected the MCC indicators that induced the largest effect modifications in the preliminary analysis of Francesco. Five are from the OECD dataset (Proportion of pop > 65 years, GDP, Poverty index, Average temperature and temperature range) and Four from the UCD dataset (Greenness in 2000 and 2014, Total built-up area in 2000 and 2015).

To measure the effect modification of each component proportions, they cannot be used directly as meta-predictors, since they sum to one. This results in a singular matrix and can also lead to spurious correlations as shown by Aitchison in his work. Thus, I use logratios. They express variations of some components relatively to others. Mathematically, this maps the ratios to a classical euclidean geometry. 
 
In a regression context, there are two possibles logratios family: i) isometric logratios (ILR), also called balances when several components are aggregated, and ii) summated logratios (SLR), also called amalgamations in the compositional data literature. To express the variation of a component vs others in the ILR case, it is expressed as the logratio between the component and $geometrical mean$ of the others. In the SLR case, it is expressed as the logratio between the component and the $sum$ of the others. When aggregating several components (such as secondary inorganic aerosols vs organic components), the ILR contrasts them as the logratios between their $geometrical means$ while the SLR simply contrasts them as the logratio between their $sums$.

In term of interpretation, the main difference is that the geometrical mean of a group is affected by the subcomposition inside the group, while the sum is not. Thus, in the former (ILR), the variations between each components of the secondary inorganic aerosols will change its logratio contrasting it to organic components, while it doesn't matter with SLR. To me, the SLR makes more sense for interpretation (and it is well shown in Greenacre 2019), but ILR is more often championed by compositional data analysts, perhaps because it is more 'elegant' mathematically (Egozcue & Pawlowski-Glahn 2005, Fiserova & Hron 2011). Note that the ILR approach was the first one I applied by default because it is more represented in the literature.

In the following, I apply both of them to see the results in three cases as discussed previously. The first one is to estimate the effect modification of each single component versus all others. The second one is to estimate the effect modification of several groups of components, as discussed between us and done in Hvidtfeldt et al. (2019), which are secondary inorganic aerosols (SO4, NH4 and NIT), organic components (BC and OM) and considering DUST and SS separately. The last one is the effect modification of traffic related components (SO4, NIT and BC). Fortunately, there are some agreement between the results of both approaches.
 
### ILR approach

This approach is implemented in the script `3d_SecondStage_balances.R`.
NB: The results here are different than what I showed before because of a mistake of mine in using the R package `compositions`. I now understand it better and am more satisfied with the results. 
The first graph below shows effect modification of each individual component and the second of each group defined above.

![3d_forestplot_all](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3d_forestplot_all.png) 
![3d_forestplot_balances](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3d_forestplot_balances.png)

Here, the single component that seems to have the largest effect modification is Ammonia (NH4), seemingly associated with livestock and fertilizer use. The grouped results indicate effect modification of organic component, while surprisingly none of its components (BC and OC) result in effect modification. Finally, I found no effect of traffic-related components.

### SLR approach

This approach is implemented in the script `3e_SecondStage_amalgamation.R`. 

![3e_forestplot_All](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3e_forestplot_All.png) 
![3e_forestplot_aggregated](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3e_forestplot_aggregated.png)

Here, we still find NH4 having effect modification, but this time there is also BC and, to a lesser degree, OC. The grouped results overall agree with the ILR approach although falling short of being significant at the 5% level. This may be because of NH4. Finally, no effect of traffic-related components where found here neither.

### PCA to compare

Below is the result of regression with the two PCA component described above. The 1st component that contrasts organic components and NH4 to SS and DUST also results in effect modification. We have a nice consistency here.

![3c_forestplot_pca](https://github.com/PierreMasselot/MCC-HetPoll/blob/master/Results/3c_forestplot_pca.png)

## Conclusion

The results suggest effect modification of NH4 first and then organic components. To date, I havn't found a paper reporting an effect of NH4 specifically, although some report effects of secondary inorganic aerosols. Effect of organic components however would be the most consensual result.
In the current version of the paper, I chose to report the ILR approach but presented in a different way. I realized that it is equivalent to the first approach proposed by Aitchison and Bacon-Shone (1984) and is thus backed by literature. Although the SLR approach makes sense to me, it is not backed in this context.
## References cited here

Aitchison, J., Bacon-Shone, J., 1984. Log contrast models for experiments with mixtures. Biometrika 71, 323–330. https://doi.org/10.1093/biomet/71.2.323

Egozcue, J.J., Pawlowsky-Glahn, V., 2005. Groups of Parts and Their Balances in Compositional Data Analysis. Math Geol 37, 795–828. https://doi.org/10.1007/s11004-005-7381-9

Fišerová, E., Hron, K., 2011. On the Interpretation of Orthonormal Coordinates for Compositional Data. Math Geosci 43, 455. https://doi.org/10.1007/s11004-011-9333-x

Greenacre, M., Grunsky, E.C., Bacon-Shone, J., 2019. A comparison of amalgamation and isometric logratios in compositional data analysis. Preprint.

Hron, K., Filzmoser, P., Thompson, K., 2012. Linear regression with compositional explanatory variables. Journal of Applied Statistics 39, 1115–1128. https://doi.org/10.1080/02664763.2011.644268

Hvidtfeldt, U.A., Geels, C., Sørensen, M., Ketzel, M., Khan, J., Tjønneland, A., Christensen, J.H., Brandt, J., Raaschou-Nielsen, O., 2019. Long-term residential exposure to PM2.5 constituents and mortality in a Danish cohort. Environment International 133, 105268. https://doi.org/10.1016/j.envint.2019.105268