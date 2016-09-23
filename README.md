<h5 align="right">
Latest version: Sept. 9, 2016
</h5>
A Brief Introduction to SpadeR (R package): Species-Richness Prediction and Diversity Estimation
=====

<h4>Anne Chao, K. H., Ma, T. C., Hsieh and Chun-Huo Chiu.
<br><br>
Institute of Statistics, National Tsing Hua University, Hsin-Chu, Taiwan 30043</h4>

### Overview

SpadeR (Species-Richness Prediction and Diversity Estimation with R) is an updated R package from the original version of SPADE. SpadeR provides simple R functions to compute various biodiversity indices and related (dis)similarity measures based on individual-based (abundance) data or sampling-unit-based (incidence) data taken from one or multiple communities/assemblages. The SpadeR package is available in [CRAN](https://cran.r-project.org/web/packages/SpadeR/index.html). We have been updating SpadeR and you can download the latest version from Github (see below) or from [Anne Chao's website](http://chao.stat.nthu.edu.tw/wordpress/software_download/). 

<li> You need to acquire basic knowledge about R to use functions supplied by SpadeR.
<li> For readers without R background, please try our SpadeR Online, an R-based online version via the link [Anne Chao's website](http://chao.stat.nthu.edu.tw/wordpress/software_download/) or https://chao.shinyapps.io/SpadeR/. Users do not need to learn/understand R to run SpadeR Online.

Both SpadeR (R package) and SpadeR Online include nearly all of the important features from the original program SPADE while also having the advantages of expanded output displays and simplified data input formats. See [SpadeR Manual](https://cran.r-project.org/web/packages/SpadeR/SpadeR.pdf) for all details of the functions supplied in the package. For numerical examples with proper interpretations, see the detailed Online [SpadeR User's Guide](http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/SpadeR_UserGuide.pdf).

This package contains six main functions: <br>
1. <font color="red">ChaoSpecies</font> (estimating species richness for one community). <br><br>
2. <font color="red">Diversity</font> (estimating a continuous diversity profile and various diversity indices in one community including species richness, Shannon diversity and Simpson diversity). This function also features plots of empirical and estimated continuous diversity profiles. <br><br>
3. <font color="red">ChaoShared</font> (estimating the number of shared species between two communities). <br><br>
4. <font color="red">SimilartyPair</font> (estimating various similarity indices between two assemblages). Both richness and abundance-based two-community similarity indices are included. <br><br>
5. <font color="red">SimilarityMult</font> (estimating various similarity indices among N communities). Both richness and abundance-based N-community similarity indices are included. <br><br>
6. <font color="red">Genetics</font> (estimating allelic dissimilarity/differentiation among sub-populations based on multiple subpopulation genetics data). <br>

Except for the Genetics function, there are at least three types of data are supported for each function.

### Data Types

It is very important to prepare your data in correct format. Data are generally classified as abundance data and incidence data and there are five types of data input formats options (datatype="abundance", "abundance_freq_count", "incidence_freq", "incidence_freq_count", "incidence_raw"). 

<li> Individual-based abundance data when a sample of individuals is taken from each community. <br> 
<font color="red">Type (1) abundance data</font> (datatype = "abundance"): Input data consist of species (in rows) by community (in columns) matrix. The entries of each row are the observed abundances of a species in N communities. <br>
<font color="red">Type (1A) abundance-frequency counts data</font> only for a single community (datatype = "abundance_freq_count"): input data are arranged as (1 f1 2 f2 ... r fr)(each number needs to be separated by at least one blank space or separated by rows), where r denotes the maximum frequency and fk denotes the number of species represented by exactly k individuals/times in the sample. Here the data (f1, f2,..., fr) are referred to as "abundance-frequency counts".

<li> Sampling-unit-based incidence data when a number of sampling units are randomly taken from each community. Only the incidence (detection/non-detection) of species is recorded in each sampling unit. There are three data formats options. <br>

<font color="red">Type (2) incidence-frequency data</font> (datatype="incidence_freq"): The first row of the input data must be the number of sampling units in each community. Beginning with the second row, input data consist of species (in rows) by community (in columns) matrix. The entries of each row are the observed incidence frequencies (the number of detections or the number of sampling units in which a species are detected) of a species in N communities. <br>
<font color="red">Type (2A) incidence-frequency counts data</font> only for a single community (datatype="incidence_ freq_count"): input data are arranged as (T 1 Q1 2 Q2 ... r Qr) (each number needs to be separated by at least one blank space or separated by rows), where Qk denotes the number of species that were detected in exactly k sampling units, while r denotes the number of sampling units in which the most frequent species were found. The first entry must be the total number of sampling units, T. The data (Q1,Q2,...,Qr) are referred to as "incidence frequency counts". <br>
<font color="red">Type (2B) incidence-raw data</font> (datatype="incidence_raw"): Data consist of a species-by-sampling-unit incidence (detection/non-detection) matrix; typically "1" means a detection and "0" means a non-detection. Each row refers to the detection/non-detection record of a species in T sampling units. Users must specify the number of sampling units in the function argument "units". The first T1 columns of the input matrix denote species detection/non-detection data based on the T1 sampling units from Community 1, and the next T2 columns denote the detection/non-detection data based on the T2 sampling units from Community 2, and so on, and the last TN columns denote the detection/non-detection data based on TN sampling units from Community N, T1+ T2+ ... + TN = T.


### Software needed 

-   Required: [R](https://cran.r-project.org/)
-   Suggested: [RStudio IDE](http://www.rstudio.com/ide/download/)

### How to install
start R(Studio) and copy-and-paste the following commands:


```r
## install the latest version from github
install.packages('devtools')
library(devtools)
install_github('AnneChao/SpadeR')
library(SpadeR)
```

Remark that in order to install `devtools` package, you should update R
to the last version. Also, to get `install_github` to work, the `httr` package should be installed.

### Run SpadeR by examples

In the package, we have included many demo datasets for illustration. To gain familiarity with the program, we suggest that users first run the demo data sets included in SpadeR package and check the output with that given in the [SpadeR User's Guide](http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/SpadeR_UserGuide.pdf). Part of the output for each example is also interpreted in the guide to help users understand the statistical results. The formulas for estimators featured in SpadeR with relevant references are also provided in the SpadeR User's Guide. 

- Part I: ChaoSpecies (estimating species richness for one community).


```r
# Data for Function ChaoSpecies(data, datatype, k = 10, conf = 0.95)

data(ChaoSpeciesData)

# Type (1) abundance data
ChaoSpecies(ChaoSpeciesData$Abu,"abundance",k=10,conf=0.95)

# Type (1A) abundance frequency counts data
ChaoSpecies(ChaoSpeciesData$Abu_count,"abundance_freq_count",k=10,conf=0.95)

# Type (2) incidence frequency data
ChaoSpecies(ChaoSpeciesData$Inci,"incidence_freq",k=10,conf=0.95)

# Type (2A) incidence frequency counts data
ChaoSpecies(ChaoSpeciesData$Inci_count,"incidence_freq_count",k=10,conf=0.95)

# Type (2B) incidence raw data
ChaoSpecies(ChaoSpeciesData$Inci_raw,"incidence_raw",k=10,conf=0.95)
```

- Part II: Diversity (estimating a continuous diversity profile and various diversity indices in one community including species richness, Shannon diversity and Simpson diversity). This function also features plots of empirical and estimated continuous diversity profiles. 

```r
# Data for Function Diversity(data, datatype, q = NULL)

data(DiversityData)

# Type (1) abundance data
Diversity(DiversityData$Abu,"abundance",q=c(0,0.5,1,1.5,2))

# Type (1A) abundance frequency counts data
Diversity(DiversityData$Abu_count,"abundance_freq_count",q=seq(0,3,by=0.5))

# Type (2) incidence frequency data
Diversity(DiversityData$Inci,"incidence_freq",q=NULL)

# Type (2A) incidence frequency counts data
Diversity(DiversityData$Inci_freq_count,"incidence_freq_count",q=NULL)

# Type (2B) incidence raw data
Diversity(DiversityData$Inci_raw,"incidence_raw",q=NULL)
```



- Part III: ChaoShared (estimating the number of shared species between two communities).

```r
# Data for Function ChaoShared(data, datatype, units, se = TRUE, nboot = 200, conf = 0.95)

data(ChaoSharedData)

# Type (1) abundance data
ChaoShared(ChaoSharedData$Abu,"abundance",se=TRUE,nboot=200,conf=0.95)

# Type (2) incidence frequency data
ChaoShared(ChaoSharedData$Inci,"incidence_freq",se=TRUE,nboot=200,conf=0.95)

# Type (2B) incidence raw data
ChaoShared(ChaoSharedData$Inci_raw,"incidence_raw",units=c(16,17),se=TRUE,nboot=200,conf=0.95)

```

- Part IV: SimilartyPair (estimating various similarity indices between two assemblages). Both richness and abundance-based two-community similarity indices are included.


```r
# Data for Function SimilarityPair(data, datatype, units, nboot = 200)

data(SimilarityPairData)

# Type (1) abundance data
SimilarityPair(SimilarityPairData$Abu,"abundance",nboot=200)

# Type (2) incidence frequency data
SimilarityPair(SimilarityPairData$Inci,"incidence_freq",nboot=200)

# Type (2B) incidence raw data
SimilarityPair(SimilarityPairData$Inci_raw,"incidence_raw",units=c(19,17),nboot=200)
```

- Part V: SimilarityMult (estimating various similarity indices among N communities). Both richness and abundance-based N-community similarity indices are included.


```r
# Data for Function SimilarityMult(data, datatype, units, q, nboot = 200, goal)

data(SimilarityMultData)

# Type (1) abundance data
SimilarityMult(SimilarityMultData$Abu,"abundance",q=2,nboot=200,"relative")

# Type (2) incidence frequency data
SimilarityMult(SimilarityMultData$Inci,"incidence_freq",q=2,nboot=200,"relative")

# Type (2B) incidence raw data
SimilarityMult(SimilarityMultData$Inci_raw,"incidence_raw",
units=c(19,17,15),q=2,nboot=200,"relative")
```

- Part VI: Genetics (estimating allelic dissimilarity/differentiation among sub-populations based on multiple subpopulation genetics data).


```r
# Data for Function Genetics(X, q, nboot = 200)
data(GeneticsDataAbu)
# Type (1) abundance data 
Genetics(GeneticsDataAbu,q=2,nboot=200)
```

### How to cite

If you publish your work based on results from `SpadeR`, please make references to our relevant methodology papers mentioned below and also use the following reference for citing SpadeR:

Chao, A., Ma, K. H., Hsieh, T. C. and Chiu, C. H. (2016). SpadeR (Species-richness Prediction And Diversity Estimation in R): an R package in CRAN. Program and User’s Guide also published at http://chao.stat.nthu.edu.tw/blog/software-download/

### References

We recommend the following recent papers for pertinent background on biodiversity measures and statistical analyses. These papers can be directly downloaded from Anne Chao’s website. 

Chao, A., and Chiu, C. H. (2012). Estimation of species richness and shared species richness. In N. Balakrishnan (ed). Methods and Applications of Statistics in the Atmospheric and Earth Sciences. p.76–111, Wiley, New York. (<font color="blue">Background on species richness and shared species richness estimation</font>)

Chao, A., and Chiu, C. H. (2016). Nonparametric estimation and comparison of species richness. Wiley Online Reference in the Life Science. In: eLS. John Wiley & Sons, Ltd: Chichester.  (<font color="blue">Background on comparing species richness across communities</font>)

Chao, A., and Chiu, C. H. (2016). Bridging the variance and diversity decomposition approaches to beta diversity via similarity and differentiation measures. Methods in Ecology and Evolution, 7, 919–928. (<font color="blue">A unified theoretical framework on similarity/differentiation measures</font>)

Chao, A., Chiu, C. H. and Jost, L. (2014). Unifying species diversity, phylogenetic diversity, functional diversity, and related similarity/differentiation measures through Hill numbers. Annual Reviews of Ecology, Evolution, and Systematics, 45, 297–324. (<font color="blue">A unified theoretical framework on diversity measures</font>)

Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K. and Ellison, A. M. (2014). Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. Ecological Monographs, 84, 45–67. (<font color="blue">Background on comparing diversity measures across communities</font>)

Chao, A. and Jost, L. (2015). Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution, 6, 873–882. (<font color="blue">A unified approach to estimating diversity in a community based on incomplete samples</font>)

Chao, A., Wang, Y. T. and Jost, L. (2013). Entropy and the species accumulation curve: a novel entropy estimator via discovery rates of new species. Methods in Ecology and Evolution, 4, 1091–1100. (<font color="blue">A nearly optimal estimator of Shannon entropy/diversity based on incomplete samples</font>)

