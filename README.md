<!-- README.md is generated from README.Rmd. Please edit that file -->



SpadeR (R package)
=====
<h4 style="text-align: right;">Most recent update time: May 12, 2016    

by Anne Chao, K. H., Ma, T. C., Hsieh and Chun-Huo Chiu.

Institute of Statistics, National Tsing Hua University, Hsin-Chu, Taiwan 30043</h4>


The program SpadeR (Species Prediction And Diversity Estimation) is an updated R package from the original version SPADE from Chao and Shen (2010).

Like [SPADE](http://chao.stat.nthu.edu.tw/wordpress/software_download/), the program SpadeR computes various biodiversity indices based on two types of sample data (abundance data and incidence data) from one or multiple communities. This user guide attempts to explain how to use this program in an easily accessible way using numerical examples and explanations.


### Software needed to run the development version

-   Required: [R](https://cran.r-project.org/)
-   Suggested: [RStudio IDE](http://www.rstudio.com/ide/download/)

### How to run
start R(Studio) and copy-and-paste the following commands:


```r
## install the latest version from github
install.packages('devtools')
library(devtools)
install_github('AnneChao/SpadeR')
library(SpadeR)
```

Remark that in order to install `devtools` package, you should update R
to the last version. Further, to get `install_github` to work, you
should install the `httr` package.

### The program is divided into six parts:

- Part I: Species (estimating species richness for one community based on five data formats options.


```r
# Example1: (abundance data)
data(ChaoSpeciesDataAbu)
ChaoSpecies(ChaoSpeciesDataAbu,"abundance",k=10,conf=0.95)

# Example2: (abundance frequency counts data)
data("ChaoSpeciesDataAbu_count")
ChaoSpecies(ChaoSpeciesDataAbu_count,"abundance_freq_count",k=10,conf=0.95)

# Example3: (incidence frequency data)
data(ChaoSpeciesDataInci)
ChaoSpecies(ChaoSpeciesDataInci,"incidence_freq",k=10,conf=0.95)

# Example4: (incidence frequency counts data)
data("ChaoSpeciesDataInci_freq_count")
ChaoSpecies(ChaoSpeciesDataInci_freq_count,"incidence_freq_count",k=10,conf=0.95)

# Example5: (incidence raw data)
data(ChaoSpeciesDataInci_raw)
ChaoSpecies(ChaoSpeciesDataInci_raw,"incidence_raw",k=10,conf=0.95)


```

- Part II: Diversity Profile Estimation (estimating a continuous diversity profile in one community including species richness, Shannon diversity and Simpson diversity) base on five dat formats options. This function also features plots of empirical and estimated continuous diversity profiles. Various estimates for Shannon entropy and the Gini-Simpson index are also provided. 

```r
# Example1: (abundance data)
data(DiversityDataAbu)
Diversity(DiversityDataAbu,"abundance",q=c(0,0.5,1,1.5,2))

# Example2: (abundance frequency counts data)
data("DiversityDataAbu_count")
Diversity(DiversityDataAbu_count,"abundance_freq_count",q=seq(0,3,by=0.5))

# Example3: (incidence frequency data)
data(DiversityDataInci)
Diversity(DiversityDataInci,"incidence_freq",q=NULL)

# Example4: (incidence frequency counts data)
data("DiversityDataInci_freq_count")
Diversity(DiversityDataInci_freq_count,"incidence_freq_count",q=NULL)

# Example5: (incidence raw data)
data(DiversityDataInci_raw)
Diversity(DiversityDataInci_raw,"incidence_raw",q=NULL)

```



- Part III: Shared Species (estimating the number of shared species between two communities) base on three dat formats options.

```r
# Example1: (abundance data)
data(ChaoSharedDataAbu)
ChaoShared(ChaoSharedDataAbu,"abundance",se=TRUE,nboot=200,conf=0.95)
# Example2: (incidence frequency data)
data(ChaoSharedDataInci)
ChaoShared(ChaoSharedDataInci,"incidence_freq",se=TRUE,nboot=200,conf=0.95)
# Example3: (incidence raw data)
data(ChaoSharedDataInci_raw)
ChaoShared(ChaoSharedDataInci_raw,"incidence_raw",units=c(16,17),se=TRUE,nboot=200,conf=0.95)

```

- Part IV: Two-Community Similarity Index (estimating various similarity indices for two assemblages) base on three dat formats options. The richness-based indices include the classic two-community Jaccard and Sorensen indices; the abundance-based indices include the Horn, Morisita-Horn, Two-community Bray-Curtis and the abundance-based Jaccard and Sorensen indices.


```r
# Example1: (abundance data)
data(SimilarityPairDataAbu)
SimilarityPair(SimilarityPairDataAbu,"abundance",nboot=200)
# Example2: (incidence frequency data)
data(SimilarityPairDataInci)
SimilarityPair(SimilarityPairDataInci,"incidence_freq",nboot=200)
# Example3: (incidence raw data)
data(SimilarityPairDataInci_raw)
SimilarityPair(SimilarityPairDataInci_raw,"incidence_raw",units=c(19,17),nboot=200)
```

- Part V: (estimating various N-community similarity indices) base on three dat formats options. The richness-based indices include the classic N-community Jaccard and Sorensen indices; the abundance-based indices include the Horn, Morisita-Horn, and the N-community Bray-Curtis indices.


```r
# Example1: (abundance data)
data("SimilarityMultDataAbu")
SimilarityMult(SimilarityMultDataAbu,"abundance",q=2,nboot=200,"relative")
# Example2: (incidence frequency data)
data("SimilarityMultDataInci")
SimilarityMult(SimilarityMultDataInci,"incidence_freq",q=2,nboot=200,"relative")
# Example3: (incidence raw data)
data("SimilarityMultDataInci_raw")
SimilarityMult(SimilarityMultDataInci_raw,"incidence_raw",units=c(19,17,15),q=2,nboot=200,"relative")
```

- Part VI: Genentic Measure (estimatinge allelic dissimilarity/differentiation among sub-populations based on multiple-population genetics data).


```r
# Example: (abundance data)
data(GeneticsDataAbu)
Genetics(GeneticsDataAbu,q=2,nboot=200)
```

### How to cite

If you publish your work based on results from `SpadeR`, please make references to our relevant papers mentioned in the following sections and also use the following reference for citing SpadeR:

> Anne Chao, K. H. Ma and T. C. Hsieh (2015). SpadeR: Species Prediction and Diversity Estimation with R. R package version 0.1.0. URL http://chao.stat.nthu.edu.tw/blog/software-download/

### License

The iNEXT package is licensed under the GPLv3. To help refine `SpadeR`, your comments or feedbacks would be welcome (please send them to [Anne Chao](chao@stat.nthu.edu.tw)).

