\name{ChaoSharedData}
\alias{ChaoSharedData}
\docType{data}
\title{
 Data for Function ChaoShared
}
\description{
 There are three data sets: \cr\cr
1. Type (1) abundance data (\code{ChaoSharedData$Abu}) \cr\cr
   The data consist of the observed bird abundances/frequencies collected from two estuaries (Chao et al. 2000). For each species (row), the entry of the first column is the observed species frequency from Estuary I, and the second column is the observed species frequency from Estuary II. The species checklist includes 201 species, so the entry data includes a matrix of 201 rows and 2 columns. \cr\cr

2. Type (2) incidence-frequency data (\code{ChaoSharedData$Inci}) \cr\cr
 The data consist of bird incidence (detection/non-detection) frequencies observed in 2015 (by 16 teams) and 2016 (by 17 teams) in the Hong Kong Bird Race. Each team is regarded as a sampling unit. Unlike the abundance data, the numbers of sampling units (16 and 17 for these data) are specified in the first row. Beginning with the second row, the entry of the first column is the observed incidence frequency (the total number of detections among all teams) of a given species in 2015, and the entry of the second column is the observed incidence frequency of the same species in 2016. A 280-species checklist was used, thus the input data consist of 281 rows (the first entry records the number of sampling units) and 2 columns. \cr\cr
3. Type (2B) incidence-raw data (\code{ChaoSharedData$Inci_raw}) \cr\cr
 The data consist of raw detection/non-detection records of bird species in 2015 (by 16 teams) and 2016 (by 17 teams) in the Hong Kong Bird Race. A 280-species checklist was used. The raw data consist of a 280 x 33 (species-by-sampling-unit) matrix with element 1's (detection) or 0's (non-detection). Each row of the matrix refers to the detection/non-detection records of the same species so that the information about shared species can be computed. The first 16 columns of the matrix denote the species detection/non-detection data by 16 teams in 2015, and the next 17 columns denote the species detection/non-detection data by 17 teams in 2016.
}
\usage{data(ChaoSharedData)}

\references{
 Chao, A., Hwang, W.-H., Chen, Y.-C. and Kuo. C.-Y. (2000). Estimating the number of shared species in two communities. Statistica Sinica, 10, 227-246. \cr\cr
 World Wildlife Fund (WWF) for Nature, Hong Kong. Bird Bird Race. \cr
 http://www.wwf.org.hk/en/getinvolved/hkbbr/. Assessed on July 26, 2016
}

\keyword{datasets}
