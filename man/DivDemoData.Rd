\name{DiversityData}
\alias{DiversityData}
\docType{data}
\title{
  Data for Function Diversity
}
\description{
 There are five data sets: \cr\cr
 1. Type (1) abundance data (\code{DiversityData$Abu}) \cr\cr
 The data include a column of the observed tree abundances/frequencies from an old-growth rain forest in Costa Rica (Chao et al. 2005, 2008). There were 69 tree species among 557 individuals. \cr\cr
 2. Type (1A) abundance-frequency counts data (\code{DiversityData$Abu_count}) \cr\cr
 The data consist of the observed beetles species abundance-frequency counts collected from the Osa old-growth forest site in Costa Rica (Janzen, 1973). There were 112 species among 237 individuals. The input abundance-frequency counts data are arranged as = (1 84 2 10 3 4 4 3 ... 42 1); each number needs to be separated by at least one blank space or separated by rows. Here the first pair (1, 84) indicates that there are 84 singletons, the second pair (2, 10) indicates there are 10 doubletons, and so on, with the last pair (42, 1) indicating that there is one species that is represented by 42 individuals. \cr\cr
 3. Type (2) incidence-frequency data (\code{DiversityData$Inci}) \cr\cr
 The single-column data include the observed incidence-based frequencies of tropical rainforest ants collected by Berlese extraction of soil samples (217 sampling units) in Costa Rica (Longino et al. 2002). In the input data, the entry in the first row denotes the number of sampling units (217); the subsequent 117 rows denote species incidence frequencies of the observed species.  \cr\cr
 4. Type (2A) incidence-frequency counts data (\code{DiversityData$Inci_freq_count}) \cr\cr
 The seed-bank data consist of the observed species incidence-based frequency counts of seedlings that germinated from soil samples (Butler and Chazdon, 1998); here each soil sample is regarded as a sampling unit. A total of 34 species of seedlings were found in the 121 soil samples. The incidence frequency counts are read as = (121 1 3 2 2 3 3 ... 61 1); each number needs to be separated by at least one blank space or by separated by rows. The first entry, indicating that there are 121 soil samples, is followed by the 18 pairs (1, 3), (2, 2), (3, 3), (4, 3), (5, 1), (6, 5), and so on, up to (61, 1). Here (1, 3) indicates that there are 3 unique species, (2, 2) indicates there are 2 duplicate species, and so on, with 
 (61, 1) indicating that there is one species found in 61 soil samples. \cr\cr
 5. Type (2B) incidence-raw data (\code{DiversityData$Inci_raw}) \cr\cr
 The data consist of raw incidence data of the seed-bank records, described above for the incidence frequency counts data. The input data include a 34 x 121 (species-by-sampling-unit) matrix. For each element of the matrix, "1" means a detection and "0" means a non-detection.
 
 
}
\usage{data(DiversityData)}

\references{
  Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters, 8, 148-159. \cr\cr
  Chao, A., Jost, L., Chiang, S.-C., Jiang, Y.-H. and Chazdon, R. L. (2008). A Two-stage probabilistic approach to multiple-community similarity indices. Biometrics, 64, 1178-1186. \cr\cr
  Janzen, D. H. (1973) Sweep samples of tropical foliage insects: description of study sites, with data on species abundances and size distributions   . Ecology, 54, 659-686. \cr\cr
  Longino, J. T., Coddington, J. A. and Colwell, R. K. (2002). The ant fauna of a tropical rain forest: estimating species richness three different   ways. Ecology, 83, 689-702. \cr\cr
  Butler, B. J., and Chazdon, R. L. (1998). Species richness, spatial variation, and abundance of the
  soil seed bank of a secondary tropical rain forest. Biotropica, 30, 214-222. \cr\cr
}

\keyword{datasets}
