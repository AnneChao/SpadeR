\name{SimilarityMultData}
\alias{SimilarityMultData}
\docType{data}
\title{
 Data for Function SimilarityMult
}
\description{
 There are three data sets: \cr\cr
1. Type (1) abundance data (\code{SimilarityMultData$Abu}) \cr\cr
The data include the observed species frequencies of three communities: seedlings (column 1), saplings (column 2) and trees (column 3) collected from an old-growth rain forest; see Chao et al. (2005, 2008). The three entries in each row are the observed frequency (or abundance) of each species from the three communities. \cr\cr
2. Type (2) incidence-frequency data (\code{SimilarityMultData$Inci}) \cr\cr
 The data include the observed incidence frequencies of tropical rainforest ants using three sampling techniques: (a) Berlese extraction of soil samples (217 samples), (b) fogging samples from canopy fogging (459 samples), and (c) Malaise trap samples for flying and crawling insects (62 samples); The data were collected in Costa-Rica (Longino et al. 2002). The three entries in the first row of the input data denote the number of sampling units (217, 459 and 62). Beginning with the second row, the three numbers in each row denotes incidence frequencies (the total number of detections) in the samples based on three sampling techniques. \cr\cr
3. Type (2B) incidence-raw data (\code{SimilarityMultData$Inci_raw}) \cr\cr
 The data include the observed soil ciliate species detection/non-detection data for a total of 51 soil samples from three areas of Namibia, Africa: Etosha Pan (19 samples), Central Namib Desert (17 samples) and Southern Namib Desert (15 samples). The raw detection/non-detection data include 365 x 51 matrix of 0's and 1's (0 denotes a non-detection and 1 denotes a detection).
}
\usage{data(SimilarityMultData)}

\references{
  Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters, 8, 148-159.\cr\cr
  Chao, A., Jost, L., Chiang, S.-C., Jiang, Y.-H. and Chazdon, R. L. (2008). A Two-stage probabilistic approach to multiple-community similarity indices. Biometrics, 64, 1178-1186. \cr\cr
  Longino, J. T., Coddington, J. A. and Colwell, R. K. (2002). The ant fauna of a tropical rain forest: estimating species richness three different ways. Ecology, 83, 689-702. \cr\cr
  Foissner, W., Agatha, S. and Berger, H. (2002). Soil Ciliates (Protozoa, Ciliophora) from Namibia (Southwest Africa), with emphasis on two contrasting environments, the Etosha Region and the Namib Desert. Denisia, 5, 1-1459.
}
\keyword{datasets}
