\name{SimilarityPairData}
\alias{SimilarityPairData}
\docType{data}
\title{
 Data for Function SimilarityPair
}
\description{
 There are three data sets: \cr\cr
 1. Type (1) abundance data (\code{SimilarityPairData$Abu})  \cr\cr
 The data include the observed species frequencies of two communities: seedlings (column 1), and trees (column 2) collected from an old-growth rain forest; see Chao et al. (2005, 2008). The two entries in each row are the observed frequency (or abundance) of each species from the two communities. (These data are subset of \code{SimilarityMultData$Abu} used in the function SimilarityMult.) \cr\cr
 2. Type (2) incidence-frequency data (\code{SimilarityPairData$Inci}) \cr\cr
 The data include the observed incidence frequencies of tropical rainforest ants based on two sampling techniques: (a) Berlese extraction of soil samples (217 samples), and (b) Malaise trap samples for flying and crawling insects (62 samples); see Longino et al. (2002). The two entries in first row of the input data denote the number of sampling units (217 and 62). Beginning with the second row, the two numbers in each row denotes incidence frequencies (the total number of detections) in the soil samples based on the two sampling techniques. (These data are subset of \code{SimilarityMultData$Inci} used in the function SimilarityMult.)  \cr\cr
 3. Type (2B) incidence-raw data (\code{SimilarityPairData$Inci_raw}) \cr\cr
 The data include the observed soil ciliate species detection/non-detection data for a total of 36 soil samples from two areas of Namibia, Africa: Etosha Pan (19 samples), and Central Namib Desert (17 samples). The raw detection/non-detection data include 365 x 36 matrix of 0's and 1's (0 denotes a non-detection and 1 denotes a detection). (These data are subset of \code{SimilarityMultData$Inci_raw} used in the function SimilarityMult.)
}
\usage{data(SimilarityPairData)}

\references{
  Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters, 8, 148-159.\cr\cr
  Chao, A., Jost, L., Chiang, S.-C., Jiang, Y.-H. and Chazdon, R. L. (2008). A Two-stage probabilistic approach to multiple-community similarity indices. Biometrics, 64, 1178-1186.\cr\cr
  Longino, J. T., Coddington, J. A. and Colwell, R. K. (2002). The ant fauna of a tropical rain forest: estimating species richness three different ways. Ecology, 83, 689-702. \cr\cr
  Foissner, W., Agatha, S. and Berger, H. (2002) Soil Ciliates (Protozoa, Ciliophora) from Namibia (Southwest Africa), with emphasis on two Contrasting environments, the Etosha Region and the Namib Desert. Denisia, 5, 1-1459.
}

\keyword{datasets}
