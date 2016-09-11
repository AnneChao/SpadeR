\name{ChaoSpeciesData}
\alias{ChaoSpeciesData}
\docType{data}
\title{
 Data for Function ChaoSpecies
}
\description{
 There are five data sets: \cr\cr
 1. Type (1) abundance data (\code{ChaoSpeciesData$Abu}) \cr\cr
 The data consist of 25 birds abundances/frequencies in a sample (Magurran, 1988, p.152). Their observed frequencies are respectively 752, 276, 194, 126, 121, 97, 95, 83, 72, 44, 39, 16, 15, 13, 9, 9, 9, 8, 7, 4, 2, 2, 1, 1, 1. \cr\cr
 2. Type (1A) abundance-frequency counts data (\code{ChaoSpeciesData$Abu_count}) \cr\cr
 The data consist of the observed species abundance distribution of endangered and rare vascular plant species in the central portion of the southern Appalachian region (Miller and Wiegert, 1989). A total of 188 species were recorded out of 1008 individuals compiled over a span of 150 years of field observations. The data are read as: (1 61 2 35 3 18 4 12 ... 67 1); each number needs to be separated by at least one blank space or by separated by rows. Here the first pair (1, 61) indicates that there are 61 singletons, the second pair (2, 35) indicates there are 35 doubletons, and so on, with the last pair (67, 1) indicating that there is one species that is represented by 67 individuals. \cr\cr
 3. Type (2) incidence-frequency data (\code{ChaoSpeciesData$Inci}) \cr\cr
 The data include seed-bank records taken from Butler and Chazdon (1998). There were 121 soil samples (each soil sample is regarded as a sampling unit) and species of seedlings that germinated from each soil sample were recorded. A total of 34 species of seedlings were found in the 121 soil samples. In the input data, the entry in the first row denotes the number of sampling units. Then, beginning with the second row, each row records the species incidence frequency (i.e., the number of soil samples in which the seedlings were found) of a given species in all 121 soil samples. The ordering of data entries does not affect the analysis. \cr\cr
 4. Type (2A) incidence-frequency counts data (\code{ChaoSpeciesData$Inci_freq_count}) \cr\cr
 The data consist of cottontail capture-recapture data provided in Edwards and Eberhardt (1967) to illustrate that species richness estimation can be applied to estimate the size of a population. An "individual" animal in capture-recapture studies corresponds to a "species" in the richness estimation. A total of 142 captures were recorded for 76 distinct rabbits in 18 trapping nights. For these data, the incidence frequency counts (\eqn{Q_1} to \eqn{Q_7}) were 43, 16, 8, 6, 0, 2, 1. The input data are read as follows: 
 (18 1 43 2 16 3 8 4 6 6 2 7 1); each number needs to be separated by at least one blank space or separated by rows. Here the pair (1, 43) indicates that there are 43 unique species, the next pair (2, 16) indicates there are 16 duplicate species, and so on. \cr\cr
 5. Type (2B) incidence-raw data (\code{ChaoSpeciesData$Inci_raw}) \cr\cr
 In the cottontail capture-recapture experiments conducted by Edwards and Eberhardt (1967), a total of 76 distinct individuals (regarded as 76 "species") were found in 18 trapping nights. The incidence-raw data consist of a capture/non-capture matrix (where "1" means a capture and "0" means a non-capture) with 76 rows and 18 columns. If we regard this capture-recapture matrix as a species-by-sampling-unit matrix, then species richness estimation can be applied to estimate the size of the cottontail population.
}
\usage{data(ChaoSpeciesData)}

\references{
  Magurran, A. E. (1988). Ecological Diversity and Its Measurement. 
  Princeton University Press, Princeton, New Jersey. \cr\cr
  Miller, R. I. and Wiegert, R. G. (1989). Documenting completeness, species-area relations, and the species-abundance distribution of a regional     flora. Ecology, 70, 16-22. \cr\cr
  Butler, B. J., and Chazdon, R. L. (1998). Species richness, spatial variation, and abundance of
  the soil seed bank of a secondary tropical rain forest. Biotropica, 30, 214-222.\cr\cr
  Edwards, W. R. and Eberhardt, L. (1967). Estimating cottontail abundance from live trapping data. The Journal of Wildlife Management, 31, 87-96. 
}

\keyword{datasets}
