\name{SpadeR-package}
\alias{SpadeR-package}
\alias{SpadeR}
\docType{package}
\title{Species-richness prediction and diversity estimation with R
}
\description{
Provides simple functions to compute various biodiversity indices and related (dis)similarity measures based on individual-based (abundance) data or sampling-unit-based (incidence) data taken from one or multiple communities/assemblages. \cr\cr

This package contains six main functions: \cr\cr
1. \code{ChaoSpecies} (estimating species richness for one community). \cr\cr
2. \code{Diversity} (estimating a continuous diversity profile and various diversity indices in one community including species richness, Shannon
diversity and Simpson diversity). This function also features plots of empirical and estimated continuous diversity profiles. \cr\cr
3. \code{ChaoShared} (estimating the number of shared species between two communities). \cr\cr
4. \code{SimilartyPair} (estimating various similarity indices between two assemblages). Both richness- and abundance-based two-community similarity indices are included. \cr\cr
5. \code{SimilarityMult} (estimating various similarity indices among \eqn{N} communities). Both richness- and abundance-based \eqn{N}-community similarity indices are included. \cr\cr
6. \code{Genetics} (estimating allelic dissimilarity/differentiation among sub-populations based on multiple-subpopulation genetics data). \cr\cr
Except for the \code{Genetics} function, there are at least three types of data are supported for each function. \cr\cr
}
\details{
Data are generally classified as abundance data and incidence data and there are five types of data input formats options (datatype="abundance", "abundance_freq_count", "incidence_freq", "incidence_freq_count", "incidence_raw"). \cr\cr

 \describe{
  \item{A.}{Individual-based abundance data when a sample of individuals is taken from each community.}
        \bold{Type (1) abundance data} (datatype = "abundance"): Input data consist of species (in rows) by community (in columns) matrix. The entries of each row are the observed abundances of a species in \eqn{N} communities. \cr\cr
        \bold{Type (1A) abundance-frequency counts data} only for a single community (datatype = "abundance_freq_count"): input data are arranged as (1 \eqn{f_1 \ 2 \ f_2 \ ... \ r \ f_r})(each number needs to be separated by at least one blank space or separated by rows), where \eqn{r} denotes the maximum frequency and \eqn{f_k} denotes the number of species represented by exactly \eqn{k} individuals/times in the sample. Here the data (\eqn{f_1, f_2, ..., f_r}) are referred to as "abundance-frequency counts".
  \item{B.}{Sampling-unit-based incidence data when a number of sampling units are randomly taken from each community. Only the incidence (detection/non-detection) of species is recorded in each sampling unit. There are three data formats options.} \cr\cr
        \bold{Type (2) incidence-frequency data} (datatype="incidence_freq"): The first row of the input data must be the number of sampling units in each community. Beginning with the second row, input data consist of species (in rows) by community (in columns) matrix. The entries of each row are the observed incidence frequencies (the number of detections or the number of sampling units in which a species are detected) of a species in \eqn{N} communities. \cr\cr
        \bold{Type (2A) incidence-frequency counts data} only for a single community (datatype="incidence \cr _freq_count"): input data are arranged as  (\eqn{T \ 1 \ Q_1 \ 2 \ Q_2 \ ... \ r \ Q_r}) (each number needs to be separated by at least one blank space or separated by rows), where \eqn{Q_k} denotes the number of species that were detected in exactly \eqn{k} sampling units, while \eqn{r} denotes the number of sampling units in which the most frequent species were found. The first entry must be the total number of sampling units, \eqn{T}. The data (\eqn{Q_1, Q_2, ..., Q_r}) are referred to as "incidence frequency counts". \cr\cr
        \bold{Type (2B) incidence-raw data} (datatype="incidence_raw"): Data consist of a species-by-sampling-unit incidence (detection/non-detection) matrix; typically "1" means a detection and "0" means a non-detection. Each row refers to the detection/non-detection record of a species in \eqn{T} sampling units. Users must specify the number of sampling units in the function argument "units". The first \eqn{T_1} columns of the input matrix denote species detection/non-detection data based on the \eqn{T_1} sampling units from Community 1, and the next \eqn{T_2} columns denote the detection/non-detection data based on the \eqn{T_2} sampling units from Community 2, and so on, and the last \eqn{T_N} columns denote the detection/non-detection data based on \eqn{T_N} sampling units from Community \eqn{N}, \eqn{T_1 + T_2 + ... + T_N = T}. \cr
}

An Online version of SpadeR is also available for users without an R background: \cr
http://chao.stat.nthu.edu.tw/wordpress/software_download/softwarespader_online/. \cr
In the detailed Online SpadeR User's Guide, we illustrate all the running procedures in an easily
accessible way through numerical examples with proper interpretations of portions of the output.
All the data of those illustrative examples are included in this package. \cr\cr

functions: ChaoSpecies, Diversity, ChaoShared, SimilarityPair, SimilarityMult, Genetics
}
\author{
Anne Chao, K. H. Ma, T. C. Hsieh and Chun-Huo Chiu 

Maintainer: Anne Chao <chao@stat.nthu.edu.tw>

}
\keyword{ package }


