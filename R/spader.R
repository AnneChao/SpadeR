#
#
###########################################
#' Estimation of species richness
#' 
#' \code{ChaoSpecies}: Estimation of species richness in one assemblage. \cr\cr
#' Five types of data (format/information) are supported: ("abundance","abundance_freq_count", "incidence_freq", "incidence_freq_count", "incidence_raw"):
#' \enumerate{
#'   \item{Individual-based abundance data: There are two kinds of input data.
#'   \describe{
#'      \item{(1)}{Abundance data (datatype="abundance"): input data is a vector of sample species abundances/frequencies for observed species in an empirical sample of individuals.}
#'      \item{(1A)}{Abundance-frequency counts (datatype="abundance_freq_count"): input data is a vector of sample species frequency counts.}
#'    }
#'   }
#'   \item{Sampling-unit-based incidence data: There are three kinds of input data.
#'   \describe{
#'      \item{(2)}{Incidence-frequency data (datatype="incidence_freq"): input data consists of species sample incidence frequencies (row sums of each incidence matrix). The first entry must be the total number of sampling units, followed by the species incidence frequencies.}
#'      \item{(2A)}{Incidence-frequency counts (datatype="incidence_freq_count"): input data is a vector of sample species incidence counts. The first entry must be the total number of sampling units.}
#'      \item{(2B)}{Incidence-raw data (datatype="incidence_raw"): input data for a reference sample consists of a species-by-sampling-unit matrix.}
#'    }
#'   }
#' }
#' @param data a matrix, data.frame of species abundances/incidences (see data format/information above).\cr
#' @param datatype type of input data, "abundance", "abundance_freq_count", "incidence_freq" or "incidence_freq_count", "incidence_raw". \cr
#' @param k cut-off point; it is a value that separates frequency counts into abundant and rare groups.
#' @param conf a positive number \eqn{\le} 1 specifying the level of confidence interval.
#' @return  a list of three objects: \cr\cr
#' \code{$Basic.Data.Information} and \code{$Rare.Species.Group}/\code{$Infreq.Species.Group} for summarizing data information. \cr\cr
#' \code{$Species.Table} for showing a table of various species richness estimates, standard errors, and the associated confidence intervals. \cr\cr
#' @examples
#' data(ChaoSpeciesDataAbu)
#' ChaoSpecies(ChaoSpeciesDataAbu, datatype="abundance", k=10, conf=0.95)
#' data("ChaoSpeciesDataAbu_count")
#' ChaoSpecies(ChaoSpeciesDataAbu_count,datatype="abundance_freq_count",k=10,conf=0.95)
#' data(ChaoSpeciesDataInci)
#' ChaoSpecies(ChaoSpeciesDataInci, datatype="incidence_freq", k=10, conf=0.95)
#' data("ChaoSpeciesDataInci_freq_count")
#' ChaoSpecies(ChaoSpeciesDataInci_freq_count, datatype="incidence_freq_count",k=10,conf=0.95)
#' data(ChaoSpeciesDataInci_raw)
#' ChaoSpecies(ChaoSpeciesDataInci_raw, datatype="incidence_raw", k=10, conf=0.95)
#' @references
#' Burnham, K. P. and Overton, W. S. (1978). Estimation of the size of a closed population when capture probabilities vary among animals. Biometrika 65, 625-633.\cr\cr
#' Chao, A. (1984). Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics 11, 265-270.\cr\cr
#' Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.\cr\cr
#' Chao, A. (2005). Species estimation and applications. Encyclopedia of Statistical Sciences, Second Edition, Vol. 12, 7907-7916 (N. Balakrishnan, C. B. Read and B. Vidakovic, Editors), Wiley, New York.\cr\cr
#' Chao, A., Hwang, W.-H., Chen, Y.-C. and Kuo. C.-Y. (2000). Estimating the number of shared species in two communities. Statistica Sinica 10, 227-246.\cr\cr
#' Chao, A. and Lee, S.-M. (1992). Estimating the number of classes via sample coverage. Journal of the American Statistical Association 87, 210-217.\cr\cr
#' Chao, A., Ma, M-C. and Yang, M. C. K. (1993). Stopping rule and estimation for recapture debugging with unequal detection rates. Biometrika 80, 193-201.\cr\cr
#' Chao, A., Shen, T.-J. and Hwang, W. H. (2006). Application of Laplace boundary-mode approximations to estimate species and shared species richness. Australian and New Zealand Journal of Statistics 48, 117-128.\cr\cr
#' Chiu, C.-H., Wang Y. T., Walther B. A. and Chao A. (2014). An improved non-parametric lower bound of species richness via the Good-Turing frequency formulas. Biometrics 70, 671-682.\cr\cr
#' Lee, S.-M. and Chao, A. (1994). Estimating population size via sample coverage for closed capture-recapture models. Biometrics 50, 88-97.\cr\cr
#' Shen, T.-J. (2003). Prediction of Biodiversity, Ph.D. dissertation, National Tsing-Hua University, HsinChu, Taiwan.\cr\cr
#' Shen, T.-J., Chao, A. and Lin, J.-F. (2003). Predicting the number of new species in further taxonomic sampling. Ecology 84, 798-804.
#' @export


ChaoSpecies <- function(data, datatype = c("abundance","abundance_freq_count", "incidence_freq", "incidence_freq_count", "incidence_raw"), k = 10, conf = 0.95)
{
  if (is.matrix(data) == T || is.data.frame(data) == T){
    #if (ncol(data) != 1 & nrow(data) != 1)
    #stop("Error: The data format is wrong.")
    if(datatype != "incidence_raw"){
      if (ncol(data) == 1){
        data <- data[, 1]
      } else {
        data <- data[1, ]
      }
    } else{
      t <- ncol(data)
      dat <- rowSums(data)
      dat <- as.integer(dat)
      t_infreq <- sum(colSums(data[which(dat<k),])>=1)
      data <- dat
      data <- c(t_infreq, t , data)
    }
    
  }
  
  if(datatype == "abundance_freq_count"){
    data <- as.integer(data)
    length4b <- length(data)
    data <- rep(data[seq(1,length4b,2)],data[seq(2,length4b,2)])
    names(data) <- paste("x",1:length(data),sep="")
    datatype <- "abundance"
  }
  if (datatype == "incidence_freq_count"){
    t <- as.integer(data[1])
    data <- data[-c(1)]
    data <- as.integer(data)
    lengthdat <- length(data)
    data <- rep(data[seq(1,lengthdat,2)],data[seq(2,lengthdat,2)])
    data <- c(t,data)
    names(data) <- c("T", paste("y",1:(length(data)-1),sep=""))
    datatype <- "incidence_freq"
  }
  method <- "all"
  if (k != round(k) || k < 0)
    stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
  if (is.numeric(conf) == FALSE || conf > 1 || conf < 0)
    stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  #  data <- as.numeric(round(data))
  
  if (datatype == "abundance"){
    f <- function(i, data){length(data[which(data == i)])}
    if (f(1, data) == sum(data)){
      stop("Error: The information of data is not enough.")}
    z <- (list(Basic.Data.Information = basicAbun(data, k)[[1]], Rare.Species.Group = RareSpeciesGroup(data, k),
               Species.Table = round(SpecAbunOut(data, method, k, conf), 3)))
  } else if (datatype == "incidence_raw"){
    dat <- data[-1]; Q <- function(i, data){length(data[which(data == i)])}
    if (Q(1, dat) == sum(dat)){
      stop("Error: The information of data is not enough.")}
    z <- (list(Basic.Data.Information = basicInci(data[-1], k)[[1]], Infreq.Species.Group = InfreqSpeciesGroup(data[-1], k),
               Species.Table = round(SpecInciOut_raw(data, method, k, conf),3)))
  } else if (datatype == "incidence_freq"){
    dat <- data[-1];
    Q <- function(i, data){length(data[which(data == i)])}
    if (Q(1, dat) == sum(dat)){
      stop("Error: The information of data is not enough.")}
    z <- (list(Basic.Data.Information = basicInci(data, k)[[1]], Infreq.Species.Group = InfreqSpeciesGroup(data, k),
               Species.Table = round(SpecInciOut(data, method, k, conf),3)))
  } 
  else{
    stop("Error: The data type is wrong.")
  }
  class(z) <- c("ChaoSpecies")
  z
}


#
#
###########################################
#' Estimation of the number of shared species in two assemblages
#'
#' \code{ChaoShared}: Estimation of shared species richness in two assemblages based on the following two sampling schemes. \cr\cr
#' Three types of data (format/information) are supported: ("abundance", "incidence_freq", "incidence_raw"):
#' \enumerate{
#'    \item{(1) Individual-based abundance data (datatype="abundance"): Input data for each assemblage/site include sample species abundances in an empirical sample of n individuals ("reference sample"). There are 2 assemblages, input data consist of an S by 2 abundance matrix, or 2 lists of species abundances.}    
#'   \item{Sampling-unit-based incidence data: There are two kinds of input data.
#'   \describe{
#'      \item{(2)}{Incidence-raw data (datatype="incidence_raw"): for each assemblage, input data for a reference sample consist of a species-by-sampling-unit matrix; there are 2 assemblages, input data consist of 2 lists of matrices, and each matrix is a species-by-sampling-unit matrix.}
#'      \item{(2B)}{Incidence-frequency data (datatype="incidence_freq"): input data for each assemblage consist of species sample incidence frequencies (row sums of each incidence matrix). There are 2 assemblages, input data consist of an S+1 by 2 matrix, or 2 lists of species incidence frequencies. The first entry of each column/list must be the total number of sampling units, followed by the species incidence frequencies.}
#'    }
#'   }
#' }
#' @param data a matrix, data.frame, lists of species abundances/incidences, or lists of incidence frequencies (see data format/information above).\cr
#' @param datatype type of input data, "abundance", "incidence_freq" or "incidence_raw". \cr
#' @param se a logical variable to calculate the bootstrap standard error and the associated confidence interval. \cr
#' @param nboot an integer specifying the number of bootstrap replications. \cr
#' @param conf a positive number \eqn{\le} 1 specifying the level of confidence interval.
#' @return a list of two objects: \cr\cr
#' \code{$BASIC_DATA_INFORMATION} for summarizing data information. \cr\cr
#' \code{$ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES} for showing a table of various shared richess estimates, standard errors, and the associated confidence intervals. \cr\cr
#' @examples
#' data(ChaoSharedDataAbu)
#' ChaoShared(ChaoSharedDataAbu, datatype="abundance", se=TRUE, nboot=200, conf=0.95)
#' data(ChaoSharedDataInci)
#' ChaoShared(ChaoSharedDataInci, datatype="incidence_freq", se=TRUE, nboot=200, conf=0.95)
#' data(ChaoSharedDataInci_raw)
#' ChaoShared(ChaoSharedDataInci_raw, datatype="incidence_raw", se=TRUE, nboot=200, conf=0.95)
#' @references
#' Chao, A., Hwang, W.-H., Chen, Y.-C. and Kuo. C.-Y. (2000). Estimating the  number of shared species in two communities. Statistica Sinica 10, 227-246.\cr\cr
#' Chao, A., Shen, T.-J. and Hwang, W.-H. (2006). Application of Laplace boundary-mode approximations to estimate species and shared species richness.  Australian and New Zealand Journal of Statistics 48, 117-128.\cr\cr
#' Pan, H. Y., Chao, A. and Foissner, W. (2009). A non-parametric lower bound for the number of species shared by multiple communities. Journal of Agricultural, Biological and Environmental Statistics 14, 452-468.
#' @export

ChaoShared <-
  function(data, datatype = c("abundance", "incidence_raw", "incidence_freq"),
           se = TRUE, nboot = 200, conf = 0.95) {

    method <- "all"
    if (se == TRUE) {
      if (nboot < 1)
        nboot <- 1
      if (nboot == 1)
        cat("Warning: When \"nboot\" =" ,nboot, ", the bootstrap s.e. and confidence interval can't be calculated.",
            "\n\n")
    }

    if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) {
      cat("Warning: \"conf\"(confidence level) must be a numerical value between 0 and 1, e.g. 0.95.",
          "\n")
      cat("          We use \"conf\" = 0.95 to calculate!",
          "\n\n")
      conf <- 0.95
    }

    datatype <- match.arg(datatype)
    if (datatype == "abundance") {
      if(class(data)=="list"){data <-cbind(data[[1]],data[[2]]) }
      x1 <- data[, 1]
      x2 <- data[, 2]
      Basic <- BasicFun(x1, x2, nboot, datatype)
      #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
      output <- ChaoShared.Ind(x1, x2, method, nboot, conf, se)
      colnames(output) <- c("Estimate", "s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
    }
    if (datatype == "incidence_freq") {
      if(class(data)=="list"){data <-cbind(data[[1]],data[[2]]) }
      y1 <- data[, 1]
      y2 <- data[, 2]
      Basic <- BasicFun(y1, y2, B=nboot, datatype)
      #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
      output <- ChaoShared.Sam(y1, y2, method, conf, se)
      colnames(output) <- c("Estimate", "s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
    }
    if (datatype=="incidence_raw"){
      a <- data[[1]]
      b <- data[[2]]
      t1 <- ncol(a) ; t2 <- ncol(b)
      data1 <- as.integer(rowSums(a)) ; data2 <- as.integer(rowSums(b))
      y1 <- c(t1, data1)
      y2 <- c(t2, data2)
      datatype = "incidence_freq"
      Basic <- BasicFun(y1, y2, B=nboot, datatype)
      output <- ChaoShared.Sam(y1, y2, method, conf, se)
      colnames(output) <- c("Estimate", "s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
    }
    out <- list(BASIC_DATA_INFORMATION=Basic,
                ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES=output)
    class(out) <- c("ChaoShared")
    return(out)
  }


#
#
###########################################
#' Estimation of species diversity
#'
#' \code{Diversity}: Estimation of various diversity indices including species richness (diversity of order 0),
#' the Shannon index/diversity (diversity of order 1), the Simpson index/diversity (diversity of order 2),
#' and Hill number (diversity of order from 0 to 3). \cr\cr
#' Five types of data (format/information) are supported: ("abundance","abundance_freq_count", "incidence_freq", "incidence_freq_count", "incidence_raw"):
#' \enumerate{
#'   \item{Individual-based abundance data: There are two kinds of input data.
#'   \describe{
#'      \item{(1)}{Abundance data (datatype="abundance"): input data is a vector of sample species abundances/frequencies for observed species in an empirical sample of individuals.}
#'      \item{(1A)}{Abundance-frequemcy counts (datatype="abundance_freq_count"): input data is a vector of sample species frequency counts.}
#'    }
#'   }
#'   \item{Sampling-unit-based incidence data: There are three kinds of input data.
#'   \describe{
#'      \item{(2)}{Incidence-frequency data (datatype="incidence_freq"):  input data is species sample incidence frequencies (row sums of each incidence matrix). The first entry must be the total number of sampling units.}
#'      \item{(2A)}{Incidence-frequency counts (datatype="incidence_freq_count"): input data is a vector of sample species incidence counts. The first entry must be the total number of sampling units.}
#'      \item{(2B)}{Incidence-raw data (datatype="incidence_raw"): input data for a reference sample consists of a species-by-sampling-unit matrix.}
#'    }
#'   }
#' }
#' @param data a matrix, data.frame of species abundances/incidences (see data format/information above).\cr
#' @param q a vector of nonnegative numbers for which the diversity order of Hill numbers will be estimated. If \code{NULL}, then
#' Hill numbers will be estimated at order q from 0 to 3 with equally-spaced 0.25.
#' @param datatype type of input data, "abundance", "abundance_freq_count", "incidence_freq", "incidence_freq_count" or "incidence_raw". \cr\cr
#' @return a list of seven objects: \cr\cr
#' \code{$BASIC.DATA} for summarizing data information. \cr\cr
#' \code{$SPECIES.RICHNESS} for showing various species richness estimates along with related statistics. \cr\cr
#' \code{$SHANNON.INDEX} and \code{$EXPONENTIAL.OF.SHANNON.INDEX} for showing various Shannon index/diversity estimates. \cr\cr
#' \code{$SIMPSON.INDEX} and \code{$INVERSE.OF.SIMPSON.INDEX} for showing various Simpson index/diversity estimates. \cr\cr
#' \code{$HILL.NUMBERS} for showing Hill number (diversity of order from 0 to 3) estimates. \cr\cr
#' @examples
#' \dontrun{
#' data(DiversityDataAbu)
#' Diversity(DiversityDataAbu,datatype="abundance")
#' data("DiversityDataAbu_count")
#' Diversity(DiversityDataAbu_count, datatype="abundance_freq_count", q=NULL)
#' data(DiversityDataInci)
#' Diversity(DiversityDataInci, datatype="incidence_freq", q=NULL)
#' data("DiversityDataInci_freq_count")
#' Diversity(DiversityDataInci_freq_count, datatype="incidence_freq_count", q=NULL)
#' data(DiversityDataInci_raw)
#' Diversity(DiversityDataInci_raw, datatype="incidence_raw", q=NULL)
#' }
#' @references
#' Chao, A. (1984). Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics 11, 265-270.\cr\cr
#' Chao, A. (1987). Estimating the population size for capture-recapture data with  unequal catchability. Biometrics 43, 783-791.\cr\cr
#' Chao, A. (2005). Species estimation and applications. Encyclopedia of Statistical Sciences, Second Edition, Vol. 12, 7907-7916 (N. Balakrishnan, C. B. Read and B. Vidakovic, Editors), Wiley, New York.\cr\cr
#' Chao, A., Hwang, W.-H., Chen, Y.-C. and Kuo. C.-Y. (2000). Estimating the number of shared species in two communities. Statistica Sinica 10, 227-246.\cr\cr
#' Chao, A. and Jost, L. (2015). Estimating diversity and entropy profiles via discovery rates of new species. Methods in Ecology and Evolution.\cr\cr
#' Chao, A. and Lee, S.-M. (1992). Estimating the number of classes via sample coverage. Journal of the American Statistical Association 87, 210-217.\cr\cr
#' Chao, A., Ma, M-C. and Yang, M. C. K. (1993). Stopping rule and estimation for recapture debugging with unequal detection rates. Biometrika 80, 193-201.\cr\cr
#' Chao, A. and Shen, T.-J. (2003). Nonparametric estimation of Shannon index of diversity when there are unseen species in sample. Environmental and Ecological Statistics 10, 429-443.\cr\cr
#' Chao, A., Shen, T.-J. and Hwang, W. H. (2006). Application of Laplace boundary-mode approximations to estimate species and shared species richness. Australian and New Zealand Journal of Statistics 48, 117-128.\cr\cr
#' Chao, A., Wang, Y. T. and Jost, L. (2013). Entropy and the species accumulation curve: a novel estimator of entropy via discovery rates of new species. Methods in Ecology and Evolution 4, 1091-1110.\cr\cr
#' Magurran, A. E. (1988). Ecological Diversity and Its Measurement. Princeton, Princeton University Press, New Jersey.\cr\cr
#' Shen, T.-J., Chao, A. and Lin, J.-F. (2003). Predicting the number of new species in further taxonomic sampling. Ecology 84, 798-804.\cr\cr
#' Zahl, S. (1977). Jackknifing an index of diversity. Ecology 58, 907-913.
#' @export

Diversity=function(data, datatype=c("abundance","abundance_freq_count", "incidence_freq", "incidence_freq_count", "incidence_raw"), q=NULL)
{
  if (is.matrix(data) == T || is.data.frame(data) == T){
    #if (ncol(data) != 1 & nrow(data) != 1)
    #stop("Error: The data format is wrong.")
    if(datatype != "incidence_raw"){
      if (ncol(data) == 1){
        data <- data[, 1]
      } else {
        data <- as.vector(data[1, ])
      }
    } else{
      t <- ncol(data)
      dat <- rowSums(data)
      dat <- as.integer(dat)
      data <- c(t , dat)
    }
    
  }
  X <- data
  if(datatype == "abundance_freq_count"){
    data <- as.integer(data)
    length4b <- length(data)
    data <- rep(data[seq(1,length4b,2)],data[seq(2,length4b,2)])
    names(data) <- paste("x",1:length(data),sep="")
    datatype <- "abundance"
    X <- data
  }
  if(datatype=="abundance"){
    type="abundance"
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))

    BASIC.DATA <- matrix(round(c(sum(X), sum(X>0), 1-sum(X==1)/sum(X), CV.Ind(X)),3), ncol = 1)
    nickname <- matrix(c("n", "D", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)

    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("    Sample size", "    Number of observed species",
                              "    Estimated sample coverage",
                              "    Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)

    table0 <- matrix(0,5,4)
    table0[1,]=c(Chao1(X)[-5])
    table0[2,]=c(Chao1_bc(X))
    table0[3,]=round(SpecAbuniChao1(X, k=10, conf=0.95)[1,],1)
    table0[4,]=round(c(SpecAbunAce(X)),1)
    table0[5,]=round(c(SpecAbunAce1(X)),1)
    colnames(table0) <- c("Estimate", "s.e.", paste(Chao1(X)[5]*100,"%Lower", sep=""), paste(Chao1(X)[5]*100,"%Upper", sep=""))
    rownames(table0) <- c("    Chao1 (Chao, 1984)","    Chao1-bc ", "    iChao1","    ACE (Chao & Lee, 1992)",
                          "    ACE-1 (Chao & Lee, 1992)")

    SHANNON=Shannon_index(X)
    table1=round(SHANNON[c(1:5),],3)
    table1=table1[-2,]              ##2016.05.09
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Jackknife",
    #                      " Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("     MLE","     Jackknife",
                          "     Chao & Shen","     Chao et al. (2013)")

    table1_exp=round(SHANNON[c(6:10),],3)
    table1_exp=table1_exp[-2,]      ##2016.05.09
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Jackknife",
    #                         " Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("     MLE","     Jackknife",
                              "     Chao & Shen","     Chao et al. (2013)")

    table2=round(Simpson_index(X)[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("     MVUE","     MLE")

    table2_recip=round(Simpson_index(X)[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("     MVUE","     MLE")

    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", q=NULL, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "abundance", q=q, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2

    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- nrow(Hill)
    rownames(Hill) <- paste("    ",1:q_hill)
    z <- list("datatype"= type,"BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0,
              "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
              "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
              "HILL.NUMBERS"= Hill)
  }
  if(datatype == "incidence_freq_count"){
    t <- as.integer(data[1])
    data <- data[-c(1)]
    data <- as.integer(data)
    lengthdat <- length(data)
    data <- rep(data[seq(1,lengthdat,2)],data[seq(2,lengthdat,2)])
    data <- c(t,data)
    names(data) <- c("T", paste("y",1:(length(data)-1),sep=""))
    datatype <- "incidence_freq"
    X <- data
  }
  if(datatype=="incidence_freq"){
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))
    type="incidence"
    U<-sum(X[-1])
    D<-sum(X[-1]>0)
    T<-X[1]
    C<-Chat.Sam(X,T)
    CV_squre<-max( D/C*T/(T-1)*sum(X[-1]*(X[-1]-1))/U^2-1, 0)
    CV<-CV_squre^0.5
    BASIC.DATA <- matrix(round(c(D,T,U, C, CV),3), ncol = 1)
    nickname <- matrix(c("D", "T","U", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)

    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("    Number of observed species", "    Number of Sampling units","    Total number of incidences",
                              "    Estimated sample coverage",
                              "    Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    #BASIC.DATA <- basicInci(X, k=10)[[1]]
    ############################################################
    table0=SpecInci(X, k=10, conf=0.95)
    rownames(table0) <- c("    Chao2 (Chao, 1987)","    Chao2-bc ", "    iChao2","    ICE (Lee & Chao, 1994)",
                          "    ICE-1 (Lee & Chao, 1994)")
    SHANNON=Shannon_Inci_index(X)
    table1=round(SHANNON[c(1,4),],3)
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("     MLE","     Chao et al. (2013)")
    table1_exp=round(SHANNON[c(5,8),],3)
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("     MLE","     Chao et al. (2013)")

    SIMPSON=Simpson_Inci_index(X)
    table2=round(SIMPSON[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("     MVUE","     MLE")

    table2_recip=round(SIMPSON[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("     MVUE","     MLE")


    ############################################################
    #Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", from=0, to=3, interval=0.25, B=50, conf=0.95))
    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence_freq", q=NULL, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence_freq", q=q, from=0, to=3, interval=0.25, B=50, conf=0.95))}

    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2

    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- nrow(Hill)
    rownames(Hill) <- paste("    ",1:q_hill)
    #z <- list("BASIC.DATA"=BASIC.DATA,"HILL.NUMBERS"= Hill)

    z <- list("datatype"= type,"BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0,
              "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
              "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
              "HILL.NUMBERS"= Hill)
  }
  if(datatype=="incidence_raw"){
    type="incidence"
    datatype = "incidence_freq"
    U<-sum(X[-1])
    D<-sum(X[-1]>0)
    T<-X[1]
    C<-Chat.Sam(X,T)
    CV_squre<-max( D/C*T/(T-1)*sum(X[-1]*(X[-1]-1))/U^2-1, 0)
    CV<-CV_squre^0.5
    BASIC.DATA <- matrix(round(c(D,T,U, C, CV),3), ncol = 1)
    nickname <- matrix(c("D", "T","U", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)

    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("    Number of observed species", "    Number of Sampling units","    Total number of incidences",
                              "    Estimated sample coverage",
                              "    Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    ############################################################
    table0=SpecInci(X, k=10, conf=0.95)
    rownames(table0) <- c("    Chao2 (Chao, 1987)","    Chao2-bc ", "    iChao2","    ICE (Lee & Chao, 1994)",
                          "    ICE-1 (Lee & Chao, 1994)")
    SHANNON=Shannon_Inci_index(X)
    table1=round(SHANNON[c(1,4),],3)
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c("     MLE","     Chao et al. (2013)")
    table1_exp=round(SHANNON[c(5,8),],3)
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c("     MLE","     Chao et al. (2013)")

    SIMPSON=Simpson_Inci_index(X)
    table2=round(SIMPSON[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c("     MVUE","     MLE")

    table2_recip=round(SIMPSON[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c("     MVUE","     MLE")


    ############################################################
    #Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", from=0, to=3, interval=0.25, B=50, conf=0.95))
    if(is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=NULL, from=0, to=3, interval=0.25, B=50, conf=0.95))}
    if(!is.null(q)){Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence", q=q, from=0, to=3, interval=0.25, B=50, conf=0.95))}

    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
    #Chao.LCL <- Hill[14:26,3] - 1.96*Hill[14:26,4]
    #Chao.UCL <- Hill[14:26,3] + 1.96*Hill[14:26,4]
    #Emperical.LCL <- Hill[1:13,3] - 1.96*Hill[1:13,4]
    #Emperical.UCL <- Hill[1:13,3] + 1.96*Hill[1:13,4]
    #Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Chao.LCL,Chao.UCL,Emperical.LCL,Emperical.UCL)
    #Hill<-round(Hill,3)
    #Hill <- data.frame(Hill)
    q_length<-length(Hill[,1])/2

    Chao.LCL <- Hill[(q_length+1):(2*q_length),3] - 1.96*Hill[(q_length+1):(2*q_length),4]
    Chao.UCL <- Hill[(q_length+1):(2*q_length),3] + 1.96*Hill[(q_length+1):(2*q_length),4]
    Emperical.LCL <- Hill[1:q_length,3] - 1.96*Hill[1:q_length,4]
    Emperical.UCL <- Hill[1:q_length,3] + 1.96*Hill[1:q_length,4]
    Hill<-cbind(Hill[1:q_length,1],Hill[(q_length+1):(2*q_length),3],Chao.LCL,Chao.UCL,Hill[1:q_length,3],Emperical.LCL,Emperical.UCL)
    Hill<-round(Hill,3)
    Hill <- data.frame(Hill)
    #colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")
    colnames(Hill)<-c("q","ChaoJost","95%Lower","95%Upper","Empirical","95%Lower","95%Upper")
    q_hill <- nrow(Hill)
    rownames(Hill) <- paste("   ",1:q_hill)
    #z <- list("BASIC.DATA"=BASIC.DATA,"HILL.NUMBERS"= Hill)

    z <- list("datatype"= type,"BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0,
              "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
              "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
              "HILL.NUMBERS"= Hill)
  }
  class(z) <- c("spadeDiv")
  return(z)
}







#
#
###########################################
#' Estimation of two-assemblage similarity index
#'
#' \code{SimilarityPair}: Estimation of various similarity indices based on either individual-based abundance data or sampling-unit-based incidence data
#' taken from two assemblages. The incidence-based indices include the classic Jaccard and Sorensen indices,
#' and the abundance-based indices include the Bray-Curtis, Morisita-Horn, Horn and abundance-based Jaccard
#' and Sorensen indices developed in Chao et al. (2005). \cr\cr
#' Three types of data (format/information) are supported: ("abundance", "incidence_freq", "incidence_raw"):
#' \enumerate{
#'   \item{(1) Individual-based abundance data (datatype="abundance"): Input data for each assemblage/site include sample species abundances in an empirical sample of n individuals ("reference sample"). There are 2 assemblages, input data consist of an S by 2 abundance matrix, or 2 lists of species abundances.}
#'   \item{Sampling-unit-based incidence data: There are two kinds of input data.
#'   \describe{
#'      \item{(2)}{Incidence-frequency data (datatype="incidence_freq"): input data for each assemblage consist of species sample incidence frequencies (row sums of each incidence matrix). There are 2 assemblages, input data consist of an S+1 by 2 matrix, or 2 lists of species incidence frequencies. The first entry of each column/list must be the total number of sampling units, followed by the species incidence frequencies.}
#'      \item{(2B)}{Incidence-raw data (datatype="incidence_raw"): for each assemblage, input data for a reference sample consist of a species-by-sampling-unit matrix; there are 2 assemblages, input data consist of 2 lists of matrices, and each matrix is a species-by-sampling-unit matrix.}
#'    }
#'   }
#' }
#' @param X a matrix, data.frame, lists of species abundances/incidences, or lists of incidence frequencies (see data format/information above).\cr
#' @param datatype type of input data, "abundance", "incidence_freq" or "incidence_raw". \cr
#' @param nboot an integer specifying the number of replications.
#' @return a list of six objects: \cr\cr
#' \code{$datatype} for abundance or incidence. \cr\cr
#' \code{$info1} and \code{$info2} for summarizing data information. \cr\cr
#' \code{$Empirical_incidence} for showing the estimation of classic richness-based similarity indices by empirical method. \cr\cr
#' \code{$Empirical_ew} for showing the estimation of similarity indices for comparing species relative abundances by empirical method. \cr \cr
#' \code{$Empirical_ee} for showing the estimation of similarity indices for comparing size-weighted species relative abundances and comparing species absolute abundances by empirical method. \cr\cr
#' \code{$estimated_incidence} for showing the estimation of classic richness-based similarity indices by estimated. \cr\cr
#' \code{$estimated_ew} for showing the estimation of similarity indices for comparing species relative abundances by estimated. \cr\cr
#' \code{$estimated_ee} for showing the estimation of similarity indices for comparing size-weighted species relative abundances and comparing species absolute abundances by estimated. \cr\cr
#' @examples
#' \dontrun{
#' data(SimilarityPairDataAbu)
#' SimilarityPair(SimilarityPairDataAbu, datatype="abundance",nboot=200)
#' data(SimilarityPairDataInci)
#' SimilarityPair(SimilarityPairDataInci, datatype="incidence_freq",nboot=200)
#' data(SimilarityPairDataInci_raw)
#' SimilarityPair(SimilarityPairDataInci_raw, datatype="incidence_raw",nboot=200)
#' }
#' @references
#' Chao, A. (1984). Nonparametric estimation of the number of classes in a population. Scandinavian Journal of Statistics 11, 265-270.\cr\cr
#' Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.\cr\cr
#' Chao, A. (2005). Species estimation and applications. Encyclopedia of Statistical Sciences, Second Edition, Vol. 12, 7907-7916 (N. Balakrishnan, C. B. Read and B. Vidakovic, Editors), Wiley, New York.\cr\cr
#' Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2005). A new statistical approach for assessing similarity of species composition with incidence and abundance data. Ecology Letters 8, 148-159.\cr\cr
#' Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2006). Abundance-based similarity indices and their estimation when there are unseen species in samples. Biometrics, 62, 361-371.\cr\cr
#' Chao, A., Hwang, W.-H., Chen, Y.-C. and Kuo. C.-Y. (2000). Estimating the number of shared species in two communities. Statistica Sinica 10, 227-246.\cr\cr
#' Chao, A. and Lee, S.-M. (1992). Estimating the number of classes via sample coverage. Journal of the American Statistical Association 87, 210-217.\cr\cr
#' Chao, A., Ma, M-C. and Yang, M. C. K. (1993). Stopping rule and estimation for recapture debugging with unequal detection rates. Biometrika 80, 193-201.\cr\cr
#' Chao, A., Shen, T.-J. and Hwang, W.-H. (2006). Application of Laplace boundary-mode approximations to estimate species and shared species richness. Australian and New Zealand Journal of Statistics 48, 117-128.\cr\cr
#' Lee, S.-M. and Chao, A. (1994). Estimating population size via sample coverage for closed capture-recapture models. Biometrics 50, 88-97.\cr\cr
#' Lennon, J. J., Koleff, P., Greenwood, J. J. D. and Gaston, K. J. (2001). The geographical structure of British bird distributions: diversity, spatial turnover and scale. J. Anim. Ecol., 70, 966-979.\cr\cr
#' Shen, T.-J. (2003). Prediction of Biodiversity, Ph.D. dissertation, National Tsing-Hua University, HsinChu, Taiwan.\cr\cr
#' Shen, T.-J., Chao, A. and Lin, J.-F. (2003). Predicting the number of new species in further taxonomic sampling. Ecology 84, 798-804.
#' @export

SimilarityPair=function(X, datatype = c("abundance","incidence_freq", "incidence_raw"),nboot=200)
{ 
  
  if(datatype=="abundance")
  {
    if(class(X)=="list"){X <- do.call(cbind,X)}
    type="abundance"
    info1 <- c("S.total"=sum(rowSums(X)>0), "n1"=sum(X[,1]), "n2"=sum(X[,2]), 
               "D1"=sum(X[,1]>0), "D2"=sum(X[,2]>0), "D12"=sum(X[,1]>0 & X[,2]>0),
               "nboot"=nboot)
    
    info2 <- c("f[11]"=sum(X[,1]==1 & X[,2]==1), 
               "f[1+]"=sum(X[,1]==1 & X[,2]>0), "f[+1]"=sum(X[,1]>0 & X[,2]==1),
               "f[2+]"=sum(X[,1]==2 & X[,2]>0), "f[+2]"=sum(X[,1]>0 & X[,2]==2),"f[22]"=sum(X[,1]==2 & X[,2]==2))
    ################################################################2016.07.08-(P.L.Lin)
    plus_CI <-function(x){
      if(x[1] >= 1) x[1] <- 1
      if(x[1] <= 0) x[1] <- 0
      c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
    }
    temp <- list()
    weight <- c(sum(X[,1])/(sum(X[,1])+sum(X[,2])), sum(X[,2])/(sum(X[,1])+sum(X[,2])))
    weight <- - sum(weight*log(weight)) / log(2)
    mat <- Jaccard_Sorensen_Abundance_equ(datatype,X[, 1],X[, 2], nboot)[, c(1, 2)]
    mat <- cbind(mat, mat[, 1]-1.96*mat[, 2],  mat[, 1]+1.96*mat[, 2])
    MLE_Jaccard <- mat[1, ]
    Est_Jaccard <- mat[2, ]
    MLE_Sorensen <- mat[3, ]
    Est_Sorensen <- mat[4, ]
    mat2 <-  Two_Horn_equ(X[,1], X[,2], method="all", weight="unequal", nboot = nboot)
    MLE_Ee_Horn <- mat2$mle
    MLE_Ee_Horn <- plus_CI(c(MLE_Ee_Horn[1],MLE_Ee_Horn[2]))
    Est_Ee_Horn <- mat2$est
    MLE_Ee_U12 <- plus_CI(c(weight*MLE_Ee_Horn[1],MLE_Ee_Horn[2]))
    Est_Ee_U12 <- plus_CI(c(weight*Est_Ee_Horn[1],Est_Ee_Horn[2]))
    mat3 <- Two_BC_equ(X[, 1],X[, 2], datatype="abundance", nboot)
    MLE_Ee_Braycurtis <- mat3$mle
    Est_Ee_Braycurtis <- mat3$est
    mat4 <- SimilarityTwo(X,2,nboot,method="unequal weight")
    MLE_Ee_C22 <- mat4$CqN[1, ]
    Est_Ee_C22 <- mat4$CqN[2, ]
    MLE_Ee_U22 <- mat4$UqN[1, ]
    Est_Ee_U22 <- mat4$UqN[2, ]
    mat5 <- Two_Horn_equ(X[,1], X[,2], method="all", weight="equal", nboot = nboot)
    MLE_ew_Horn <- mat5$mle
    Est_ew_Horn <- mat5$est
    mat6 <- SimilarityTwo(X,2,nboot,method="equal weight")
    MLE_ew_C22 <- mat6$CqN[1, ]
    Est_ew_C22 <- mat6$CqN[2, ]
    MLE_ew_U22 <- mat6$UqN[1, ]
    Est_ew_U22 <- mat6$UqN[2, ]
    #MLE_ew_Braycurtis <- plus_CI(MLE_Braycurtis_equ(X[,1],X[,2],w1=0.5))
    #Est_ew_Braycurtis <- plus_CI(KH_Braycurtis_equ(X[,1],X[,2],w1=0.5))
    MLE_ew_ChaoSoresen <- mat[11,]
    Est_ew_ChaoSoresen <- mat[12, ]
    MLE_ew_ChaoJaccard <- mat[9, ]
    Est_ew_ChaoJaccard <- mat[10, ]
    temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
    rownames(temp[[1]]) <- c("C02(q=0,Sorensen)","U02(q=0,Jaccard)") 
    temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22, MLE_ew_ChaoJaccard, MLE_ew_ChaoSoresen)
    rownames(temp[[2]]) <- c("C12=U12(q=1,Horn)","C22(q=2,Morisita)","U22(q=2,Regional overlap)",
                             "ChaoJaccard","ChaoSorensen")  
    temp[[3]] <- rbind(MLE_Ee_Horn, MLE_Ee_U12, MLE_Ee_C22, MLE_Ee_U22, MLE_Ee_Braycurtis)
    rownames(temp[[3]]) <- c("Horn(q=1)","C12=U12(q=1)","C22(Morisita)", "U22(Regional overlap)","Bray-Curtis")  
    temp[[4]] <- rbind(Est_Sorensen, Est_Jaccard)
    rownames(temp[[4]]) <- c("C02(q=0,Sorensen)","U02(q=0,Jaccard)") 
    temp[[5]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22, Est_ew_ChaoJaccard, Est_ew_ChaoSoresen)
    rownames(temp[[5]]) <- c("C12=U12(q=1,Horn)","C22(q=2,Morisita)","U22(q=2,Regional overlap)",
                             "ChaoJaccard","ChaoSorensen")  
    temp[[6]] <- rbind(Est_Ee_Horn, Est_Ee_U12, Est_Ee_C22, Est_Ee_U22, Est_Ee_Braycurtis)
    rownames(temp[[6]]) <- c("Horn(q=1)","C12=U12(q=1)","C22(Morisita)", "U22(Regional overlap)","Bray-Curtis")  
    z <- list("datatype"=type,"info1"=info1, "info2"=info2, "Empirical_incidence"=temp[[1]], "Empirical_ew"=temp[[2]], "Empirical_ee"=temp[[3]],
              "estimated_incidence"=temp[[4]], "estimated_ew"=temp[[5]], "estimated_ee"=temp[[6]]) 
  }      
  ##---------------------------------------------------------------
  if(datatype=="incidence_raw"){
    a <- X[[1]]
    b <- X[[2]]
    t1 <- ncol(a) ; t2 <- ncol(b)
    data1 <- as.integer(rowSums(a)) ; data2 <- as.integer(rowSums(b))
    y1 <- c(t1, data1)
    y2 <- c(t2, data2)
    X <- cbind(y1, y2)
    type <- "incidence_freq"
    X <- as.data.frame(X)
  }
  if(datatype=="incidence_freq") type <- "incidence_freq" 
  if(datatype=="incidence_freq" | type == "incidence_freq")
  { 
    if(class(X)=="list"){X <- do.call(cbind,X)}
    no.assemblage=length(X[1,])
    Y=X[-1,]  
    type="incidence"
    info1 <- c("S.total"=sum(rowSums(Y)>0), "w"=X[1,1], "z"=sum(X[1,2]), "U1"=sum(Y[,1]), "U2"=sum(Y[,2]), 
               "D1"=sum(Y[,1]>0), "D2"=sum(Y[,2]>0), "D12"=sum(Y[,1]>0 & Y[,2]>0),
               "nboot"=nboot)
    
    info2 <- c("Q[11]"=sum(Y[,1]==1 & Y[,2]==1), 
               "Q[1+]"=sum(Y[,1]==1 & Y[,2]>0), "Q[+1]"=sum(Y[,1]>0 & Y[,2]==1),
               "Q[2+]"=sum(Y[,1]==2 & Y[,2]>0), "Q[+2]"=sum(Y[,1]>0 & Y[,2]==2),  "Q[22]"=sum(Y[,1]==2 & Y[,2]==2))
    
    plus_CI <-function(x){
      if(x[1] >= 1) x[1] <- 1
      if(x[1] <= 0) x[1] <- 0
      c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
    }
    temp <- list()
    weight <- c(sum(Y[,1])/(sum(Y[,1])+sum(Y[,2])), sum(Y[,2])/(sum(Y[,1])+sum(Y[,2])))
    weight <- - sum(weight*log(weight)) / log(2)
    mat <- Jaccard_Sorensen_Abundance_equ(datatype="incidence",X[, 1],X[, 2], nboot)[, c(1, 2)]
    mat <- cbind(mat, mat[, 1]-1.96*mat[, 2],  mat[, 1]+1.96*mat[, 2])
    MLE_Jaccard <- mat[1, ]
    Est_Jaccard <- mat[2, ]
    MLE_Sorensen <- mat[3, ]
    Est_Sorensen <- mat[4, ]
    mat2 <-  Two_Horn_equ(X[,1], X[,2], datatype = "incidence", method="all", weight="unequal", nboot)
    MLE_Ee_Horn <- mat2$mle
    MLE_Ee_Horn <- plus_CI(c(MLE_Ee_Horn[1],MLE_Ee_Horn[2]))
    Est_Ee_Horn <- mat2$est
    MLE_Ee_U12 <- plus_CI(c(weight*MLE_Ee_Horn[1],MLE_Ee_Horn[2]))
    Est_Ee_U12 <- plus_CI(c(weight*Est_Ee_Horn[1],Est_Ee_Horn[2]))
    mat3 <- C2N_ee_se_inc(X, nboot)
    MLE_Ee_C22 <- plus_CI(mat3[1,])
    Est_Ee_C22 <- plus_CI(mat3[3,])
    MLE_Ee_U22 <- plus_CI(mat3[2,])
    Est_Ee_U22 <- plus_CI(mat3[4,])
    mat4 <- Two_Horn_equ(X[,1], X[,2], datatype = "incidence", method="all", weight="equal", nboot)
    MLE_ew_Horn <- mat4$mle
    Est_ew_Horn <- mat4$est
    mat5 <- SimilarityTwo(X, 2, nboot, method="equal weight", datatype="incidence")
    MLE_ew_C22 <- mat5$CqN[1, ]
    Est_ew_C22 <- mat5$CqN[2, ]
    MLE_ew_U22 <- mat5$UqN[1, ]
    Est_ew_U22 <- mat5$UqN[2, ]
    MLE_ew_ChaoSoresen <- mat[11,]
    Est_ew_ChaoSoresen <- mat[12, ]
    MLE_ew_ChaoJaccard <- mat[9, ]
    Est_ew_ChaoJaccard <- mat[10, ]
    mat5 <- Two_BC_equ(X[, 1],X[, 2], datatype="incidence", nboot)
    MLE_Ee_Braycurtis <- mat5$mle
    Est_Ee_Braycurtis <- mat5$est
    temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
    rownames(temp[[1]]) <- c("C02(q=0,Sorensen)","U02(q=0,Jaccard)") 
    temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22, MLE_ew_ChaoJaccard, MLE_ew_ChaoSoresen)
    rownames(temp[[2]]) <- c("C12=U12(q=1,Horn)","C22(q=2,Morisita)","U22(q=2,Regional overlap)",
                             "ChaoJaccard","ChaoSorensen")  
    temp[[3]] <- rbind(MLE_Ee_Horn, MLE_Ee_U12, MLE_Ee_C22, MLE_Ee_U22, MLE_Ee_Braycurtis)
    rownames(temp[[3]]) <- c("Horn(q=1)","C12=U12(q=1)","C22(Morisita)", "U22(Regional overlap)","Bray-Curtis")  
    temp[[4]] <- rbind(Est_Sorensen, Est_Jaccard)
    rownames(temp[[4]]) <- c("C02(q=0,Sorensen)","U02(q=0,Jaccard)") 
    temp[[5]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22, Est_ew_ChaoJaccard, Est_ew_ChaoSoresen)
    rownames(temp[[5]]) <- c("C12=U12(q=1,Horn)","C22(q=2,Morisita)","U22(q=2,Regional overlap)",
                             "ChaoJaccard","ChaoSorensen")  
    temp[[6]] <- rbind(Est_Ee_Horn, Est_Ee_U12, Est_Ee_C22, Est_Ee_U22, Est_Ee_Braycurtis)
    rownames(temp[[6]]) <- c("Horn(q=1)","C12=U12(q=1)","C22(Morisita)", "U22(Regional overlap)","Bray-Curtis")
    z <- list("datatype"=datatype,"info1"=info1, "info2"=info2, "Empirical_incidence"=temp[[1]], "Empirical_ew"=temp[[2]], "Empirical_ee"=temp[[3]],
              "estimated_incidence"=temp[[4]], "estimated_ew"=temp[[5]], "estimated_ee"=temp[[6]]) 
    
  }  
  class(z) <- c("spadeTwo")
  return(z)   
}


#
#
###########################################
#' Estimation of multiple-assemblage similarity measure
#'
#' \code{SimilarityMult}: Estimation of the generalized Sorensen, Horn, and Morisita similarity/dissimilarity indices for comparing frequency or abundance data taken from more than two assemblages. \cr\cr
#' Three types of data (format/information) are supported: ("abundance", "incidence_freq", "incidence_raw"):
#' \enumerate{
#'   \item{(1) Individual-based abundance data (datatype="abundance"): Input data for each assemblage/site include sample species abundances in an empirical sample of n individuals ("reference sample"). When there are N assemblages, input data consist of an S by N abundance matrix, or N lists of species abundances.}
#'   \item{Sampling-unit-based incidence data: There are two kinds of input data.
#'   \describe{
#'      \item{(2)}{Incidence-frequency data (datatype="incidence_freq"): input data for each assemblage consist of species sample incidence frequencies (row sums of each incidence matrix). When there are N assemblages, input data consist of an S+1 by N matrix, or N lists of species incidence frequencies. The first entry of each column/list must be the total number of sampling units, followed by the species incidence frequencies.}
#'      \item{(2B)}{Incidence-raw data (datatype="incidence_raw"): for each assemblage, input data for a reference sample consist of a species-by-sampling-unit matrix; when there are N assemblages, input data consist of N lists of matrices, and each matrix is a species-by-sampling-unit matrix.}
#'    }
#'   }
#' }
#' @param X a matrix, data.frame, lists of species abundances/incidences, or lists of incidence frequencies (see data format/information above).\cr
#' @param datatype type of input data, "abundance", "incidence_freq" or "incidence_raw". \cr
#' @param q set diversity order \code{q = 0}, \code{SimilarityMult} computes the estimates of Sorensen index for pairwise assemblages;
#' set diversity order \code{q = 1}, computes the estimates of Horn index for pairwise assemblages;
#' set diversity order \code{q = 2}, computes the estimates of Morisita index for pairwise assemblages.
#' For diversity order q = 0, 1, 2, \code{SimilarityMult}  computes the overlap estimates among all assemblages.
#' @param nboot an integer specifying the number of bootstrap replications.
#' @param method there are two methods: absolute or relative could be choosen for calculating pairwise similarity indices estimates depended on which diversity order q you choose; But there is no difference for q=0. 
#' @return a list of eleven objects: \cr\cr
#' \code{$info} for summarizing data information.\cr\cr 
#' \code{$Empirical_incidence} for showing the estimation of classic richness-based similarity indices by empirical method. \cr\cr
#' \code{$Empirical_ew} for showing the estimation of similarity indices for comparing species relative abundances by empirical method. \cr \cr
#' \code{$Empirical_ee} for showing the estimation of similarity indices for comparing size-weighted species relative abundances and comparing species absolute abundances by empirical method. \cr\cr
#' \code{$estimated_incidence} for showing the estimation of classic richness-based similarity indices by estimated. \cr\cr
#' \code{$estimated_ew} for showing the estimation of similarity indices for comparing species relative abundances by estimated. \cr\cr
#' \code{$estimated_ee} for showing the estimation of similarity indices for comparing size-weighted species relative abundances and comparing species absolute abundances by estimated. \cr\cr
#' \code{$pairwise} and {$similarity.matrix} for showing pairwise similarity indices estimates (with related statistics) depended on which diversity order \code{q} you choose. \cr\cr
#' \code{$method} absolute or relative. \cr\cr
#' \code{$q} for showing which diversity order \code{q} you choose.
#' @examples
#' \dontrun{
#' data("SimilarityMultDataAbu")
#' SimilarityMult(SimilarityMultDataAbu, datatype="abundance" ,q=2, nboot=200,"relative")
#' data("SimilarityMultDataInci")
#' SimilarityMult(SimilarityMultDataInci, datatype="incidence_freq", q=2,nboot=200,"relative")
#' data("SimilarityMultDataInci_raw")
#' SimilarityMult(SimilarityMultDataInci_raw, datatype="incidence_raw", q=2,nboot=200,"relative")
#' }
#' @references
#' Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two- stage probabilistic approach to multiple-community similarity indices. Biometrics, 64, 1178-1186.\cr\cr
#' Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular Ecology, 17, 4015-4026.
#' @export


SimilarityMult=function(X,datatype=c("abundance","incidence_freq", "incidence_raw"),q=2,nboot=200,method="relative")
{ 
  if(datatype=="abundance"){
    if(class(X)=="list"){X <- do.call(cbind,X)}
    type <- "abundance"
    N <- no.community <- ncol(X)
    temp <- c("N"=ncol(X), "S.total"=sum(rowSums(X)>0))
    n <- apply(X,2,sum)
    D <- apply(X,2,function(x)sum(x>0))
    
    if(N > 2){
      temp1 <- temp2 <- rep(0, N*(N-1)/2)
      k <- 1
      for(i in 1:(N-1)){     
        for(j in (i+1):N){
          temp1[k] <- paste('D',i,j,sep="")
          temp2[k] <- sum(X[,i]>0 & X[,j]>0)
          k <- k + 1
        }
      }
    }
    names(temp2) <- temp1
    names(n) <- paste('n',1:N, sep="")
    names(D) <- paste('D',1:N, sep="")
    info <- c(temp, n, D, temp2)
    if(N == 3) info <- c(temp, n, D, temp2, D123=sum(X[,1]>0 & X[,2]>0 & X[,3]>0))
    info <- c(info, nboot=nboot)
    temp <- list()
    plus_CI <-function(x){
      if(x[1] >= 1) x[1] <- 1
      if(x[1] <= 0) x[1] <- 0
      c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
    }
    n <- apply(X = X, MARGIN = 2, FUN = sum)
    weight <- n/sum(n)
    weight <- - sum(weight*log(weight)) / log(N)
    mat <- SimilarityMul(X, 0, nboot, method ="unequal weight")
    MLE_Jaccard <- mat$UqN[1, ]
    Est_Jaccard <- mat$UqN[2, ]
    MLE_Sorensen <- mat$CqN[1, ]
    Est_Sorensen <- mat$CqN[2, ]
    mat2 <- Horn_Multi_equ(X, datatype="abundance", nboot, method=c("unequal"))
    MLE_Ee_Horn <- mat2$mle
    Est_Ee_Horn <- mat2$est
    Est_Ee_U12 <- plus_CI(c(weight*Est_Ee_Horn[1], Est_Ee_Horn[2]))
    MLE_Ee_U12 <- plus_CI(c(weight*MLE_Ee_Horn[1], MLE_Ee_Horn[2]))
    mat3 <- BC_equ(X, datatype="abundance", nboot)
    MLE_Ee_Braycurtis <- mat3$mle
    Est_Ee_Braycurtis <- mat3$est
    mat4 <- SimilarityMul(X,2,nboot,method="unequal weight")
    MLE_Ee_C22 <- mat4$CqN[1, ]
    Est_Ee_C22 <- mat4$CqN[2, ]
    MLE_Ee_U22 <- mat4$UqN[1, ]
    Est_Ee_U22 <- mat4$UqN[2, ]
    mat5 <- Horn_Multi_equ(X, datatype="abundance", nboot, method=c("unequal"))
    MLE_ew_Horn <- mat5$mle
    Est_ew_Horn <- mat5$est
    mat6 <- SimilarityMul(X,2,nboot,method="equal weight")
    MLE_ew_C22 <- mat6$CqN[1, ]
    Est_ew_C22 <- mat6$CqN[2, ]
    MLE_ew_U22 <- mat6$UqN[1, ]
    Est_ew_U22 <- mat6$UqN[2, ]
    temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
    rownames(temp[[1]]) <- c("C0N(q=0,Sorensen)","U0N(q=0,Jaccard)") 
    temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22)
    rownames(temp[[2]]) <- c("C1N=U1N(q=1,Horn)","C2N(q=2,Morisita)","U2N(q=2,Regional overlap)")  
    temp[[3]] <- rbind(MLE_Ee_Horn, MLE_Ee_U12, MLE_Ee_C22, MLE_Ee_U22, MLE_Ee_Braycurtis)
    rownames(temp[[3]]) <- c("Horn(q=1)","C1N=U1N(q=1)","C2N(Morisita)", "U2N(Regional overlap)","Bray-Curtis")  
    temp[[4]] <- rbind(Est_Sorensen, Est_Jaccard)
    rownames(temp[[4]]) <- c("C0N(q=0,Sorensen)","U0N(q=0,Jaccard)") 
    temp[[5]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22)
    rownames(temp[[5]]) <- c("C1N=U1N(q=1,Horn)","C2N(q=2,Morisita)","U2N(q=2,Regional overlap)")  
    temp[[6]] <- rbind(Est_Ee_Horn, Est_Ee_U12, Est_Ee_C22, Est_Ee_U22, Est_Ee_Braycurtis)
    rownames(temp[[6]]) <- c("Horn(q=1)","C1N=U1N(q=1)","C2N(Morisita)", "U2N(Regional overlap)","Bray-Curtis")  
    Cqn_ee=matrix(0,choose(no.community,2),4)
    Uqn_ee=matrix(0,choose(no.community,2),4)
    Cqn_ew=matrix(0,choose(no.community,2),4)
    Uqn_ew=matrix(0,choose(no.community,2),4)
    k=1
    temp_PC <- temp_PD <- rep(0, N*(N-1)/2)
    for(i in 1:(N-1)){  
      for(j in (i+1):N){
        mat <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal effort')
        Cqn_ee[k,] <- mat[1, ]
        Uqn_ee[k,] <- mat[2, ]
        if(method == "absolute" & q == 1){
          mat2 <- matrix(0, nrow = 2, ncol = 4) 
        }else{
          mat2 <- Cq2_est_equ(X[,c(i,j)], q, nboot, method='equal weight') 
        }
        Cqn_ew[k,] <- mat2[1, ]
        Uqn_ew[k,] <- mat2[2, ]
        temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
        temp_PD[k] <- paste("1-C",q,"2(",i,",",j,")", sep="")
        k <- k+1
      }
    }
    colnames(Cqn_ee) <- colnames( Uqn_ee )<- colnames(Cqn_ew) <- colnames( Uqn_ew ) <-c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
    rownames(Cqn_ee) <- rownames( Uqn_ee )<- rownames(Cqn_ew) <- rownames( Uqn_ew ) <- temp_PC
    Cqn_PC <- list("CqN_ee"=Cqn_ee, "UqN_ee"=Uqn_ee, "CqN_ew"=Cqn_ew, "UqN_ew"=Uqn_ew)
    C_SM_1=matrix(1,N,N)
    C_SM_2=matrix(1,N,N)
    C_SM_3=matrix(1,N,N)
    C_SM_4=matrix(1,N,N)
    k <- 1
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        C_SM_1[i,j] <- C_SM_1[j,i] <- Cqn_ee[k,1]
        C_SM_2[i,j] <- C_SM_2[j,i] <- Uqn_ee[k,1]
        C_SM_3[i,j] <- C_SM_3[j,i] <- Cqn_ew[k,1]
        C_SM_4[i,j] <- C_SM_4[j,i] <- Uqn_ew[k,1]
        k <- k+1
      }
    }
    C_SM <- list("CqN_ee"=C_SM_1, "UqN_ee"=C_SM_2, "CqN_ew"=C_SM_3, "UqN_ew"=C_SM_4)
    
    #z <- list("datatype"=type, "info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
    
    z <- list("datatype"=datatype,"info"=info, "Empirical_incidence"=temp[[1]], "Empirical_ew"=temp[[2]], "Empirical_ee"=temp[[3]],
              "estimated_incidence"=temp[[4]], "estimated_ew"=temp[[5]], "estimated_ee"=temp[[6]], "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
  }
  if(datatype == "incidence_raw"){
    t <- unlist(lapply(X, ncol))
    X <- lapply(X, FUN = function(x){
      rowSums(x)
    })
    X <- do.call(cbind, X)
    X <- rbind(t,X)
    X <- as.data.frame(X)
    if(ncol(X) <= 2) stop("Multiple Commumity measures is only for the data which has three community or more")
    type = "incidence_freq"
  }
  if(datatype=="incidence_freq") type <- "incidence_freq" 
  if(datatype=="incidence_freq" | type == "incidence_freq"){
    if(class(X)=="list"){X <- do.call(cbind,X)}
    Y <- X
    X <- X[-1,]
    t <- as.vector(Y[1,])
    N <- no.community <- ncol(X)
    temp <- c("N"=ncol(X), "S.total"=sum(rowSums(X)>0))
    n <- apply(X,2,sum)
    D <- apply(X,2,function(x)sum(x>0))
    if(N > 2){
      temp1 <- temp2 <- rep(0, N*(N-1)/2)
      k <- 1
      for(i in 1:(N-1)){     
        for(j in (i+1):N){
          temp1[k] <- paste('D',i,j,sep="")
          temp2[k] <- sum(X[,i]>0 & X[,j]>0)
          k <- k + 1
        }
      }
    }
    names(temp2) <- temp1
    names(t) <- paste('t',1:N, sep="")
    names(n) <- paste('u',1:N, sep="")
    names(D) <- paste('D',1:N, sep="")
    info <- c(temp, t, n, D, temp2)
    if(N == 3) info <- c(temp, t, n, D, temp2, D123=sum(X[,1]>0 & X[,2]>0 & X[,3]>0))
    info <- unlist(c(info, nboot=nboot))
    temp <- list()
    plus_CI <-function(x){
      if(x[1] >= 1) x[1] <- 1
      if(x[1] <= 0) x[1] <- 0
      c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
    }
    n <- apply(X = X, MARGIN = 2, FUN = sum)
    weight <- n/sum(n)
    weight <- - sum(weight*log(weight)) / log(N)
    mat <- SimilarityMul(Y, 0, nboot, method ="unequal weight", datatype="incidence")
    MLE_Jaccard <- mat$UqN[1, ]
    Est_Jaccard <- mat$UqN[2, ]
    MLE_Sorensen <- mat$CqN[1, ]
    Est_Sorensen <- mat$CqN[2, ]
    mat2 <- Horn_Multi_equ(Y, datatype="incidence", nboot, method=c("unequal"))
    MLE_Ee_Horn <- mat2$mle
    Est_Ee_Horn <- mat2$est
    Est_Ee_U12 <- plus_CI(c(weight*Est_Ee_Horn[1], Est_Ee_Horn[2]))
    MLE_Ee_U12 <- plus_CI(c(weight*MLE_Ee_Horn[1], MLE_Ee_Horn[2]))
    mat3 <- BC_equ(Y, datatype="incidence", nboot)
    MLE_Ee_Braycurtis <- mat3$mle
    Est_Ee_Braycurtis <- mat3$est
    mat4 <- C2N_ee_se_inc(Y, nboot)
    MLE_Ee_C22 <- plus_CI(mat4[1, ])
    Est_Ee_C22 <- plus_CI(mat4[3, ])
    MLE_Ee_U22 <- plus_CI(mat4[2, ])
    Est_Ee_U22 <- plus_CI(mat4[4, ])
    mat5 <- Horn_Multi_equ(Y, datatype="incidence", nboot, method=c("equal"))
    MLE_ew_Horn <- mat5$mle
    Est_ew_Horn <- mat5$est
    mat6 <- SimilarityMul(Y, 2, nboot, datatype = "incidence", method="equal weight")
    MLE_ew_C22 <- mat6$CqN[1, ]
    Est_ew_C22 <- mat6$CqN[2, ]
    MLE_ew_U22 <- mat6$UqN[1, ]
    Est_ew_U22 <- mat6$UqN[2, ]
    temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
    rownames(temp[[1]]) <- c("C0N(q=0,Sorensen)","U0N(q=0,Jaccard)") 
    temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22)
    rownames(temp[[2]]) <- c("C1N=U1N(q=1,Horn)","C2N(q=2,Morisita)","U2N(q=2,Regional overlap)")  
    temp[[3]] <- rbind(MLE_Ee_Horn, MLE_Ee_U12, MLE_Ee_C22, MLE_Ee_U22, MLE_Ee_Braycurtis)
    rownames(temp[[3]]) <- c("Horn(q=1)","C1N=U1N(q=1)","C2N(Morisita)", "U2N(Regional overlap)","Bray-Curtis")  
    temp[[4]] <- rbind(Est_Sorensen, Est_Jaccard)
    rownames(temp[[4]]) <- c("C0N(q=0,Sorensen)","U0N(q=0,Jaccard)") 
    temp[[5]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22)
    rownames(temp[[5]]) <- c("C1N=U1N(q=1,Horn)","C2N(q=2,Morisita)","U2N(q=2,Regional overlap)")  
    temp[[6]] <- rbind(Est_Ee_Horn, Est_Ee_U12, Est_Ee_C22, Est_Ee_U22, Est_Ee_Braycurtis)
    rownames(temp[[6]]) <- c("Horn(q=1)","C1N=U1N(q=1)","C2N(Morisita)", "U2N(Regional overlap)","Bray-Curtis")  
    
    Cqn_ee=matrix(0,choose(no.community,2),4)
    Uqn_ee=matrix(0,choose(no.community,2),4)
    Cqn_ew=matrix(0,choose(no.community,2),4)
    Uqn_ew=matrix(0,choose(no.community,2),4)
    k=1
    temp_PC <- temp_PD <- rep(0, N*(N-1)/2)
    for(i in 1:(N-1)){  
      for(j in (i+1):N){
        mat <- Cq2_est_equ(Y[,c(i,j)], q, nboot, datatype="incidence", method='equal effort')
        Cqn_ee[k,] <- mat[1, ]
        Uqn_ee[k,] <- mat[2, ]
        if(method == "absolute" & q == 1){
          mat2 <- matrix(0, nrow = 2, ncol = 4) 
        }else{
          mat2 <- Cq2_est_equ(Y[,c(i,j)], q, nboot, datatype="incidence",method='equal weight') 
        }
        Cqn_ew[k,] <- mat2[1, ]
        Uqn_ew[k,] <- mat2[2, ]
        temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
        temp_PD[k] <- paste("1-C",q,"2(",i,",",j,")", sep="")
        k <- k+1
      }
    }
    colnames(Cqn_ee) <- colnames( Uqn_ee )<- colnames(Cqn_ew) <- colnames( Uqn_ew ) <-c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
    rownames(Cqn_ee) <- rownames( Uqn_ee )<- rownames(Cqn_ew) <- rownames( Uqn_ew ) <- temp_PC
    Cqn_PC <- list("CqN_ee"=Cqn_ee, "UqN_ee"=Uqn_ee, "CqN_ew"=Cqn_ew, "UqN_ew"=Uqn_ew)
    C_SM_1=matrix(1,N,N)
    C_SM_2=matrix(1,N,N)
    C_SM_3=matrix(1,N,N)
    C_SM_4=matrix(1,N,N)
    k <- 1
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        C_SM_1[i,j] <- C_SM_1[j,i] <- Cqn_ee[k,1]
        C_SM_2[i,j] <- C_SM_2[j,i] <- Uqn_ee[k,1]
        C_SM_3[i,j] <- C_SM_3[j,i] <- Cqn_ew[k,1]
        C_SM_4[i,j] <- C_SM_4[j,i] <- Uqn_ew[k,1]
        k <- k+1
      }
    }
    C_SM <- list("CqN_ee"=C_SM_1, "UqN_ee"=C_SM_2, "CqN_ew"=C_SM_3, "UqN_ew"=C_SM_4)
    
    z <- list("datatype"=datatype,"info"=info, "Empirical_incidence"=temp[[1]], "Empirical_ew"=temp[[2]], "Empirical_ee"=temp[[3]],
              "estimated_incidence"=temp[[4]], "estimated_ew"=temp[[5]], "estimated_ee"=temp[[6]], "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
  }
  class(z) <- c("spadeMult")
  z
}


#
#
###########################################
#' Estimation of Genetic measure
#'
#' \code{Genetic}: Estimation of the multiple community similarity C_qN and dissimilarity 1-C_qN including Jost (2008) differentiation measure D. \cr\cr
#' One types of data (format/information) is supported: 
#' \enumerate{
#'   \item{Individual-based abundance data \code{(datatype="abundance")}: Input data for each assemblage/site include sample species abundances in an empirical sample of n individuals ("reference sample"). When there are N assemblages, input data consist of an S by N abundance matrix, or N lists of species abundances.}
#' }
#' @param X a matrix, data.frame, lists of species abundances (see data format/information above).\cr
#' @param q set diversity order \code{q = 0}, \code{Genetic} computes the estimates of Sorensen index for pairwise assemblages.
#' set diversity order \code{q = 1}, computes the estimates of Horn index for pairwise assemblages;
#' set diversity order \code{q = 2}, computes the estimates of Morisita index for pairwise assemblages.
#' For diversity order q = 0, 1, 2, \code{Genetic}  computes the overlap estimates among all assemblages.
#' @param nboot an integer specifying the number of bootstrap replications.
#' @return a list of ten objects: \cr\cr
#' \code{$info} for summarizing data information.\cr\cr 
#' \code{$Empirical_incidence} for showing the estimation of classic richness-based dis-similarity indices by empirical method. \cr\cr
#' \code{$Empirical_ew} for showing the estimation of dis-similarity indices for comparing species relative abundances by empirical method. \cr \cr
#' \code{$Empirical_ee} for showing the estimation of dis-similarity indices for comparing size-weighted species relative abundances and comparing species absolute abundances by empirical method. \cr\cr
#' \code{$estimated_incidence} for showing the estimation of classic richness-based dis-similarity indices by estimated. \cr\cr
#' \code{$estimated_ew} for showing the estimation of dis-similarity indices for comparing species relative abundances by estimated. \cr\cr
#' \code{$estimated_ee} for showing the estimation of dis-similarity indices for comparing size-weighted species relative abundances and comparing species absolute abundances by estimated. \cr\cr
#' \code{$pairwise} and {$similarity.matrix} for showing pairwise dis-similarity indices estimates (with related statistics) depended on which diversity order \code{q} you choose. \cr\cr
#' \code{$q} for showing which diversity order \code{q} you choose.
#' @examples
#' data(GeneticDataAbu)
#' Genetic(GeneticDataAbu, q=2, nboot=200)
#' @references
#' Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two-stage probabilistic approach to multiple-community similarity indices. Biometrics, 64, 1178-1186.\cr\cr
#' Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular Ecology, 17, 4015-4026.
#' @export


Genetic=function(X,q=2,nboot=200)
{ 
  if(class(X)=="list"){X <- do.call(cbind,X)}
  type <- "abundance"
  N <- no.community <- ncol(X)
  temp <- c("N"=ncol(X), "S.total"=sum(rowSums(X)>0))
  n <- apply(X,2,sum)
  D <- apply(X,2,function(x)sum(x>0))
  
  if(N > 2){
    temp1 <- temp2 <- rep(0, N*(N-1)/2)
    k <- 1
    for(i in 1:(N-1)){     
      for(j in (i+1):N){
        temp1[k] <- paste('D',i,j,sep="")
        temp2[k] <- sum(X[,i]>0 & X[,j]>0)
        k <- k + 1
      }
    }
  }
  names(temp2) <- temp1
  names(n) <- paste('n',1:N, sep="")
  names(D) <- paste('D',1:N, sep="")
  info <- c(temp, n, D, temp2)
  if(N == 3) info <- c(temp, n, D, temp2, D123=sum(X[,1]>0 & X[,2]>0 & X[,3]>0))
  info <- c(info, nboot=nboot)
  temp <- list()
  n <- apply(X = X, MARGIN = 2, FUN = sum)
  weight <- n/sum(n)
  weight <- - sum(weight*log(weight)) / log(N)
  plus_CI <-function(x){
    if(x[1] >= 1) x[1] <- 1
    if(x[1] <= 0) x[1] <- 0
    c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
  }
  mat2 <- GST_se_equ(X,nboot)
  MLE_ew_Gst <- mat2[1, ]
  Est_ew_Gst <- mat2[2, ]
  mat <- SimilarityMul(X,0,nboot,method="unequal weight")
  MLE_Jaccard <- plus_CI(c(1-mat$UqN[1, 1],mat$UqN[1, 2]))
  Est_Jaccard <- plus_CI(c(1-mat$UqN[2, 1],mat$UqN[2, 2]))
  MLE_Sorensen <- plus_CI(c(1-mat$CqN[1, 1],mat$CqN[1, 2]))
  Est_Sorensen <- plus_CI(c(1-mat$CqN[2, 1],mat$CqN[2, 2]))
  mat3 <- Horn_Multi_equ(X, datatype="abundance", nboot, method=c("unequal"))
  MLE_Ee_Horn <- mat3$mle
  Est_Ee_Horn <- mat3$est
  mat4 <- SimilarityMul(X,2,nboot,method="equal weight")
  mat5 <- Horn_Multi_equ(X, datatype="abundance", nboot, method=c("unequal"))
  MLE_ew_Horn <- mat5$mle
  Est_ew_Horn <- mat5$est
  MLE_ew_Horn <- plus_CI(c(1-MLE_ew_Horn[1],MLE_ew_Horn[2]))
  Est_ew_Horn <- plus_CI(c(1-Est_ew_Horn[1],Est_ew_Horn[2]))
  MLE_ew_C22 <- plus_CI(c(1-mat4$CqN[1, 1],mat4$CqN[1, 2]))
  Est_ew_C22 <- plus_CI(c(1-mat4$CqN[2, 1],mat4$CqN[2, 2]))
  MLE_ew_U22 <- plus_CI(c(1-mat4$UqN[1, 1],mat4$UqN[1, 2]))
  Est_ew_U22 <- plus_CI(c(1-mat4$UqN[2, 1],mat4$UqN[2, 2]))
  temp[[1]] <- rbind(MLE_Sorensen, MLE_Jaccard)
  rownames(temp[[1]]) <- c("1-C0N(q=0,Sorensen)","1-U0N(q=0,Jaccard)") 
  temp[[2]] <- rbind(MLE_ew_Horn, MLE_ew_C22, MLE_ew_U22,MLE_ew_Gst)
  rownames(temp[[2]]) <- c("1-C1N=1-U1N(q=1,Horn)","1-C2N(q=2,Morisita)","1-U2N(q=2,Regional overlap)","Gst")  
  temp[[3]] <- t(as.matrix(MLE_Ee_Horn))
  rownames(temp[[3]]) <- c("Horn(q=1)")  
  temp[[4]] <- rbind(Est_Sorensen, Est_Jaccard)
  rownames(temp[[4]]) <- c("1-C0N(q=0,Sorensen)","1-U0N(q=0,Jaccard)") 
  temp[[5]] <- rbind(Est_ew_Horn, Est_ew_C22, Est_ew_U22, Est_ew_Gst)
  rownames(temp[[5]]) <- c("1-C1N=1-U1N(q=1,Horn)","1-C2N(q=2,Morisita)","1-U2N(q=2,Regional overlap)","Gst")  
  temp[[6]] <- t(as.matrix(Est_Ee_Horn))
  rownames(temp[[6]]) <- c("Horn(q=1)")  

  Cqn_ee=matrix(0,choose(no.community,2),4)
  Cqn_ew=matrix(0,choose(no.community,2),4)
  Uqn_ew=matrix(0,choose(no.community,2),4)
  k=1
  temp_PC <- temp_PD <- rep(0, N*(N-1)/2)
  for(i in 1:(N-1)){  
    for(j in (i+1):N){
      if(q==0 & sum( X[,i]>0 & X[,j]>0)==0){
        mat <- rbind(c(0, 0), c(0 ,0))
        mat2 <- rbind(c(0, 0), c(0 ,0))
      }else{
        mat <- Cq2_est_equ(X[,c(i,j)], q, nboot,method='equal effort')
        mat2 <- Cq2_est_equ(X[,c(i,j)], q, nboot,method='equal weight')
      } 
      Cqn_ee[k,] <- plus_CI(c(mat[1, 1],mat[1, 2]))
      Cqn_ew[k,] <- plus_CI(c(1-mat2[1, 1],mat2[1, 2]))
      Uqn_ew[k,] <- plus_CI(c(1-mat2[2, 1],mat2[2, 2]))
      temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
      temp_PD[k] <- paste("1-C",q,"2(",i,",",j,")", sep="")
      k <- k+1
    }
  }
  colnames(Cqn_ee) <- colnames(Cqn_ew) <- colnames( Uqn_ew ) <-c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
  rownames(Cqn_ee) <- rownames(Cqn_ew) <- rownames( Uqn_ew ) <- temp_PD
  Cqn_PC <- list("CqN_ee"=Cqn_ee, "CqN_ew"=Cqn_ew, "UqN_ew"=Uqn_ew)
  C_SM_1=matrix(1,N,N)
  C_SM_3=matrix(1,N,N)
  C_SM_4=matrix(1,N,N)
  k <- 1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      C_SM_1[i,j] <- C_SM_1[j,i] <- Cqn_ee[k,1]
      C_SM_3[i,j] <- C_SM_3[j,i] <- Cqn_ew[k,1]
      C_SM_4[i,j] <- C_SM_4[j,i] <- Uqn_ew[k,1]
      k <- k+1
    }
  }
  C_SM <- list("CqN_ee"=C_SM_1, "CqN_ew"=C_SM_3, "UqN_ew"=C_SM_4)
  
  z <- list("info"=info, "Empirical_incidence"=temp[[1]], "Empirical_ew"=temp[[2]], "Empirical_ee"=temp[[3]],
            "estimated_incidence"=temp[[4]], "estimated_ew"=temp[[5]], "estimated_ee"=temp[[6]], "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "q"=q)
  
  class(z) <- c("spadeGenetic")
  z
}