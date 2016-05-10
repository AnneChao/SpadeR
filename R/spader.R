#
#
###########################################
#' Estimation of species richness 
#' 
#' \code{ChaoSpecies}: estimation of species richness in one assemblage
#' 
#' @param data a column vector of species abundances or incidence frequencies.\cr 
#' If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units, 
#' followed by species incidence frequencies. (See example \code{data(ChaoSpeciesDataInci)}).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or 
#' sampling-unit-based incidence data (\code{datatype = "incidence"}).
#' @param k cut-off point; it is a value that separates frequency counts into abundant and rare groups.
#' @param conf a positive number \eqn{\le} 1 specifying the level of confidence interval.
#' @return  a list of three objects: 
#' \code{$Basic.Data.Information} and \code{$Rare.Species.Group}/\code{$Infreq.Species.Group} 
#' for summarizing data information; \code{$Species.Table} for showing a table of various species richness estimates, standard errors, 
#' and the associated confidence intervals.
#' @examples
#' data(ChaoSpeciesDataAbu)
#' ChaoSpecies(ChaoSpeciesDataAbu, datatype="abundance", k = 10, conf = 0.95)
#' data(ChaoSpeciesDataInci)
#' ChaoSpecies(ChaoSpeciesDataInci, datatype="incidence", k = 10, conf = 0.95)
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


ChaoSpecies <- function(data, datatype = c("abundance", "incidence"), k = 10, conf = 0.95){
    method <- "all"
    if (k != round(k) || k < 0) 
      stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
    if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) 
      stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
    
    if (is.matrix(data) == T || is.data.frame(data) == T){
      if (ncol(data) != 1 & nrow(data) != 1)
        stop("Error: The data format is wrong.")
      if (ncol(data) == 1){
        data <- data[, 1]
      } else {
        data <- data[1, ]
      }
    }
    #  data <- as.numeric(round(data))
    
    if (datatype == "abundance"){
      f <- function(i, data){length(data[which(data == i)])}
      if (f(1, data) == sum(data)){
        stop("Error: The information of data is not enough.")}
      z <- (list(Basic.Data.Information = basicAbun(data, k)[[1]], Rare.Species.Group = RareSpeciesGroup(data, k), 
                 Species.Table = round(SpecAbunOut(data, method, k, conf), 3)))
    } else {
      dat <- data[-1]; 
      Q <- function(i, data){length(data[which(data == i)])}
      if (Q(1, dat) == sum(dat)){
        stop("Error: The information of data is not enough.")}
      z <- (list(Basic.Data.Information = basicInci(data, k)[[1]], Infreq.Species.Group = InfreqSpeciesGroup(data, k),
                 Species.Table = round(SpecInciOut(data, method, k, conf),3)))
    }
    class(z) <- c("ChaoSpecies")
    z
  }


#
#
###########################################
#' Estimation of the number of shared species in two assemblages
#' 
#' \code{ChaoShared}: estimation of shared species richness in two assemblages based on the following two sampling schemes. \cr
#' \enumerate{
#'   \item{(Abundance data) In each assemblage, a random sample of individuals is taken and species abundances/frequencies are recorded.}
#'   \item{(Incidence data) Each assemblage is sampled several times and species presence/absence data for multiple sampling units are recorded.}
#' }
#' @param data a matrix, data.frame (species by sites) of species abundances or incidence frequencies with two columns. 
#' If \code{datatype = "incidence"}, then the first entry of the input data must be the total number of sampling units, 
#' followed by species incidence frequencies in each column (See example \code{data(ChaoSharedDataInci)}).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or 
#' sampling-unit-based incidence data (\code{datatype = "incidence"}).
#' @param se a logical variable to calculate the bootstrap standard error and the associated confidence interval.
#' @param nboot an integer specifying the number of bootstrap replications.
#' @param conf a positive number \eqn{\le} 1 specifying the level of confidence interval.
#' @return a list of two objects: \code{$BASIC_DATA_INFORMATION} for summarizing data information;\cr
#' \code{$ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES} for showing a table of various shared richess estimates, standard errors, and the associated confidence intervals.
#' @examples
#' data(ChaoSharedDataAbu)
#' ChaoShared(ChaoSharedDataAbu, datatype="abundance", se=TRUE, nboot=200, conf=0.95)
#' data(ChaoSharedDataInci)
#' ChaoShared(ChaoSharedDataInci, datatype="incidence", se=TRUE, nboot=200, conf=0.95)
#' @references 
#' Chao, A., Hwang, W.-H., Chen, Y.-C. and Kuo. C.-Y. (2000). Estimating the  number of shared species in two communities. Statistica Sinica 10, 227-246.\cr\cr
#' Chao, A., Shen, T.-J. and Hwang, W.-H. (2006). Application of Laplace boundary-mode approximations to estimate species and shared species richness.  Australian and New Zealand Journal of Statistics 48, 117-128.\cr\cr
#' Pan, H. Y., Chao, A. and Foissner, W. (2009). A non-parametric lower bound for the number of species shared by multiple communities. Journal of Agricultural, Biological and Environmental Statistics 14, 452-468.
#' @export

ChaoShared <-
  function(data, datatype = c("abundance", "incidence"), 
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
      x1 <- data[, 1]
      x2 <- data[, 2]
      Basic <- BasicFun(x1, x2, nboot, datatype)
      #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
      output <- ChaoShared.Ind(x1, x2, method, nboot, conf, se)
      colnames(output) <- c("Estimate", "s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
    }
    if (datatype == "incidence") {
      y1 <- data[, 1]
      y2 <- data[, 2]
      Basic <- BasicFun(y1, y2, B=nboot, datatype)
      #     cat("(2)  ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES: ", "\n")
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
#' \code{Diversity}: estimation of various diversity indices including species richness (diversity of order 0), 
#' the Shannon index/diversity (diversity of order 1), the Simpson index/diversity (diversity of order 2), 
#' and Hill number (diversity of order from 0 to 3) .
#' @param data a column vector of species abundances or incidence frequencies.\cr 
#' If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units, 
#' followed by species incidence frequencies. (See example \code{data(ChaoSpeciesDataInci)}).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or 
#' sampling-unit-based incidence data (\code{datatype = "incidence"}).
#' @param q a vector of nonnegative numbers for which the diversity order of Hill numbers will be estimated. If \code{NULL}, then
#' Hill numbers will be estimated at order q from 0 to 3 with equally-spaced 0.25. 
#' @return a list of seven objects: 
#' \code{$BASIC.DATA} for summarizing data information; \code{$SPECIES.RICHNESS} for showing various species richness estimates along with related statistics; 
#' \code{$SHANNON.INDEX} and \code{$EXPONENTIAL.OF.SHANNON.INDEX} for showing various Shannon index/diversity estimates; 
#' \code{$SIMPSON.INDEX} and \code{$INVERSE.OF.SIMPSON.INDEX} for showing various Simpson index/diversity estimates; 
#' \code{$HILL.NUMBERS} for showing Hill number (diversity of order from 0 to 3) estimates. \cr\cr
#' @examples
#' \dontrun{
#' data(DiversityDataAbu)
#' Diversity(DiversityDataAbu,datatype="abundance")
#' data(SimilarityPairDataInci)
#' Diversity(SimilarityPairDataInci[,1],datatype="incidence")
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

Diversity=function(data, datatype=c("abundance","incidence"), q=NULL)
{
  X <- data
  if(datatype=="abundance"){
    type="abundance"
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))
    
    BASIC.DATA <- matrix(round(c(sum(X), sum(X>0), 1-sum(X==1)/sum(X), CV.Ind(X)),3), ncol = 1)
    nickname <- matrix(c("n", "D", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    
    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("Number of observed individuals", "Number of observed species",
                              "Estimated sample coverage",
                              "Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    
    table0 <- matrix(0,5,4)
    table0[1,]=c(Chao1(X)[-5])
    table0[2,]=c(Chao1_bc(X))
    table0[3,]=round(SpecAbuniChao1(X, k=10, conf=0.95)[1,],1)
    table0[4,]=round(c(SpecAbunAce(X)),1)
    table0[5,]=round(c(SpecAbunAce1(X)),1)
    colnames(table0) <- c("Estimate", "s.e.", paste(Chao1(X)[5]*100,"%Lower", sep=""), paste(Chao1(X)[5]*100,"%Upper", sep=""))
    rownames(table0) <- c(" Chao1 (Chao, 1984)"," Chao1-bc ", " iChao1"," ACE (Chao & Lee, 1992)",
                          " ACE-1 (Chao & Lee, 1992)")
    
    SHANNON=Shannon_index(X)
    table1=round(SHANNON[c(1:5),],3)
    table1=table1[-2,]              ##2016.05.09
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Jackknife",
    #                      " Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c(" MLE"," Jackknife",
                      " Chao & Shen"," Chao et al. (2013)")
    
    table1_exp=round(SHANNON[c(6:10),],3)
    table1_exp=table1_exp[-2,]      ##2016.05.09
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Jackknife",
     #                         " Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c(" MLE"," Jackknife",
                              " Chao & Shen"," Chao et al. (2013)")
    
    table2=round(Simpson_index(X)[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c(" MVUE"," MLE")
    
    table2_recip=round(Simpson_index(X)[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c(" MVUE"," MLE")
    
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
    
    z <- list("datatype"= type,"BASIC.DATA"=BASIC.DATA,"SPECIES.RICHNESS"=table0, 
              "SHANNON.INDEX"=table1,"EXPONENTIAL.OF.SHANNON.INDEX"=table1_exp,
              "SIMPSON.INDEX"=table2,"INVERSE.OF.SIMPSON.INDEX"=table2_recip,
              "HILL.NUMBERS"= Hill)
  }else if(datatype=="incidence"){
    if(!is.vector(X)) X <- as.numeric(unlist(c(X)))
    type="incidence"
    U<-sum(X[-1])
    D<-sum(X[-1]>0)
    T<-X[1]
    C<-Chat.Sam(X,T)
    CV_squre<-max( D/C*T/(T-1)*sum(X[-1]*(X[-1]-1))/U^2-1, 0)
    CV<-CV_squre^0.5
    BASIC.DATA <- matrix(round(c(D,T, C, CV),3), ncol = 1)
    nickname <- matrix(c("D", "T", "C", "CV"), ncol = 1)
    BASIC.DATA <- cbind(nickname, BASIC.DATA)
    
    colnames(BASIC.DATA) <- c("Variable", "Value")
    rownames(BASIC.DATA) <- c("Number of observed species", "Number of Sampling units",
                              "Estimated sample coverage",
                              "Estimated CV")
    BASIC.DATA <- data.frame(BASIC.DATA)
    #BASIC.DATA <- basicInci(X, k=10)[[1]]
    ############################################################
    table0=SpecInci(X, k=10, conf=0.95)
    SHANNON=Shannon_Inci_index(X)
    table1=round(SHANNON[c(1,4),],3)
    colnames(table1) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1) <- c(" MLE"," Chao et al. (2013)")
    table1_exp=round(SHANNON[c(5,8),],3)
    colnames(table1_exp) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    #rownames(table1_exp) <- c(" MLE"," MLE_bc"," Chao & Shen"," Chao et al. (2013)")
    rownames(table1_exp) <- c(" MLE"," Chao et al. (2013)")
    
    SIMPSON=Simpson_Inci_index(X)
    table2=round(SIMPSON[c(1:2),],5)
    colnames(table2) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2) <- c(" MVUE"," MLE")
    
    table2_recip=round(SIMPSON[c(3:4),],5)
    colnames(table2_recip) <- c("Estimate", "s.e.", paste("95%Lower"), paste("95%Upper"))
    rownames(table2_recip) <- c(" MVUE"," MLE")
    
    
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
#' \code{SimilarityPair}: estimation of various similarity indices based on either individual-based abundance data or sampling-unit-based incidence data 
#' taken from two assemblages. The incidence-based indices include the classic Jaccard, Sorensen and Lennon et al. (2001) indices, 
#' and the abundance-based indices include the Bray-Curtis, Morisita-Horn, Horn and abundance-based Jaccard 
#' and Sorensen indices developed in Chao et al. (2005).
#' 
#' @param X a matrix (species by sites) of species abundances, or incidence frequencies with two columns.
#' If \code{datatype = "incidence"}, then the first entry of the input data must be total number of sampling units in each column.
#' (See example \code{data(SimilarityPairDataInci)}).
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or 
#' sampling-unit-based incidence data (\code{datatype = "incidence"}).
#' @param nboot an integer specifying the number of replications.
#' @return a list of six objects: \code{$datatype} for showing data type; 
#' \code{$info1} and \code{$info2} for summarizing data information;
#' \code{$similarity} for showing the various similarity index estimates along with related statistics;
#' \code{$assemblage1} and {$assemblage2} for showing species richness estimates with related statistics of each assemblage.
#' @examples
#' data(SimilarityPairDataAbu)
#' SimilarityPair(SimilarityPairDataAbu, datatype="abundance")
#' data(SimilarityPairDataInci)
#' SimilarityPair(SimilarityPairDataInci,datatype="incidence")
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

SimilarityPair=function(X, datatype = c("abundance","incidence"),nboot=200)
{
  if(datatype=="abundance")
  {
    type="abundance"
    info1 <- c("S.total"=sum(rowSums(X)>0), "n1"=sum(X[,1]), "n2"=sum(X[,2]), 
               "D1"=sum(X[,1]>0), "D2"=sum(X[,2]>0), "D12"=sum(X[,1]>0 & X[,2]>0),
               "nboot"=nboot)
    
    info2 <- c("f[11]"=sum(X[,1]==1 & X[,2]==1), 
               "f[1+]"=sum(X[,1]==1 & X[,2]>0), "f[+1]"=sum(X[,1]>0 & X[,2]==1),
               "f[2+]"=sum(X[,1]==2 & X[,2]>0), "f[+2]"=sum(X[,1]>0 & X[,2]==2))
    
    temp <- list()
    temp[[1]] <- Jaccard_Sorensen_Abundance_equ(datatype,X[,1],X[,2],nboot)
    temp[[1]] <- matrix(sapply(c(temp[[1]]), function(x) ifelse(x==0,"",x)),12)
    subset=round(C1n_equ(method="absolute",X[,c(1,2)],nboot),4)[1:2]
    subset[1]=1- subset[1]
    #X1=as.numeric(X[,1])
    #X2=as.numeric(X[,2])
    #subset=New_C12(X1, X2)
    #subset=temp[[1]][13,c(1,2)]
    Horn_MLE=Two_horn_MLE_equ(X[,1],X[,2])
    BC_estimated=KH_Braycurtis_equ(X[,1],X[,2],w1=sum(X[,1])/(sum(X[,1])+sum(X[,2])))
    MLE_ew_Braycurtis=MLE_Braycurtis_equ(X[,1],X[,2],w1=0.5)
    BC_ew_estimated=KH_Braycurtis_equ(X[,1],X[,2],w1=0.5)
    temp[[1]] <- rbind(temp[[1]][1:8,],
                       #c(round(C1n_equ(method="relative",X[,c(1,2)],nboot),4)[1:2],"","","",""),
                       NA,
                       c(subset,"","","",""),
                       temp[[1]][9:12,],c(Horn_MLE,"","","",""),c(BC_estimated,"","","",""),
                       c(MLE_ew_Braycurtis,"","","",""),c(BC_ew_estimated,"","","",""))
    
    colnames(temp[[1]]) <- c("Estimate", "s.e.", "U_hat*", "U_hat* se.", "V_hat**", "V_hat** se.")
    #rownames(temp[[1]]) <- c("Jaccard incidence", "Sorensen incidence", "Lennon et al. (2001)",
    #                         "Bray-Curtis", "Morisita-Horn", "Morisita Original",
    #                         "Horn (relative)", "Horn (absolute)","Jaccard Abundance (unadjusted)",
     #                        "Jaccard Abundance (adjusted)", "Sorensen Abundance(unadjusted)", "Sorensen Abundance(adjusted)")
    rownames(temp[[1]]) <- c("Jaccard incidence (observed)","Jaccard incidence (estimated)", "Sorensen incidence (observed)","Sorensen incidence (estimated)", "Lennon et al. (2001)",
                             "Bray-Curtis (observed)", "Morisita-Horn", "Morisita Original","Horn (relative)", "Horn (estimated)",
                             "Jaccard Abundance (unadjusted)",
                           "Jaccard Abundance (adjusted)", "Sorensen Abundance(unadjusted)", "Sorensen Abundance(adjusted)","Horn (observed)","Bray-Curtis (estimated)","Bray-curtis (observed-ew)","Bray-Curtis (estimated-ew)")
    
    
    temp[[1]] <- temp[[1]][-9,]
    temp[[1]] <- as.data.frame(temp[[1]])
    
    #temp[[2]] <- rbind(Chao1_equ(X[,1],conf=0.95)[-5],Chao1_bc_equ(X[,1],conf=0.95),
     #                  SpecAbunAce_equ(X[,1], k=10, conf=0.95), SpecAbunAce1_equ(X[,1] ,k=10, conf=0.95))
    #colnames(temp[[2]]) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
    #rownames(temp[[2]]) <- c("Chao1", "Chao1-bc", "ACE", "ACE-1")
    #temp[[2]] <- as.data.frame(temp[[2]])
    
    #temp[[3]] <- rbind(Chao1_equ(X[,2],conf=0.95)[-5],Chao1_bc_equ(X[,2],conf=0.95),
     #                  SpecAbunAce_equ(X[,2], k=10, conf=0.95), SpecAbunAce1_equ(X[,2] ,k=10, conf=0.95))
    #colnames(temp[[3]]) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
    #rownames(temp[[3]]) <- c("Chao1", "Chao1-bc", "ACE", "ACE-1")
    #temp[[3]] <- as.data.frame(temp[[3]])
    #z <- list("datatype"=type,"info1"=info1, "info2"=info2, "similarity"=temp[[1]], "assemblage1"=temp[[2]], "assemblage2"=temp[[3]])
    z <- list("datatype"=type,"info1"=info1, "info2"=info2, "similarity"=temp[[1]]) 
  }      
  ##---------------------------------------------------------------
  if(datatype=="incidence")
  {
    no.assemblage=length(X[1,])
    Y=X[-1,]  
    "type"="incidence"
    info1 <- c("S.total"=sum(rowSums(Y)>0), "w"=X[1,1], "z"=sum(X[1,2]), 
               "D1"=sum(Y[,1]>0), "D2"=sum(Y[,2]>0), "D12"=sum(Y[,1]>0 & Y[,2]>0),
               "nboot"=nboot)
    
    info2 <- c("Q[11]"=sum(Y[,1]==1 & Y[,2]==1), 
               "Q[1+]"=sum(Y[,1]==1 & Y[,2]>0), "Q[+1]"=sum(Y[,1]>0 & Y[,2]==1),
               "Q[2+]"=sum(Y[,1]==2 & Y[,2]>0), "Q[+2]"=sum(Y[,1]>0 & Y[,2]==2))
    
    temp <- list()
    temp[[1]] <- Jaccard_Sorensen_Abundance_equ(datatype,X[,1],X[,2],nboot)
    temp[[1]] <- matrix(sapply(c(temp[[1]]), function(x) ifelse(x==0,"",x)),12)
    colnames(temp[[1]]) <- c("Estimate", "s.e.", "U_hat*", "U_hat* se.", "V_hat**", "V_hat** se.")
    rownames(temp[[1]]) <- c("Jaccard incidence (observed)","Jaccard incidence (estimated)", "Sorensen incidence (observed)","Sorensen incidence (estimated)","Lennon et al. (2001)",
                             "Bray-Curtis", "Morisita-Horn", "Morisita Original",
                             "Incidence-based Jaccard (unadjusted)", "Incidence-based Jaccard (adjusted)", 
                             "Incidence-based Sorensen(unadjusted)", "Incidence-based Sorensen(adjusted)")
    temp[[1]] <- temp[[1]]
    
    temp[[2]] <- rbind(SpecInciChao2(X[,1],conf=0.95)[-5],SpecInciChao2bc(X[,1],conf=0.95),
                       SpecInciModelh(X[,1], k=10, conf=0.95), SpecInciModelh1(X[,1] ,k=10, conf=0.95))
    colnames(temp[[2]]) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
    rownames(temp[[2]]) <- c("Chao2", "Chao2-bc", "ICE", "ICE-1")
    temp[[2]] <- as.data.frame(temp[[2]])
    
    temp[[3]] <- rbind(SpecInciChao2(X[,2],conf=0.95)[-5],SpecInciChao2bc(X[,2],conf=0.95),
                       SpecInciModelh(X[,2], k=10, conf=0.95), SpecInciModelh1(X[,2] ,k=10, conf=0.95))
    colnames(temp[[3]]) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
    rownames(temp[[3]]) <- c("Chao2", "Chao2-bc", "ICE", "ICE-1")
    temp[[3]] <- as.data.frame(temp[[3]])
    z <- list("datatype"=type, "info1"=info1, "info2"=info2, "similarity"=as.data.frame(temp[[1]]), "assemblage1"=temp[[2]], "assemblage2"=temp[[3]])
  }  
  class(z) <- c("spadeTwo")
  return(z)   
}


#
#
###########################################
#' Estimation of multiple-assemblage similarity measure 
#' 
#' \code{SimilarityMult}: estimation of the generalized Sorensen, Horn, and Morisita similarity/dissimilarity indices for comparing frequency or abundance data taken from more than two assemblages.
#' @param X a matrix (species by sites) of species abundances.
#' @param q set diversity order \code{q = 0}, \code{SimilarityMult} computes the estimates of Sorensen index for pairwise assemblages; 
#' set diversity order \code{q = 1}, computes the estimates of Horn index for pairwise assemblages; 
#' set diversity order \code{q = 2}, computes the estimates of Morisita index for pairwise assemblages. 
#' For diversity order q = 0, 1, 2, \code{SimilarityMult}  computes the overlap estimates among all assemblages.
#' @param nboot an integer specifying the number of bootstrap replications.
#' @return a list of six objects: \code{$info} for summarizing data information; \code{$overlap} for showing overlap estimates along with related statistics among all assemblages;
#' \code{$pairwise} and {$similarity.matrix} for showing pairwise similarity estimates (with related statistics) depended on which diversity order \code{q} you choose;
#' \code{$method}: fixed \code{method="absolute"}; \code{$q} for showing which diversity order \code{q} you choose.
#' @examples
#' data(SimilarityMultDataAbu)
#' SimilarityMult(SimilarityMultDataAbu,q=2,nboot=200)
#' @references
#' Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two- stage probabilistic approach to multiple-community similarity indices. Biometrics, 64, 1178-1186.\cr\cr
#' Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular Ecology, 17, 4015-4026.
#' @export


SimilarityMult=function(X,q=2,nboot=200)
{ 
  type <- "abundance"
  method<-"absolute"
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
  
  
  Cqn=rbind(Cqn_se_equ(X,q=0,nboot),
            #C1n_equ(method="relative",X,nboot),
            NA,
            C1n_equ(method="absolute",X,nboot), 
            Cqn_se_equ(X,q=2,nboot)[1:4])
  if(N==3){Cqn <- rbind(Cqn, C33_se_equ(X,nboot)[1:4])}
  colnames(Cqn) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
  rownames(Cqn) <- c(paste("C0",N," (Sorensen)",sep=""),paste("C1",N,"(Horn)",sep=""),paste("C1",N,"*","(Horn)",sep=""),paste("C2",N," (Morisita)",sep=""),if(N==3) "C33")
  Cqn <- Cqn[-2,]
  
  if(q==0 || q==1){Cqn_PC=matrix(0,choose(no.community,2),4)}
  if(q==2)        {Cqn_PC=matrix(0,choose(no.community,2),6)}
  k=1
  temp_PC <- temp_PD <- rep(0, N*(N-1)/2)
  for(i in 1:(N-1)){  
    for(j in (i+1):N){
      Cqn_PC[k,] <- Cqn_se_equ(X[,c(i,j)],q,nboot,method)
      temp_PC[k] <- paste("C22(",i,",",j,")", sep="")
      temp_PD[k] <- paste("1-C22(",i,",",j,")", sep="")
      k <- k+1
    }
  }
  if(q==0 || q==1){
    colnames(Cqn_PC) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL")
    rownames(Cqn_PC) <- temp_PC
  }
  if(q==2){
    colnames(Cqn_PC) <- c("Estimate", "s.e.", "95%.LCL", "95%.UCL", "D.95%.LCL", "D.95%.UCL")
    rownames(Cqn_PC) <- temp_PC
  }
  
  C_SM=matrix(1,N,N)
  k <- 1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      C_SM[i,j] <- C_SM[j,i] <- Cqn_PC[k,1]
      k <- k+1
    }
  }
  
  #z <- list("datatype"=type, "info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
  z <- list("info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
  class(z) <- c("spadeMult")
  z
}


#
#
###########################################
#' Estimation of Genetic measure 
#' 
#' \code{Genetic}: estimation of the multiple community similarity C_qN and dissimilarity 1-C_qN including Jost (2008) differentiation measure D.
#' @param X a matrix (species by sites) of species abundances.
#' @param q set diversity order \code{q = 0}, \code{Genetic} computes the estimates of Sorensen index for pairwise assemblages; 
#' set diversity order \code{q = 1}, computes the estimates of Horn index for pairwise assemblages; 
#' set diversity order \code{q = 2}, computes the estimates of Morisita index for pairwise assemblages. 
#' For diversity order q = 0, 1, 2, \code{Genetic}  computes the overlap estimates among all assemblages.
#' @param nboot an integer specifying the number of bootstrap replications.
#' @return a list of five objects: \code{$info} for summarizing data information;
#' \code{$overlap} for showing the C0N, C1N and C2N of this data;
#' \code{$pairwise} for showing the pairwise data CqN output;
#' \code{$similarity.matrix} for showing pairwise similarity of this data. \code{q} you choose.
#' @examples
#' data(GeneticDataAbu)
#' Genetic(GeneticDataAbu, q=2, nboot=200)
#' @references
#' Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two-stage probabilistic approach to multiple-community similarity indices. Biometrics, 64, 1178-1186.\cr\cr
#' Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular Ecology, 17, 4015-4026.
#' @export


Genetic=function(X,q=2,nboot=200)
{ 
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
  
  
  Cqn=rbind(Cqn_se_equ(X,q=0,nboot),
            #C1n_equ(method="relative",X,boot),
            NA,
            C1n_equ(method="absolute",X,nboot), 
            Cqn_se_equ(X,q=2,nboot)[1:4])
  #if(N==3){Cqn <- rbind(Cqn, C33_se_equ(X,boot))}
  colnames(Cqn) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL")
  rownames(Cqn) <- c(paste("C0",N,sep=""),paste("C1",N,sep=""),paste("C1",N,"*",sep=""),paste("C2",N,sep=""),if(N==3) "C33")
  
  
  if(q==0 || q==1){Cqn_PC=matrix(0,choose(no.community,2),4)}
  if(q==2)        {Cqn_PC=matrix(0,choose(no.community,2),6)}
  k=1
  temp_PC <- temp_PD <- rep(0, N*(N-1)/2)
  for(i in 1:(N-1)){  
    for(j in (i+1):N){
      Cqn_PC[k,] <- Cqn_se_equ(X[,c(i,j)],q,nboot,method="absolute")
      temp_PC[k] <- paste("C",q,"2(",i,",",j,")", sep="")
      temp_PD[k] <- paste("1-C",q,"2(",i,",",j,")", sep="")
      k <- k+1
    }
  }
  if(q==0 || q==1){
    colnames(Cqn_PC) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL")
    rownames(Cqn_PC) <- temp_PC
  }
  if(q==2){
    colnames(Cqn_PC) <- c("Estimate", "Est_s.e.", "95%.LCL", "95%.UCL", "D.95%.LCL", "D.95%.UCL")
    rownames(Cqn_PC) <- temp_PC
  }
  
  C_SM=matrix(1,N,N)
  k <- 1
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      C_SM[i,j] <- C_SM[j,i] <- Cqn_PC[k,1]
      k <- k+1
    }
  }
  #z <- list("datatype"=type, "info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "method"=method, "q"=q)
  z <- list("info"=info, "overlap"=Cqn, "pairwise"=Cqn_PC, "similarity.matrix"=C_SM, "q"=q)
  class(z) <- c("spadeGenetic")
  z
}
