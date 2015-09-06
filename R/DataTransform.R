DataTransform <-
function(data, type = c("FreqCount", "MatrixInci", "MatrixAbun", "InciCount")){
  Freq2Abun <- function(data){
    dat <- data[-(1:2)]  
    j <- 1:length(dat)
    data.abun <- as.numeric(rep(dat[which(j %% 2 == 1)], dat[which(j %% 2 == 0)]))
    return(data.abun)
  }
  Mat2Inci <- function(data){
    t <- length(data[1,])
    dat <- apply(data, 1, sum)
    data.inci <- c(t, dat)
    return(data.inci)  
  }
  Count2Inci <- function(data){
    t <- data[1]
    dat <- data[-(1:3)]
    j <- 1:length(dat)
    data.inci <- as.numeric(c(t, as.numeric(rep(dat[which(j %% 2 == 1)], dat[which(j %% 2 == 0)]))))
    return(data.inci)
  }
  
  if (type == "FreqCount"){
    data <- Freq2Abun(data)
  } else if (type == "MatrixInci"){
    data <- Mat2Inci(data)
  } else if (type == "MatrixAbun") {
    data <- apply(data, 1, sum)
  } else {
    data <- Count2Inci(data)
  }
  return(data)
}
