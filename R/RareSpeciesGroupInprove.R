RareSpeciesGroup <- function(data, k){
  if (is.matrix(data) == T || is.data.frame(data) == T){
    if (ncol(data) != 1 & nrow(data) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(data) == 1){
      data <- data[, 1]
    } else {
      data <- data[1, ]
    }
  }
  
  data <- as.numeric(data)
  f <- function(i, data){length(data[which(data == i)])}
  
  x <- data[which(data != 0)]
  r <- c(1:k)
  Rare.Species.Group <- matrix(sapply(r, function(r)f(r, x)), 1, k)
  rownames(Rare.Species.Group) <- c("Rare.Species.Group")
  colnames(Rare.Species.Group) <- paste("f", r, sep="")
  return(Rare.Species.Group)
}
