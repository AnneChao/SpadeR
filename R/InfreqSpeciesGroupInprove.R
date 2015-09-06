InfreqSpeciesGroup <- function(data, k){
  data <- data[-1]
  data <- as.numeric(data)
  Q <- function(i, data){length(data[which(data == i)])}
  
  x <- data[which(data != 0)]
  r <- c(1:k)
  Rare.Species.Group <- matrix(sapply(r, function(r)Q(r, x)), 1, k)
  rownames(Rare.Species.Group) <- c("Infreq.Species.Group")
  colnames(Rare.Species.Group) <- paste("Q", r, sep="")
  return(Rare.Species.Group)
}
