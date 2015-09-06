RareSpeciesGroup <-
function(data, k){
  data <- as.numeric(data)
  f <- function(i, data){length(data[which(data == i)])}

  x <- data[which(data != 0)]
  r <- c(1:k)
  Rare.Species.Group <- matrix(sapply(r, function(r)f(r, x)), 1, k)
  rownames(Rare.Species.Group) <- c("f(i)")
  colnames(Rare.Species.Group) <- c(1:k)
  return(Rare.Species.Group)
}
