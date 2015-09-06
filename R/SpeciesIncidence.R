SpeciesIncidence <-
function(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"), k = 10, conf = 0.95){
  method <- match.arg(method)
  return(list(BASIC.DATA.INFORMATION = basicInci(data, k)[[1]], SPECIES.TABLE = SpecInciOut(data, method, k, conf)))
}
