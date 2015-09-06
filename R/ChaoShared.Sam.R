ChaoShared.Sam <-
function(y1, y2, method = c("all",
                                              "Chao2-shared", 
                                              "Chao2-shared-bc",
                                              "Lower-bound",
                                              "Lower-bound-bc"), 
                           conf = 0.95, se = TRUE) {
  method <- match.arg(method)
  
  if (method == "all") {
    a <- Chao2_sharedFun(y1, y2, conf)
    b <- Chao2_bcFun(y1, y2, conf)
    c <- PanFun.Sam(y1, y2, conf)
    d <- PanbcFun.Sam(y1, y2, conf)
    out <- rbind(a, b, c, d)
  }
  if (method == "Chao2-shared")
    out <- Chao2_sharedFun(y1, y2, conf)
  if (method == "Chao2-shared-bc")
    out <- Chao2_bcFun(y1, y2, conf)
  if (method == "Lower-bound") 
    out <- PanFun.Sam(y1, y2, conf)
  if (method == "Lower-bound-bc") 
    out <- PanbcFun.Sam(y1, y2, conf)
  
  if (se == FALSE) {
    out <- data.frame(Estimator = out[, 1], row.names = rownames(out))
  } 
  return(out)
}
