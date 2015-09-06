ChaoShared.Ind <-
function(x1, x2, method = c("all", "Homogeneous", 
                                              "Heterogeneous(ACE-shared)",
                                              "Chao1(shared)", 
                                              "Chao1-shared-bc", 
                                              "Lower-bound",
                                              "Lower-bound-bc"), 
                           B = 200, conf = 0.95, se = TRUE) {
  method <- match.arg(method)
  if (se == FALSE)
    B <- 1
  if (method == "all") {
    a <- HomoFun(x1, x2, B, conf)
    b <- HeteroFun(x1, x2, B, conf)
    c <- Chao1_sharedFun(x1, x2, conf)
    d <- Chao1_bcFun(x1, x2, conf)
    e <- PanFun(x1, x2, conf)
    f <- PanbcFun(x1, x2, conf)
    out <- rbind(a, b, c, d, e, f)
  }
  if (method == "Homogeneous")
    out <- HomoFun(x1, x2, B, conf)
  if (method == "Heterogeneous(ACE-shared)")
    out <- HeteroFun(x1, x2, B, conf)
  if (method == "Chao1(shared)")
    out <- Chao1_sharedFun(x1, x2, conf)
  if (method == "Chao1-shared-bc")
    out <- Chao1_bcFun(x1, x2, conf)
  if (method == "Lower-bound)")
    out <- PanFun(x1, x2, conf)
  if (method == "Lower-bound-bc")
    out <- PanbcFun(x1, x2, conf)
  
  if (se == FALSE) {
    out <- data.frame(Estimator = out[, 1], row.names = rownames(out))
  } 
  return(out)
}
