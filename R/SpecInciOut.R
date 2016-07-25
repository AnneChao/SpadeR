SpecInciOut <-
function(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"), k, conf){
    a <- as.numeric(round(SpecInciHomo(data, k, conf), 3))
    b <- as.numeric(round(SpecInciChao2(data, k, conf), 3))
    c <- as.numeric(round(SpecInciChao2bc(data, k, conf), 3))
    d <- as.numeric(round(SpecInciiChao2(data, k, conf), 3))
    e <- as.numeric(round(SpecInciModelh(data, k, conf)[1:4], 3))
    f <- as.numeric(round(SpecInciModelh1(data, k, conf)[1:4], 3))
    g <- as.numeric(round(SpecInciJack1(data, k, conf), 3))
    h <- as.numeric(round(SpecInciJack2(data, k, conf), 3))
    est.cv <- data.frame(c("", "", "", "", round(SpecInciModelh(data, k, conf)[5], 3),
                           round(SpecInciModelh1(data, k, conf)[5], 3), "", ""))
    colnames(est.cv) <- "Est.CV(rare)"
  if (method == "all") {
    #out <- data.frame(rbind(a, b, c, d, e, f, g, h), est.cv)
	out <- data.frame(rbind(a, b, c, d, e, f, g, h))
	rownames(out) <- c("    Homogeneous Model",
						"    Chao2 (Chao, 1987)",
						"    Chao2-bc",
						"    iChao2 (Chiu et al. 2014)",
						"    ICE (Lee & Chao, 1994)",
						"    ICE-1 (Lee & Chao, 1994)",
						"    1st order jackknife",
						"    2nd order jackknife")
	colnames(out) <- c("Estimate", "s.e.", "95%Lower", "95%Upper")
  }

  if (method == "Homogeneous")
    out <- a
  if (method == "Chao")
    out <- rbind(b, c, d)
  if (method == "CE"){
    est.cv <- matrix(c(SpecInciModelh(data, k, conf)[5], SpecInciModelh1(data, k, conf)[5]), ncol = 1)
    est.cv <- round(est.cv, 3)
    colnames(est.cv) <- "Est.CV(rare)"
    out1 <- rbind(e, f)
    out <- cbind(out1, est.cv)
  }
  if (method == "Jackknife")
    out <- rbind(g, h)
  colnames(out) <- c("Estimate", "s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
  return(out)
}

#####J.H Lin 2016.7.11 #####add raw data input
SpecInciOut_raw <-
  function(data, method = c("all", "Homogeneous", "Chao", "CE", "Jackknife"), k, conf){
    a <- round(SpecInciHomo(data[-2], k, conf), 3)
    b <- round(SpecInciChao2(data[-1], k, conf), 3)
    c <- round(SpecInciChao2bc(data[-1], k, conf), 3)
    d <- round(SpecInciiChao2(data[-1], k, conf), 3)
    e <- round(SpecInciModelh(data[-2], k, conf)[1:4], 3)
    f <- round(SpecInciModelh1(data[-2], k, conf)[1:4], 3)
    g <- round(SpecInciJack1(data[-1], k, conf), 3)
    h <- round(SpecInciJack2(data[-1], k, conf), 3)
    est.cv <- data.frame(c("", "", "", "", round(SpecInciModelh(data[-2], k, conf)[5], 3),
                           round(SpecInciModelh1(data[-2], k, conf)[5], 3), "", ""))
    colnames(est.cv) <- "Est.CV(rare)"
    if (method == "all") {
      #out <- data.frame(rbind(a, b, c, d, e, f, g, h), est.cv)
      out <- data.frame(rbind(a, b, c, d, e, f, g, h), row.names = NULL)
      rownames(out) <- c("    Homogeneous Model",
                         "    Chao2 (Chao, 1987)",
                         "    Chao2-bc",
                         "    iChao2 (Chiu et al. 2014)",
                         "    ICE (Lee & Chao, 1994)",
                         "    ICE-1 (Lee & Chao, 1994)",
                         "    1st order jackknife",
                         "    2nd order jackknife")
      colnames(out) <- c("Estimator", "Est_s.e.", "95% Lower Bound", "95% Upper Bound")
    }

    if (method == "Homogeneous")
      out <- a
    if (method == "Chao")
      out <- rbind(b, c, d)
    if (method == "CE"){
      est.cv <- matrix(c(SpecInciModelh(data[-2], k, conf)[5], SpecInciModelh1(data[-2], k, conf)[5]), ncol = 1)
      est.cv <- round(est.cv, 3)
      colnames(est.cv) <- "Est.CV(rare)"
      out1 <- rbind(e, f)
      out <- cbind(out1, est.cv)
    }
    if (method == "Jackknife")
      out <- rbind(g, h)
    colnames(out) <- c("Estimate", "   s.e.", paste(conf*100,"%Lower",sep=""), paste(conf*100,"%Upper",sep=""))
    return(out)
  }
