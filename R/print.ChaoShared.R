print.ChaoShared <- function(x, ...){
	cat('\n(1) BASIC DATA INFORMATION:\n\n')
	
	if(nrow(x[[2]])==4){
	  cat("    Sample size in Community 1                      n1  = ", x[[1]]$n1, "\n")
	  cat("    Sample size in Community 2                      n2  = ", x[[1]]$n2, "\n")
	  cat("    Number of observed species in Community 1       D1  = ", x[[1]]$D1, "\n")
	  cat("    Number of observed species in Community 2       D2  = ", x[[1]]$D2, "\n")
	  cat("    Number of observed shared species               D12 = ", x[[1]]$D12, "\n")
	  cat("    Bootstrap replications for s.e. estimate              ", x[[1]]$B,   "\n\n")
	  cat("    \"Entire\" Shared Species Group:", "\n")
	  cat("         Some statistics:", "\n")
	  cat("         ---------------------------------------------------------------------------", "\n")
	  cat("         f[11] =", x[[1]]$f11, "; ", "f[1+] =", x[[1]]$f1.plus, ";", "f[+1] =", x[[1]]$fplus.1, "; ", "f[2+] =", x[[1]]$f2.plus, "; ", "f[+2] =", x[[1]]$fplus.2, ";" , "f[22] =", x[[1]]$f22 , "\n")
	  cat("         ---------------------------------------------------------------------------", "\n\n")
    cat("    \"Rare\" Shared Species Group: (Both frequencies can only up to 10)", "\n")
	  cat("         Some statistics:", "\n")
	  cat("         -------------------------------------------------------------------", "\n")
	  cat("         f[1+]_rare =", x[[1]]$f1.plus.rare, ";", "f[+1]_rare =", x[[1]]$fplus.1.rare, "; ", "f[2+]_rare =", x[[1]]$f2.plus.rare, "; ", "f[+2]_rare =", x[[1]]$fplus.2.rare, "\n")
	  cat("         -------------------------------------------------------------------", "\n")
	  cat("    Number of observed individuals in Community 1     n1_rare   = ", x[[1]]$n1_rare, "\n")
	  cat("    Number of observed individuals in Community 2     n2_rare   = ", x[[1]]$n2_rare, "\n")
	  cat("    Number of observed shared species                 D12_rare  = ", x[[1]]$D12_rare, "\n")
	  cat("    Estimated sample coverage                         C12_rare  = ", round(x[[1]]$C12_rare,3), "\n")
	  cat("    Estimated CCVs                                    CCV_1     = ", round(x[[1]]$CCV_1,3), "\n")
	  cat("                                                      CCV_2     = ", round(x[[1]]$CCV_2,3), "\n")
	  cat("                                                      CCV_12    = ", round(x[[1]]$CCV_12,3), "\n")
  
    }else{
      #cat("(1)  BASIC DATA INFORMATION:", "\n\n")
      cat("    Number of sampling units in Community 1               T1  = ", x[[1]]$T1, "\n")
      cat("    Number of sampling units in Community 2               T2  = ", x[[1]]$T2, "\n")
      cat("    Number of total incidences in Community 1             U1  = ", x[[1]]$U1, "\n")
      cat("    Number of total incidences in Community 2             U2  = ", x[[1]]$U2, "\n")
      cat("    Number of observed species in Community 1             D1  = ", x[[1]]$D1, "\n")
      cat("    Number of observed species in Community 2             D2  = ", x[[1]]$D2, "\n")
      cat("    Number of observed shared species in two communities  D12 = ", x[[1]]$D12, "\n")
      cat("    Bootstrap replications for s.e. estimate                    ", x[[1]]$B,   "\n\n")
      cat("     Some statistics:", "\n")
      cat("         --------------------------------------------------------------------------", "\n")
      cat("         Q[11] =", x[[1]]$Q11, "; ", "Q[1+] =", x[[1]]$Q1.plus, ";", "Q[+1] =", x[[1]]$Qplus.1, "; ", "Q[2+] =", x[[1]]$Q2.plus, "; ", "Q[+2] =", x[[1]]$Qplus.2, ";" , "Q[22] =", x[[1]]$Q22, "\n")
      cat("         --------------------------------------------------------------------------", "\n")
      
	}
  
	cat('\n')
	cat('\n(2) ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES:\n\n')
	print(round(x$Estimation_results,3))
	cat('\n')
	if(nrow(x[[2]])==4){
	cat('
(3) DESCRIPTION OF MODELS FOR ESTIMATING SHARED SPECIES RICHNESS:

Homogeneous: This model assumes that the shared species in each community have the same discovery probabilities; see the Eq. (3.11a) of Chao et al. (2000).

Heterogeneous (ACE-shared): This model allows for heterogeneous discovery probabilities among shared species; see Eq. (3.11b) of Chao et al. (2000). It is an extension of the ACE estimator to two communities. It is replaced by Chao1-shared when the estimated sample coverage for rare shared species group (C12_rare in the output) is zero.

Chao1-shared: An extension of the Chao1 estimator to estimate shared species richness between two communities. It provides a lower bound of shared species richness. See Eq. (3.6) of Pan et al. (2009). It is replaced by Chao1-shared-bc for the case f[2+]=0 or f[+2]=0.
   
Chao1-shared-bc: A bias-corrected form of Chao1-shared estimator; See Pan et al. (2009).
	')
	}else{
	cat('
(3) DESCRIPTION OF MODELS FOR ESTIMATING SHARED SPECIES RICHNESS:

Chao2-shared: An extension of the Chao2 estimator to estimate shared species richness between two communities. It provides a lower bound of shared species richness. See Pan et al. (2009). It is replaced by Chao2-shared-bc for the case Q[2+]=0 or Q[+2]=0.

Chao2-shared-bc: A bias-corrected form of Chao2-shared. See Pan et al. (2009).
	')}
}

