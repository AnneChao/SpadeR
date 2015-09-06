print.ChaoShared <- function(x, ...){
	cat('\n(1) BASIC DATA INFORMATION:\n\n')
	
	if(nrow(x[[2]])==6){
	  cat("              (Number of observed individuals in community 1)   n1  = ", x[[1]]$n1, "\n")
	  cat("              (Number of observed individuals in community 2)   n2  = ", x[[1]]$n2, "\n")
	  cat("              (Number of observed species in community 1)       D1  = ", x[[1]]$D1, "\n")
	  cat("              (Number of observed species in community 2)       D2  = ", x[[1]]$D2, "\n")
	  cat("              (Number of observed shared species)               D12 = ", x[[1]]$D12, "\n")
	  cat("              (Bootstrap replications for s.e. estimate)       ", x[[1]]$B,   "\n\n")
	  cat("     \"Rare\" Shared Species Group: (Both Frequencies can only up to 10)", "\n")
	  cat("         Some Statistics:", "\n")
	  cat("         --------------------------------------------------------------", "\n")
	  cat("         f[11] =", x[[1]]$f11, "; ", "f[1+] =", x[[1]]$f1.plus, ";", "f[+1] =", x[[1]]$fplus.1, "; ", "f[2+] =", x[[1]]$f2.plus, "; ", "f[+2] =", x[[1]]$fplus.2, "\n")
	  cat("         --------------------------------------------------------------", "\n")
	  cat("            (Number of observed individuals in community 1)     n1_rare   = ", x[[1]]$n1_rare, "\n")
	  cat("            (Number of observed individuals in community 2)     n2_rare   = ", x[[1]]$n2_rare, "\n")
	  cat("            (Number of observed shared species)                 D12_rare  = ", x[[1]]$D12_rare, "\n")
	  cat("            (Estimated sample coverage)                         C12_rare  = ", x[[1]]$C12_rare, "\n")
	  cat("            (Estimated CCVs)                                    CCV_1     = ", x[[1]]$CCV_1, "\n")
	  cat("                                                                CCV_2     = ", x[[1]]$CCV_2, "\n")
	  cat("                                                                CCV_12    = ", x[[1]]$CCV_12, "\n")
  
    }else{
      #cat("(1)  BASIC DATA INFORMATION:", "\n\n")
      cat("                      (Number of samples from community 1)   t1  = ", x[[1]]$t1, "\n")
      cat("                      (Number of samples from community 2)   t2  = ", x[[1]]$t2, "\n")
      cat("               (Number of observed species in community 1)   D1  = ", x[[1]]$D1, "\n")
      cat("               (Number of observed species in community 2)   D2  = ", x[[1]]$D2, "\n")
      cat("     (Number of observed shared species in two communities)  D12 = ", x[[1]]$D12, "\n")
      cat("              (Bootstrap replications for s.e. estimate)    ", x[[1]]$B,   "\n\n")
      cat("     Some Statistics:", "\n")
      cat("         --------------------------------------------------------------", "\n")
      cat("         Q[11] =", x[[1]]$Q11, "; ", "Q[1+] =", x[[1]]$Q1.plus, ";", "Q[+1] =", x[[1]]$Qplus.1, "; ", "Q[2+] =", x[[1]]$Q2.plus, "; ", "Q[+2] =", x[[1]]$Qplus.2, "\n")
      cat("         --------------------------------------------------------------", "\n")
      
	}
  
	cat('\n')
	cat('\n(2) ESTIMATION RESULTS OF THE NUMBER OF SHARED SPECIES:\n\n')
	print(round(x$ESTIMATION_RESULTS_OF_THE_NUMBER_OF_SHARED_SPECIES,3))
	cat('\n')
	if(nrow(x[[2]])==6){
	cat('
(3) DESCRIPTION OF MODELS FOR ESTIMATING SHARED SPECIES RICHNESS:
Homogeneous: This model assumes that the shared species in each community have the same discovery probabilities; see the Eq.(3.11a) of Chao et al. (2000).

Heterogeneous (ACE-shared): This model allows for heterogeneous discovery probabilities among shared species; see Eq.(3.11b) of Chao et al. (2000). It is an extension of ACE to two communities. It is replaced by Chao1-shared when the estimated sample coverage (C12_rare in the output) is zero.

Chao1-shared: An extension of the Chao1 estimator to two communities. It provides an estimate of shared species richness. See Eq. (11) of Chao, Shen and Hwang (2006). It is replaced by Chao1-shared-bc (see below) for the case f[2+]=0 or f[+2]=0.

Chao1-shared-bc: A bias-corrected form of the Chao1-shared. See Eq. (12) of Chao, Shen and Hwang (2006).

Lower-bound: See Pan et al. (2009)
   
Lower-bound-bc: See Pan et al. (2009)
	')
	}else{
	cat('
(3) DESCRIPTION OF MODELS FOR ESTIMATING SHARED SPECIES RICHNESS:

Chao2-shared: An extension of the Chao2 estimator to two communities. It provides an estimate of shared species richness. See Eq. (13) of Chao, Shen and Hwang (2006). It is replaced by Chao2-shared-bc (see below) for the case Q[2+]=0 or Q[+2]=0.

Chao2-shared-bc:  A bias-corrected form of Chao2-shared. See Eq. (14) of Chao, Shen and Hwang (2006).

Lower-bound: See Pan et al. (2009)

Lower-bound-bc: See Pan et al. (2009)
	')}
}

