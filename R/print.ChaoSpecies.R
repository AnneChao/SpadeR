print.ChaoSpecies <- function(x, ...){
  cat('\n(1) BASIC DATA INFORMATION:\n\n')
  if(names(x)[2]=="Rare.species.group"){
    print(x$Basic_data_information[1:4,])
    cat('\n')
    print(x$Basic_data_information[-c(1:4),])
    cat('\n')
    print(x$Rare_species.group)
    cat('\n')
    cat('\n(2) SPECIES RICHNESS ESTIMATORS TABLE:\n\n')
    print(x$Species_table)
    cat('\n')
    cat('\n(3) DESCRIPTION OF ESTIMATORS/MODELS:

Homogeneous Model: This model assumes that all species have the same abundances or discovery probabilities. See Eq. (2.3) of Chao and Lee (1992) or Eq. (7a) of Chao and Chiu (2016b).

Homogeneous (MLE): An approximate maximum likelihood estimate under homogeneous model. See Eq. (1.1) and Eq. (1.2) of Chao and Lee (1992) or Eq. (3) of Chao and Chiu (2016b).

Chao1 (Chao, 1984): This approach uses the numbers of singletons and doubletons to estimate the number of undetected species because undetected species information is mostly concentrated on those low frequency counts; see Chao (1984), and Chao and Chiu (2012, 2016a,b).

Chao1-bc: A bias-corrected form for the Chao1 estimator; see Chao (2005) or Eq. (6b) of Chao and Chiu (2016b).

iChao1: An improved Chao1 estimtor; see Chiu et al. (2014).

ACE (Abundance-based Coverage Estimator): A non-parametric estimator proposed by Chao and Lee (1992) and Chao, Ma and Yang (1993). The observed species are separated as rare and abundant groups; only data in the rare group is used to estimate the number of undetected species. The estimated CV of the species in rare group characterizes the degree of heterogeneity among species discovery probabilities. See Eq. (2.14) in Chao and Lee (1992) or Eq. (7b) of Chao and Chiu (2016b).

ACE-1: A modified ACE for highly-heterogeneous communities when CV of the entire dataset > 2 and species richness > 1000. See Eq. (2.15) in Chao and Lee (1992).

1st order jackknife: It uses the number of singletons to estimate the number of undetected species; see Burnham and Overton (1978).

2nd order jackknife: It uses the numbers of singletons and doubletons to estimate the number of undetected species; see Burnham and Overton (1978).

95% Confidence interval: A log-transformation is used for all estimators so that the lower bound of the resulting interval is at least the number of observed species. See Chao (1987).
')}
	if(names(x)[2]!="Rare.species.group"){
	  print(x$Basic_data_information[1:5,])
	  cat('\n')
	  print(x$Basic_data_information[-c(1:5),])
	  cat('\n')
	  print(x$Infreq_species_group)
	  cat('\n')
	  cat('\n(2) SPECIES RICHNESS ESTIMATORS TABLE:\n\n')
	  print(x$Species_table)
	  cat('\n')
	  cat('\n(3) DESCRIPTION OF ESTIMATORS/MODELS:

Homogeneous Model: This model assumes that all species have the same incidence or detection probabilities. See Eq. (3.2) of Lee and Chao (1994) or Eq. (12a) in Chao and Chiu (2016b).

Chao2 (Chao, 1987): This approach uses the frequencies of uniques and duplicates to estimate the number of undetected species; see Chao (1987) or Eq. (11a) in Chao and Chiu (2016b).
     
Chao2-bc: A bias-corrected form for the Chao2 estimator; see Chao (2005).
  
iChao2: An improved Chao2 estimator; see Chiu et al. (2014).

ICE (Incidence-based Coverage Estimator): A non-parametric estimator originally proposed by Lee and Chao (1994) in the context of capture-recapture data analysis. The observed species are separated as frequent and infrequent species groups; only data in the infrequent group are used to estimate the number of undetected species. The estimated CV for species in the infrequent group characterizes the degree of heterogeneity among species incidence probabilities. See Eq. (12b) of Chao and Chiu (2016b), which is an improved version of Eq. (3.18) in Lee and Chao (1994). This model is also called Model(h) in capture-recapture literature where h denotes "heterogeneity".

ICE-1: A modified ICE for highly-heterogeneous cases.

1st order jackknife: It uses the frequency of uniques to estimate the number of undetected species; see Burnham and Overton (1978).

2nd order jackknife: It uses the frequencies of uniques and duplicates to estimate the number of undetected species; see Burnham and Overton (1978).

95% Confidence interval: A log-transformation is used for all estimators so that the lower bound of the resulting interval is at least the number of observed species. See Chao (1987).
')
	}
}
	

	
