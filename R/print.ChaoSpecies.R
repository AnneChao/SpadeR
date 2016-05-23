print.ChaoSpecies <- function(x, ...){
	cat('\n(1) BASIC DATA INFORMATION:\n\n')
	print(x$Basic.Data.Information)
  
	cat('\n')
	if(names(x)[2]=="Rare.Species.Group"){
		print(x$Rare.Species.Group)
		cat('\n')
		cat('\n(2) SPECIES RICHNESS ESTIMATORS TABLE:\n\n')
		print(x$Species.Table)
		cat('\n')
		cat('\n(3) DESCRIPTION OF ESTIMATORS/MODELS:

Homogeneous Model: This model assumes that all species have the same abundances or discovery probabilities. See Eq. (2.3) of Chao and Lee (1992) or Eq. (2.1) of Chao et al. (2000).

Homogeneous (MLE): An approximate maximum likelihood estimate under homogeneous model. See Eq. (1.1) and Eq. (1.2) of Chao and Lee (1992).

Chao1 (Chao, 1984): This approach uses the numbers of singletons and doubletons to estimate the number of undetected species because undetected species information is mostly concentrated on those low frequency counts; see Chao (1984), Shen, Chao and Lin (2003) and Chao and Chiu (2012).

Chao1-bc: A bias-corrected form for the Chao1 estimator; see Chao (2005).

iChao1: An improved Chao1 estimtor; see Chiu et al. (2014).

ACE (Abundance-based Coverage Estimator): A non-parametric estimator proposed by Chao and Lee (1992) and Chao, Ma and Yang (1993). The observed species are separated as rare and abundant groups; only the rare group is used to estimate the number of undetected species. The estimated CV is used to characterize the degree of heterogeneity among species discovery probabilities. See Eq. (2.14) in Chao and Lee (1992) or Eq. (2.2) of Chao et al. (2000).

ACE-1: A modified ACE for highly-heterogeneous communities. See Eq. (2.15) in Chao and Lee (1992).

1st order jackknife: It uses the number of singletons to estimate the number of undetected species; see Burnham and Overton (1978).

2nd order jackknife: It uses the numbers of singletons and doubletons to estimate the number of undetected species; see Burnham and Overton (1978).

95% Confidence interval: A log-transformation is used for all estimators so that the lower bound of the resulting interval is at least the number of observed species. See Chao (1987).
')}
	if(names(x)[2]!="Rare.Species.Group"){
		print(x$Infreq.Species.Group)
		cat('\n')
		cat('\n(2) SPECIES RICHNESS ESTIMATORS TABLE:\n\n')
		print(x$Species.Table)
		cat('\n')
		cat('\n(3) DESCRIPTION OF ESTIMATORS/MODELS:

Homogeneous Model: This model assumes that all species have the same incidence or detection probabilities. See Eq. (3.2) of Lee and Chao (1994).

Chao2 (Chao, 1987): This approach uses the frequencies of uniques and duplicates to estimate the number of undetected species; see Chao (1987).
     
Chao2-bc: A bias-corrected form for the Chao2 estimator; see Chao (2005).
  
iChao2: An improved Chao2 estimator; see Chiu et al. (2014).

ICE (Incidence-based Coverage Estimator): This estimator assumes that the incidence or detection probabilities are heterogeneous among species, but each probability is kept as a constant over sampling units. The estimator given in Gotelli and Chao (2013), which is an improved version of Eq. (3.18) in Lee and Chao (1994).

ICE-1: A modified ICE for highly-heterogeneous cases.

1st order jackknife: It uses the frequency of uniques to estimate the number of undetected species; see Burnham and Overton (1978).

2nd order jackknife: It uses the frequencies of uniques and duplicates to estimate the number of undetected species; see Burnham and Overton (1978).

95% Confidence interval: A log-transformation is used for all estimators so that the lower bound of the resulting interval is at least the number of observed species. See Chao (1987).
')
	}
}
	

	