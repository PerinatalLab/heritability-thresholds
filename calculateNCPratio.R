## this part converts observed lambdas etc. to NCP'/NCP ratios
## two options for calculating lambda2:
##  can either use e^beta/(e^beta + gamma_0),
##  or calculate from heritability under the alternative case definition.
## seems with b high prevalence trumps everything,
## so using a

obs = read.table("~/Documents/binaryh2/h2_estimate_fullout.csv", h=T)


# Epi(1-pi) * m
maf = runif(100)
mmaf = sum(maf * (1-maf)) # 1/6 * nsnps given uniform distr
mmaf = 100 / 6

# find sigma^2_gg (variance of genetic betas on the lin scale)
obs$sigma2_gg = obs$h2l / mmaf

# find sigma^2_uu (variance of genetic betas on the log scale), used only for ref
obs$sigma2_uu = obs$sigma2_gg * dnorm( qnorm( 1-obs$K ) )^2

# estimate of exp(alpha)
obs$alphaest = obs$K / (1-obs$K)


## for each pair of definitions that are being compared:
# dp = 1
integral1 = function(p, l1, l2){
	((1 + p * (l1-1)) / (1 + p * (l2-1)))^2
}
# dBeta1 = dnorm, kappa = sigma_bb2/sigma_bb1
integral2 = Vectorize(function(b1, sigma_uu2, sigma_uu1){
	if(b1==0) b1 = 1e-6
	l1 = exp(b1)
	l2 = exp(b1 * sigma_uu2 / sigma_uu1)
	((l2-1) / (l1-1)) ^ 2 * integrate(integral1, 0, 1, l1=l1, l2=l2, rel.tol=1e-7)$value * dnorm(b1, 0, sigma_uu1)
}, vectorize.args="b1")


## calculates actual NCP'/NCP for each thing against ref
## (vectorize later)

moba = read.table("~/data/geno/newQC/mothers_clean_gestage.pheno")

obs$ncpratio = NA
for(i in 1:nrow(obs)){
	# choose reference
	if(obs$tu[i] <= 259){
		r = which(obs$tl == 160 & obs$tu == 259)
	} else {
		r = which(obs$tl == 160 & obs$tu == 280)
	}
	
	# get reference params
	K1 = obs$K[r]
	N1 = sum(moba$V3>obs$tl[r] & moba$V3<=obs$tu[r])
	v1 = N1/nrow(moba)
	
	# get alternative params
	N2 = sum(moba$V3>obs$tl[i] & moba$V3<=obs$tu[i])
	v2 = N2/nrow(moba)
	
	if(N2 == 0){
		obs$ncpratio[i] = 1
		next
	}
	
	# N2 must be subset of N1
	gamma_0 = (N1 - N2) / (nrow(moba) - N2)
	# bias factor = kappa
	biasfactor = obs$alphaest[i]/(obs$alphaest[i] + gamma_0)
	# sigma2 assumed to be correct, sigma1 assumed to be biased
	sigma_uu2 = sqrt(max(obs$sigma2_uu[i], 1e-7))
	sigma_uu1 = biasfactor * sigma_uu2
	
	## main calculation
	phin = v2 * (1-v2) * N2 * (1-K1)^2
	phid = v1 * (1-v1) * N1 * (1-obs$K[i])^2
	
	fullIntegral = integrate(integral2, -10*sigma_uu2, 10*sigma_uu2, sigma_uu2=sigma_uu2, sigma_uu1,
							 rel.tol=1e-10, subdivisions=1000)$value
	
	obs$phi[i] = phin / phid
	obs$fullint[i] = fullIntegral
	obs$ncpratio[i] = fullIntegral * phin / phid
}

write.table(obs, "~/Documents/binaryh2/ncp_estimates_100snpsUnifp.csv", quote=F, col.names=T, row.names=F)

mobacc = mutate(moba, V3 = (V3<=235 & V3>220) + 1)
write.table(mobacc, "~/Documents/binaryh2/cases220-235.pheno", col.names=F, row.names=F, quote=F)
