simulate = function(n, ngenes, maf, alpha, sigmab, sigmae, tau, gamma0){
	betas = rnorm(ngenes, 0, sigmab)

	genos = sapply(maf, function(x) rbinom(n, 1, x))
	effects = genos %*% betas + alpha

	l = rnorm(n, effects, sigmae)
	l2 = rnorm(n, effects, sigmae)
	l = l/sd(l)
	l2 = l2/sd(l2)
	y = l>tau
	y2 = l2>tau
	# yf = y
	# yf2 = y2
	# flip = which(!y)[1:(gamma0*sum(!y))]
	# yf[flip] = TRUE
	# flip = which(!y2)[1:(gamma0*sum(!y2))]
	# yf2[flip] = TRUE

	return(data.frame(y, y2, l, l2, genos))
}

## calculate h2 on the liability scale
h2l = function(df){
	K = mean(df$y)
	Kgr = mean(df$y2[df$y])
	x_pop = qnorm(1-K)
	x_gr = qnorm(1-Kgr)
	i_pop = dnorm(qnorm(1-K)) / K
	#(x_pop - x_gr) / i_pop
	(x_pop - x_gr*sqrt(1-(x_pop^2 - x_gr^2)*(1- x_pop/i_pop))) / (i_pop + x_gr^2*(i_pop - x_pop))
}

s## calculate h2*var(y) on the observed 1/0 scale
sigmau2 = function(df){
	K = mean(df$y)
	Kgr = mean(df$y2[df$y])
	x_pop = qnorm(1-K)
	x_gr = qnorm(1-Kgr)
	z = dnorm(qnorm(1-K))
	(x_pop - x_gr) * z * K
}

#### CLEAN SIMULATION
## to show that h2l can be used to predict sigmab
ngenes = 50
sigmab = 0.1
sigmae = 1.0
maf = rep(0.2, ngenes)

mmaf = sum(maf*(1-maf))
predratio = NULL
for(i in 1:300){
	df = simulate(5000, ngenes, maf, 0, sigmab, sigmae, qnorm(1-0.10), 0.05)

	h = h2l(df)
	sigmab_pred = sqrt((sigmae^2 * h / (1 - h) / mmaf))  # abs - yes or no?
	predratio = c(predratio, sigmab_pred/sigmab)
	print(sprintf("iteration %i: sigmab predicted is %f", i, sigmab_pred))
}
mean(predratio, na.rm=T)
qplot(predratio)


### NOW ON TO 1/0 SCALE
sigmau2(df)

ngenes = 50
sigmab = 0.3
sigmae = 1.0
maf = rep(0.2, ngenes)
mmaf = sum(maf*(1-maf))
sqrt(mmaf * sigmab^2 / (sigmab^2 * mmaf + sigmae^2))

df = simulate(5000, ngenes, maf, 0, sigmab, sigmae, qnorm(1-0.05), 0.05)
betalin = betalog = alog = NULL
for(c in 5:ncol(df)){
	betalin = c(betalin, coef(lm(df$l ~ df[,c]))[2])
	betalog = c(betalog, coef(glm(df$y ~ df[,c], family="binomial"))[2])
	alog = c(alog, coef(glm(df$y ~ df[,c], family="binomial"))[1])
}
mean(betalog)
summary(lm(betalog ~ betalin))
sd(df$y)

plot(betalog, betalin)
z = dnorm(qnorm(1-mean(df$y)))
sd(betalog)
var(betalog)
var(betalin)

h.2l(df) * z^.2 / mmaf
## below: IDIOT
## confirmed the equation u = c + zg, and u seems ~= 0.
## so log regr betas are ~ N(0, z^2*sigma_b^2).