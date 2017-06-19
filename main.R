## initialize

options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)

# TODO: replace this with proper cleanup script
load("/home/dominika/mfr.Rdata")

## load full mfr (ali's) for sib-pair identification
mfull = read.table("~/data/swed/ali/mfr14_aliids_onlyga.csv", h=T, sep=";")
mfull = filter(mfull, lpnr_barn>0, lpnr_mor>0, !is.na(lpnr_far))

sisters = filter(mfull, KON==2) %>%
	inner_join(., ., by=c("lpnr_mor", "lpnr_far")) %>%
	filter(lpnr_barn.x<lpnr_barn.y) %>%
	select(starts_with("lpnr_b"))
brothers = filter(mfull, KON==1) %>%
	inner_join(., ., by=c("lpnr_mor", "lpnr_far")) %>%
	filter(lpnr_barn.x<lpnr_barn.y) %>%
	select(starts_with("lpnr_b"))
colnames(sisters) = colnames(brothers) = c("sib1", "sib2")

# basic cleanup of MFR
mfr = filter(mfr, GRDBS>150, mult_pregn==0)

# leave only one sibling
onesib = group_by(mfr, lpnr_mor) %>% sample_n(size=1) %>% ungroup()
onesib = onesib[, c("lpnr_barn", "lpnr_mor", "lpnr_far", "GRDBS")]

# get mfr rows for maternal cousin pairs
matcousins = sisters %>%
	inner_join(onesib, by=c("sib1"="lpnr_mor")) %>%
	inner_join(onesib, by=c("sib2"="lpnr_mor"))

# get mfr rows for paternal cousin pairs
patcousins = brothers %>%
	inner_join(onesib, by=c("sib1"="lpnr_far")) %>%
	inner_join(onesib, by=c("sib2"="lpnr_far"))


# calculates heritability using Falconer's method
# args:
# 1. df - data frame with GRDBS.x and GRDBS.y columns, corresponding to GAs in some relative pairs
# 2. r - coefficient of relationship
# 3. thr - GA threshold to define cases, inclusive
heritability = function(df, r, thr)
{
	df = mutate(df, cases = GRDBS.x<=thr)
	gr1 = filter(df, cases) %>% mutate(cases2 = GRDBS.y<=thr)
	
	incidence_pop = mean(df$cases)
	incidence_gr = mean(gr1$cases2)
	
	# safety check for definitions that result in no cases
	if(incidence_pop == 0 | incidence_gr == 0) return(NA)
	
	x_pop = qnorm(1 - incidence_pop) 
	x_gr = qnorm(1 - incidence_gr)
	i_pop = dnorm(qnorm(1-incidence_pop)) / incidence_pop
	
	# h = (x_pop-x_gr)/i_pop/r - assuming equal variances in pop and gr
	h = (x_pop - x_gr*sqrt(1-(x_pop^2 - x_gr^2)*(1- x_pop/i_pop))) / r / (i_pop + x_gr^2*(i_pop - x_pop))
	return(h)
	
}

## STANDARD PTD h2 ESTIMATES
# h2 assuming only maternal genetic effects:
heritability(matcousins, 0.5, 37*7) 

# h2 assuming only paternal genetic effects:
heritability(patcousins, 0.5, 37*7) 

# h2 assuming only fetal genetic effects:
heritability(matcousins, 0.125, 37*7)
heritability(patcousins, 0.125, 37*7) 


## BOOTSTRAPPED ESTIMATES, STANDARD DEFINITION
bmat = boot(matcousins, function(x, i) heritability(x[i,], 0.5, 37*7), 50)
bpat = boot(patcousins, function(x, i) heritability(x[i,], 0.5, 37*7), 50)
bfet = boot(matcousins, function(x, i) heritability(x[i,], 0.125, 37*7), 50)

bind_rows("maternal"=data.frame(h2=bmat$t),
		  "paternal"=data.frame(h2=bpat$t),
		  "fetal"=data.frame(h2=bfet$t), .id="genome") %>%
	ggplot() + geom_boxplot(aes(x=genome, y=h2)) +
	theme_minimal() +
	geom_hline(yintercept = 0, col="red") +
	geom_hline(yintercept = 1, col="red") +
	scale_y_continuous(breaks=seq(-0.2, 1.6, by=0.1), minor_breaks=NULL)
summary(bmat$t)


## BOOTSTRAPPED ESTIMATES, SLIDING UPPER THRESHOLD
out = NULL
for(t in 160:300){
	print(sprintf("calculating h2 with threshold GA <= %i", t))
	kpop = mean(matcousins$GRDBS.x<=t)
	kgr = mean(matcousins$GRDBS.x[matcousins$GRDBS.y<=t]<=t)
	print(sprintf("K pop. = %.4f, K gr. = %.4f", kpop, kgr))
	booted = boot(matcousins, function(x, i) heritability(x[i,], 0.5, t), 100)
	out = bind_rows(out, data.frame(t=t, horig=booted$t0, hboot=booted$t))
}
write.table(out, "~/Documents/binaryh2/h2_estimates_upperthr.csv", col.names=T, row.names=F, quote=F)

out = NULL
for(c in 1:100){
	print(sprintf("iteration number %i", c))
	
	df = sample_frac(matcousins, replace=T)
	
	hs = data.frame(t=160:280, i=c)
	hs$h = sapply(hs$t, function(x) heritability(df, 0.5, x))
	out = bind_rows(out, hs)
}

group_by(out, i) %>%
	top_n(1, h) %>%
	ggplot(aes(x=t)) + geom_point(aes(y=h))

mean(onesib$GRDBS<=259)


## TODO:
# 1. find smallest heritability
# 2. use that to find sigma_b "true"
# 3. use K1 and K2 to find gamma_0
# 4. calculate phi factor for NCP'/NCP
# 5. assume some distribution for p
# 6. plug p, lambda_1 ~ N(0, z*sigma_b), lambda_2 = f(lambda_1, gamma_0) into NCP'/NCP
# 7. evaluate whether power was improved