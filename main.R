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

# get mfr rows for maternal cousin pairss
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
	
	x_pop = qnorm(1 - incidence_pop) 
	x_gr = qnorm(1 - incidence_gr)
	i_pop = dnorm(qnorm(1-incidence_pop)) / incidence_pop
	
	h = (x_pop-x_gr)/i_pop/r
	return(h)
	
}

# h2 assuming only maternal genetic effects:
heritability(matcousins, 0.5, 37*7) 

# h2 assuming only paternal genetic effects:
heritability(patcousins, 0.5, 37*7) 

# h2 assuming only fetal genetic effects:
heritability(matcousins, 0.125, 37*7)
heritability(patcousins, 0.125, 37*7) 

# bootstrap h2 assuming maternal eff. only:
bmat = boot(matcousins, function(x, i) heritability(x[i,], 0.5, 37*7), 50)
# bootstrap h2 assuming paternal eff. only:
bpat = boot(patcousins, function(x, i) heritability(x[i,], 0.5, 37*7), 50)
# bootstrap h2 assuming fetal eff. only:
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

