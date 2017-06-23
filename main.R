## initialize

options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)

## raw mfr
mfull = read.table("~/data/swed/ver/mfr_150202_recoded.csv", sep=",", h=T, quote="")
par = read.table("~/data/swed/ver/parents.csv", h=T, sep=",")

# remove missing individual ids
mfull = mfull %>% filter(lpnr_barn!=0, lpnr_mor!=0, !is.na(lpnr_far)) 
mfull = inner_join(mfull, par, by=c("lpnr_BARN"="LopnrBarn", "lpnr_mor"="LopnrMor"))

# maternal nationality, leave only european mothers
mfr = mfull %>% filter(MFODLAND=="BELGIEN" | MFODLAND=="BOSNIEN-HERCEGOVINA" |  MFODLAND=="BULGARIEN" | MFODLAND=="DANMARK" | MFODLAND=="ESTLAND" | MFODLAND=="FINLAND" | MFODLAND=="FRANKRIKE"
                       | MFODLAND=="GREKLAND" | MFODLAND=="IRLAND" | MFODLAND=="ISLAND" | MFODLAND=="ITALIEN" | MFODLAND=="JUGOSLAVIEN" | MFODLAND=="LITAUEN" | MFODLAND=="LUXEMBURG"
                       | MFODLAND=="MOLDAVIEN" | MFODLAND=="NEDERLANDERNA" | MFODLAND=="NORGE"  | MFODLAND=="OSTERRIKE" | MFODLAND=="POLEN" | MFODLAND=="PORTUGAL" |
                         MFODLAND=="RUMANIEN" | MFODLAND=="SERBIEN" | MFODLAND=="SLOVAKIEN" | MFODLAND=="SLOVENIEN" | MFODLAND=="SPANIEN" | MFODLAND=="STORBRITANNIEN OCH NORDIRLAND" | 
                         MFODLAND=="SVERIGE" | MFODLAND=="SWAZILAND" | MFODLAND=="TURKIET" | MFODLAND=="TYSKA DEM REP (DDR)" | MFODLAND=="TYSKLAND" | MFODLAND=="UKRAINA")

# identify sibpairs
sisters = filter(mfull, KON==2) %>%
	inner_join(., ., by=c("lpnr_mor", "LopnrFar")) %>%
	filter(lpnr_BARN.x<lpnr_BARN.y) %>%
	select(starts_with("lpnr_B"))
brothers = filter(mfull, KON==1) %>%
	inner_join(., ., by=c("lpnr_mor", "LopnrFar")) %>%
	filter(lpnr_BARN.x<lpnr_BARN.y) %>%
	select(starts_with("lpnr_B"))
colnames(sisters) = colnames(brothers) = c("sib1", "sib2")

# basic cleanup of MFR
mfr = filter(mfr, GRDBS>150, BORDF2==1, GRMETOD!=0)

# adjusting for different GA timing methods
mod = lm(GRDBS ~ factor(GRMETOD), data=mfr)
mfr$GAadj = resid(mod) + coef(mod)[1]

# leave only one sibling
onesib = group_by(mfr, lpnr_mor) %>% sample_n(size=1) %>% ungroup()
onesib = onesib[, c("lpnr_BARN", "lpnr_mor", "LopnrFar", "GAadj")]

# get mfr rows for maternal cousin pairs
matcousins = sisters %>%
	inner_join(onesib, by=c("sib1"="lpnr_mor")) %>%
	inner_join(onesib, by=c("sib2"="lpnr_mor"))

# get mfr rows for paternal cousin pairs
patcousins = brothers %>%
	inner_join(onesib, by=c("sib1"="LopnrFar")) %>%
	inner_join(onesib, by=c("sib2"="LopnrFar"))


# calculates heritability using Falconer's method
# args:
# 1. df - data frame with a "cases" column
# 2. gr - data frame of the relatives of the "cases" from df, with a "cases2" column
# 3. r - coefficient of relationship
heritability = function(df, gr, r)
{
	incidence_pop = mean(df$cases)
	incidence_gr = mean(gr$cases2)
	
	# safety check for definitions that result in no cases
	if(incidence_pop == 0 | incidence_gr == 0) return(NA)
	
	x_pop = qnorm(1 - incidence_pop) 
	x_gr = qnorm(1 - incidence_gr)
	i_pop = dnorm(qnorm(1-incidence_pop)) / incidence_pop
	
	#h = (x_pop-x_gr)/i_pop/r #assuming equal variances in pop and gr
	h = (x_pop - x_gr*sqrt(1-(x_pop^2 - x_gr^2)*(1- x_pop/i_pop))) / r / (i_pop + x_gr^2*(i_pop - x_pop))
	return(h)
}

# args:
# 1. df - paired observation data frame with GAadj.x and GAadj.y columns
# 2. r - coefficient of relationship
# 3. thr - upper limit, inclusive, of GA to define cases
defineCasesUpp = function(df, r, thr){
	df = mutate(df, cases = GAadj.x<=thr)
	gr = filter(df, cases) %>% mutate(cases2 = GAadj.y<=thr)
	heritability(df, gr, r)
}
defineCasesLow = function(df, r, thr){
	df = mutate(df, cases = GAadj.x>thr & GAadj.x<=280)
	gr = filter(df, cases) %>% mutate(cases2 = GAadj.y>thr & GAadj.y<=280)
	heritability(df, gr, r)
}
defineCasesBoth = function(df, r, thr, thr2){
	df = mutate(df, cases = GAadj.x>thr & GAadj.x<=thr2)
	gr = filter(df, cases) %>% mutate(cases2 = GAadj.y>thr & GAadj.y<=thr2)
	heritability(df, gr, r)
}


### WHEN CASES ARE DEFINED AS GA<=THR

## STANDARD PTD h2 ESTIMATES
# h2 assuming only maternal genetic effects:
defineCasesUpp(matcousins, 0.5, 37*7) 

# h2 assuming only paternal genetic effects:
defineCasesUpp(patcousins, 0.5, 37*7) 

# h2 assuming only fetal genetic effects:
defineCasesUpp(matcousins, 0.125, 37*7)
defineCasesUpp(patcousins, 0.125, 37*7) 


## BOOTSTRAPPED ESTIMATES, STANDARD DEFINITION
bmat = boot(matcousins, function(x, i) defineCasesUpp(x[i,], 0.5, 37*7), 50)
bpat = boot(patcousins, function(x, i) defineCasesUpp(x[i,], 0.5, 37*7), 50)
bfet = boot(matcousins, function(x, i) defineCasesUpp(x[i,], 0.125, 37*7), 50)

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
	
	kpop = mean(matcousins$GAadj.x<=t)
	kgr = mean(matcousins$GAadj.y[matcousins$GAadj.x<=t]<=t)
	print(sprintf("K pop. = %.4f, K gr. = %.4f", kpop, kgr))
	if(kpop>0 & kgr>0){
		booted = boot(matcousins, function(x, i) defineCasesUpp(x[i,], 0.5, t), 100)
		out = bind_rows(out, data.frame(t=t, horig=booted$t0, hboot=booted$t))
	}
}
write.table(out, "~/Documents/binaryh2/h2_estimates_upperthr.csv", col.names=T, row.names=F, quote=F)


## BOOTSTRAPPED ESTIMATES, SLIDING LOWER THRESHOLD
out = NULL
for(t in 160:300){
	print(sprintf("calculating h2 with threshold %i < GA <= 280", t))
	
	caseindices = matcousins$GAadj.x>t & matcousins$GAadj.x<=280
	kpop = mean(caseindices)
	kgr = mean(matcousins$GAadj.y[caseindices]>t & matcousins$GAadj.y[caseindices]<=280)
	print(sprintf("K pop. = %.4f, K gr. = %.4f", kpop, kgr))
	if(kpop>0 & kgr>0){
		booted = boot(matcousins, function(x, i) defineCasesLow(x[i,], 0.5, t), 100)
		out = bind_rows(out, data.frame(t=t, horig=booted$t0, hboot=booted$t))
	}
}
write.table(out, "~/Documents/binaryh2/h2_estimates_lowerthr.csv", col.names=T, row.names=F, quote=F)

## BOOTSTRAPPED ESTIMATES, SLIDING BOTH THRESHOLDS
out = NULL
for(t2 in seq(160, 280, by=3)){
	for(t in seq(160, 280, by=3)){
		print(sprintf("calculating h2 with threshold %i < GA <= %i", t, t2))
		
		caseindices = matcousins$GAadj.x>t & matcousins$GAadj.x<=t2
		kpop = mean(caseindices)
		kgr = mean(matcousins$GAadj.y[caseindices]>t & matcousins$GAadj.y[caseindices]<=t2)
		print(sprintf("K pop. = %.4f, K gr. = %.4f", kpop, kgr))
		if(kpop>0 & kgr>0){
			booted = boot(matcousins, function(x, i) defineCasesBoth(x[i,], 0.5, t, t2),
						  100, parallel="multicore", ncpus=7)
			out = bind_rows(out, data.frame(tl=t, tu=t2, horig=booted$t0, hboot=booted$t))
		}
	}	
}

write.table(out, "~/Documents/binaryh2/h2_estimates_boththr2.csv", col.names=T, row.names=F, quote=F)


## temp
tmp=NULL
for(t2 in seq(160, 280, by=3)){
	for(t in seq(160, 280, by=3)){
		print(sprintf("calculating h2 with threshold %i < GA <= %i", t, t2))
		k = mean(matcousins$GAadj.x>t & matcousins$GAadj.x<=t2)
		
		tmp = bind_rows(tmp, data.frame(tl=t, tu=t2, k=k))
	}	
}

obs = inner_join(out, tmp, by=c("tl", "tu")) %>%
	group_by(tl, tu) %>%
	summarize(h2l=min(horig), K=min(k)) %>%
	mutate(N = K*3000*10, v=K*10)
write.table(obs, "~/Documents/binaryh2/h2_estimate_fullout.csv", col.names=T, row.names=F, quote=F)

## TODO:
# 1. find smallest heritability
# 2. use that to find sigma_b "true"
# 3. use K1 and K2 to find gamma_0, or calculate lambda_2 like lambda_1
# 4. calculate phi factor for NCP'/NCP
# 5. assume some distribution for p
# 6. plug p, lambda_1 ~ N(0, z*sigma_b), lambda_2 = f(lambda_1, gamma_0) into NCP'/NCP
# 7. evaluate whether power was improved
