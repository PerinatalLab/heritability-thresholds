options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)

out = read.table("~/Documents/binaryh2/h2_estimates_upperthr.csv", h=T)
out = filter(out, t<281, !is.na(hboot))

ggplot(out, aes(x=t)) + geom_boxplot(aes(y=hboot, group=t))
ggplot(out, aes(x=t)) + geom_point(aes(y=hboot)) + geom_smooth(method="loess", aes(y=hboot))
group_by(out, t) %>% summarize(m=mean(hboot), s=sd(hboot)) %>%
	ggplot(aes(x=t)) + geom_pointrange(aes(y=m, ymin=m-s, ymax=m+s))


out = read.table("~/Documents/binaryh2/h2_estimates_lowerthr.csv", h=T)
out = filter(out, t<281, !is.na(hboot))

ggplot(out, aes(x=t)) + geom_boxplot(aes(y=hboot, group=t))
ggplot(out, aes(x=t)) + geom_point(aes(y=hboot)) + geom_smooth(method="loess", aes(y=hboot))
group_by(out, t) %>% summarize(m=mean(hboot), s=sd(hboot)) %>%
	ggplot(aes(x=t)) + geom_pointrange(aes(y=m, ymin=m-s, ymax=m+s))


out = read.table("~/Documents/binaryh2/h2_estimates_boththr2.csv", h=T)
group_by(out, tu, tl) %>% summarize(m=max(mean(horig, na.rm=T), 0), s=sd(hboot, na.rm=T)) %>%
	ggplot(aes(x=tl, y=tu)) + geom_point(aes(col=m), size=5, shape = 15) +
	scale_color_continuous(low="red", high="green")

filter(out, tu==232) %>% group_by(tl) %>% summarize(m=mean(horig, na.rm=T), s=sd(hboot, na.rm=T)) %>%
	ggplot(aes(x=tl)) + geom_pointrange(aes(y=m, ymin=m-s, ymax=m+s))

obs = read.table("/mnt/HUNT/ncp_estimates_100snpsUnifp.csv", h=T)

ggplot(obs) + geom_point(aes(x=tl, y=tu, col=ncpratio), shape=15, size=5) +
	scale_color_continuous(low="red", high="green")

ggplot(obs) + geom_point(aes(x=tl, y=tu, col=h2l), shape=15, size=5) +
	scale_color_continuous(low="red", high="green")

filter(obs, tu==235) %>%
	ggplot(aes(x=tl)) + geom_point(aes(y=h2l))


gwas0 = read.table("~/Documents/binaryh2/casesALL.assoc", h=T)
gwas1 = read.table("~/Documents/binaryh2/cases220-235.assoc", h=T)

gwas = inner_join(gwas0, gwas1, by="SNP") %>%
	select(one_of(c("SNP", "P.x", "P.y")))
colnames(gwas) = c("SNP", "Pnull", "Pdef")

gwasplot = arrange(gwas, Pnull) %>%
	mutate(Punif = sort(runif(nrow(gwas))))

ggplot(gwasplot[1:1000,]) + geom_point(aes(x=-log(Punif, 10), y=-log(Pnull, 10))) +
	geom_abline(slope=1, intercept=0, col="red")

gwasplot = arrange(gwas, Pdef) %>%
	mutate(Punif = sort(runif(nrow(gwas))))

ggplot(gwasplot[1:1000,]) + geom_point(aes(x=-log(Punif, 10), y=-log(Pdef, 10))) +
	geom_abline(slope=1, intercept=0, col="red")

ggplot(filter(gwasplot, Pnull<0.05 | Pdef<0.05)) +
	geom_point(aes(x=-log(Pnull, 10), y=-log(Pdef, 10))) +
	geom_abline(slope=1, intercept=0, col="red")
