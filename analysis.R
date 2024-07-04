library(ggplot2)
theme_set(theme_bw(base_size=20))

dat <- read.csv("results/results.csv", header=T)

#defaults
Drate <- 0.5
Lrate <- 0.5
Rrate <- 0.5

dat$allD <- dat$correctD+dat$locationD+dat$typeD+dat$wrongD
dat$propD <- dat$correctD/dat$allD
dat$proplocationD <- (dat$correctD+dat$locationD)/dat$allD
dat$proptypeD <- (dat$correctD+dat$typeD)/dat$allD

dat$allS <- dat$correctS+dat$locationS+dat$typeS+dat$wrongS
dat$propS <- dat$correctS/dat$allS
dat$proplocationS <- (dat$correctS+dat$locationS)/dat$allS
dat$proptypeS <- (dat$correctS+dat$typeS)/dat$allS

dat$allR <- dat$correctR+dat$locationR+dat$typeR+dat$wrongR
dat$propR <- dat$correctR/dat$allR
dat$proplocationR <- (dat$correctR+dat$locationR)/dat$allR
dat$proptypeR <- (dat$correctR+dat$typeR)/dat$allR

#D/L
DLdat <- dat[dat$Rrate == Rrate,]

DLtotals <- aggregate(DLdat, by = list(DLdat$Drate), FUN = sum, na.rm=T)
DLmeans <- aggregate(DLdat, by = list(DLdat$Drate), FUN = mean, na.rm=T)
DLsds <- aggregate(DLdat, by = list(DLdat$Drate), FUN = sd, na.rm=T)

#plot(DLmeans$Drate, DLmeans$propD)
#points(DLmeans$Drate, DLmeans$proplocationD,col="red")
#points(DLmeans$Drate, DLmeans$proptypeD,col="blue")
pdf("alg_perecon/figures/DL-D.pdf")
ggplot(data=DLmeans) + geom_point(aes(Drate,propD,col="red")) + 
	geom_point(aes(Drate,proplocationD, col="darkgreen")) + geom_point(aes(Drate,proptypeD,col="blue")) + 
	geom_errorbar(aes(Drate,ymin=propD-2/sqrt(1000)*DLsds$propD,ymax=propD+2/sqrt(1000)*DLsds$propD),col="red",width=0.02) +
	geom_errorbar(aes(Drate,ymin=proplocationD-2/sqrt(1000)*DLsds$proplocationD,ymax=proplocationD+2/sqrt(1000)*DLsds$proplocationD),col="darkgreen",width=0.02) +
	geom_errorbar(aes(Drate,ymin=proptypeD-2/sqrt(1000)*DLsds$proptypeD,ymax=proptypeD+2/sqrt(1000)*DLsds$proptypeD),col="blue",width=0.02) +
	xlab("D") + ylab("Correct proportion") +
	scale_colour_manual(values = c('red'='red','blue'='blue','darkgreen'='darkgreen'), labels = c('Correct type and location','Correct type','Correct location')) +
	theme(legend.title=element_blank(), legend.position = c(.05, .05), legend.justification = c("left", "bottom"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6))
dev.off()

#plot(DLmeans$Drate, DLmeans$propS)
#points(DLmeans$Drate, DLmeans$proplocationS,col="red")
#points(DLmeans$Drate, DLmeans$proptypeS,col="blue")
pdf("alg_perecon/figures/DL-S.pdf")
ggplot(data=DLmeans) + geom_point(aes(Drate,propS,col="red")) + 
	geom_point(aes(Drate,proplocationS, col="darkgreen")) + geom_point(aes(Drate,proptypeS,col="blue")) + 
	geom_errorbar(aes(Drate,ymin=propS-2/sqrt(1000)*DLsds$propS,ymax=propS+2/sqrt(1000)*DLsds$propS),col="red",width=0.02) +
	geom_errorbar(aes(Drate,ymin=proplocationS-2/sqrt(1000)*DLsds$proplocationS,ymax=proplocationS+2/sqrt(1000)*DLsds$proplocationS),col="darkgreen",width=0.02) +
	geom_errorbar(aes(Drate,ymin=proptypeS-2/sqrt(1000)*DLsds$proptypeS,ymax=proptypeS+2/sqrt(1000)*DLsds$proptypeS),col="blue",width=0.02) +
	xlab("D") + ylab("Correct proportion") +
	scale_colour_manual(values = c('red'='red','blue'='blue','darkgreen'='darkgreen'), labels = c('Correct type and location','Correct type','Correct location')) +
	theme(legend.title=element_blank(), legend.position = c(.05, .05), legend.justification = c("left", "bottom"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	coord_cartesian(ylim = c(0.9985,1.0005))
dev.off()

#plot(DLmeans$Drate, DLmeans$propR)
pdf("alg_perecon/figures/DL-R.pdf")
ggplot(data=DLmeans) + geom_point(aes(Drate,propR),col="red") + 
	geom_errorbar(aes(Drate,ymin=propR-2/sqrt(1000)*DLsds$propR,ymax=propR+2/sqrt(1000)*DLsds$propR),col="red",width=0.02) +
	xlab("D") + ylab("Correct proportion")
dev.off()

#plot(DLmeans$Drate, DLmeans$correctParalogy)
pdf("alg_perecon/figures/DL-paralogy.pdf")
ggplot(data=DLmeans) + geom_point(aes(Drate,correctParalogy)) + 
	geom_errorbar(aes(Drate,ymin=correctParalogy-2/sqrt(1000)*DLsds$correctParalogy,ymax=correctParalogy+2/sqrt(1000)*DLsds$correctParalogy),col="black",width=0.02) +
	xlab("D") + ylab("Correct paralogy proportion")
dev.off()

#D/R
DRdat <- dat[dat$Lrate == Lrate,]

DRtotals <- aggregate(DRdat, by = list(DRdat$Drate), FUN = sum, na.rm=T)
DRmeans <- aggregate(DRdat, by = list(DRdat$Drate), FUN = mean, na.rm=T)
DRsds <- aggregate(DRdat, by = list(DRdat$Drate), FUN = sd, na.rm=T)

#plot(DRmeans$Drate, DRmeans$propD)
#points(DRmeans$Drate, DRmeans$proplocationD,col="red")
#points(DRmeans$Drate, DRmeans$proptypeD,col="blue")
pdf("alg_perecon/figures/DR-D.pdf")
ggplot(data=DRmeans) + geom_point(aes(Drate,propD,col="red")) + 
	geom_point(aes(Drate,proplocationD, col="darkgreen")) + geom_point(aes(Drate,proptypeD,col="blue")) + 
	geom_errorbar(aes(Drate,ymin=propD-2/sqrt(1000)*DRsds$propD,ymax=propD+2/sqrt(1000)*DRsds$propD),col="red",width=0.02) +
	geom_errorbar(aes(Drate,ymin=proplocationD-2/sqrt(1000)*DRsds$proplocationD,ymax=proplocationD+2/sqrt(1000)*DRsds$proplocationD),col="darkgreen",width=0.02) +
	geom_errorbar(aes(Drate,ymin=proptypeD-2/sqrt(1000)*DRsds$proptypeD,ymax=proptypeD+2/sqrt(1000)*DRsds$proptypeD),col="blue",width=0.02) +
	xlab("D") + ylab("Correct proportion") +
	scale_colour_manual(values = c('red'='red','blue'='blue','darkgreen'='darkgreen'), labels = c('Correct type and location','Correct type','Correct location')) +
	theme(legend.title=element_blank(), legend.position = c(.95, .05), legend.justification = c("right", "bottom"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	coord_cartesian(ylim = c(0.9,1))	
dev.off()

#plot(DRmeans$Drate, DRmeans$propS)
#points(DRmeans$Drate, DRmeans$proplocationS,col="red")
#points(DRmeans$Drate, DRmeans$proptypeS,col="blue")
pdf("alg_perecon/figures/DR-S.pdf")
ggplot(data=DRmeans) + geom_point(aes(Drate,propS,col="red")) + 
	geom_point(aes(Drate,proplocationS, col="darkgreen")) + geom_point(aes(Drate,proptypeS,col="blue")) + 
	geom_errorbar(aes(Drate,ymin=propS-2/sqrt(1000)*DRsds$propS,ymax=propS+2/sqrt(1000)*DRsds$propS),col="red",width=0.02) +
	geom_errorbar(aes(Drate,ymin=proplocationS-2/sqrt(1000)*DRsds$proplocationS,ymax=proplocationS+2/sqrt(1000)*DRsds$proplocationS),col="darkgreen",width=0.02) +
	geom_errorbar(aes(Drate,ymin=proptypeS-2/sqrt(1000)*DRsds$proptypeS,ymax=proptypeS+2/sqrt(1000)*DRsds$proptypeS),col="blue",width=0.02) +
	xlab("D") + ylab("Correct proportion") +
	scale_colour_manual(values = c('red'='red','blue'='blue','darkgreen'='darkgreen'), labels = c('Correct type and location','Correct type','Correct location')) +
	theme(legend.title=element_blank(), legend.position = c(.05, .05), legend.justification = c("left", "bottom"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	coord_cartesian(ylim = c(0.993,1.0005))
dev.off()

#plot(DRmeans$Drate, DRmeans$propR)
pdf("alg_perecon/figures/DR-R.pdf")
ggplot(data=DRmeans) + geom_point(aes(Drate,propR),col="red") + 
	geom_errorbar(aes(Drate,ymin=propR-2/sqrt(1000)*DRsds$propR,ymax=propR+2/sqrt(1000)*DRsds$propR),col="red",width=0.02) +
	xlab("D") + ylab("Correct proportion") + coord_cartesian(ylim = c(0.895,1))	
dev.off()

#plot(DRmeans$Drate, DRmeans$correctParalogy)
pdf("alg_perecon/figures/DR-paralogy.pdf")
ggplot(data=DRmeans) + geom_point(aes(Drate,correctParalogy)) + 
	geom_errorbar(aes(Drate,ymin=correctParalogy-2/sqrt(1000)*DRsds$correctParalogy,ymax=correctParalogy+2/sqrt(1000)*DRsds$correctParalogy),col="black",width=0.02) +
	xlab("D") + ylab("Correct paralogy proportion")
dev.off()
