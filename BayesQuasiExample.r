lung <- read.table("http://faculty.washington.edu/jonno/book/MNlung.txt", header=TRUE, sep="\t")
radon <- read.table("http://faculty.washington.edu/jonno/book/MNradon.txt", header=TRUE)
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
	rad.avg[i]<-mean(radon[radon$county==i,2])
}
x <- rad.avg
rad.avg[26] <- 0
rad.avg[63] <- 0
x[26] <- NA
x[63] <- NA
#
# Fig 2.1
#
#pdf("lungcancerlogass.pdf",width=4,height=3.5)
plot(log(Obs/Exp)~x,xlab="Average Radon (pCi/liter)",ylab="Log SMR")
lines(lowess(log(Obs/Exp)[is.na(x)==F]~x[is.na(x)==F]))
#dev.off()
#
# Lung cancer MLE example in Section 2.4
#
poismod <- glm(Obs~offset(log(Exp))+x,family="poisson",na.action=na.omit)
quasimod <- glm(Obs~offset(log(Exp))+x,family=quasipoisson(link=log),
            na.action=na.omit)
coef(summary( poismod))[,2] # estimated standard errors of beta-hats
coef(summary(quasimod))[,2] # estimated standard errors of beta-hats
summary(quasimod, corr=TRUE)$corr # estimated correlation of beta-hats
confint.default(quasimod)

useme <- !is.na(x) ## nonmissing entries
#where the alpha estimate comes from:
sum(resid(poismod, type="response")^2/fitted(poismod))/(sum(!is.na(x))-2)

#for comparison, the version we use in the Bayes analog:
sum(resid(poismod, type="response")^2/fitted(poismod))/(sum(!is.na(x)))



poismod$coeff[2]-1.96*sqrt(vcov(poismod)[2,2]) # 95%
poismod$coeff[2]+1.96*sqrt(vcov(poismod)[2,2]) # CI for beta1
exp(quasimod$coeff[2]-1.96*sqrt(vcov(poismod)[2,2]))  # 95%
exp(quasimod$coeff[2]+1.96*sqrt(vcov(poismod)[2,2]))  # CI for exp(beta1)
#
# Lung cancer QMLEs in Section 2.5
#
quasimod <- glm(Obs~offset(log(Exp))+x,family=quasipoisson(link="log"))
exp(quasimod$coeff[2]-1.96*sqrt(vcov(quasimod)[2,2]))  # 95%
exp(quasimod$coeff[2]+1.96*sqrt(vcov(quasimod)[2,2]))  # CI for exp(beta1)

mean(x, na.rm=TRUE)
summary( glm(Obs~offset(log(Exp))+x,family=quasipoisson(link=log), na.action=na.omit) , corr=TRUE)
summary( glm(Obs~offset(log(Exp))+I(x-5.288),family=quasipoisson(link=log), na.action=na.omit) , corr=TRUE)

setwd("C:/Users/kenrice/Desktop/intervalset")
# may need to run as administrator, to write/execute external files

# set up an external model file, using cat()
cat(file="qpoissonexample.txt", "model{ 
for(i in 1:n){
	y[i] ~ dpois(mu[i]) 
	mu[i] <- exp(log(Exp[i]) + beta[1] + beta[2]*x[i])
	pearson.resid2[i] <- (y[i]-mu[i])^2/mu[i]
	}
	alpha <- mean(pearson.resid2[])
	beta[1]~dnorm(0,0.1)
	beta[2]~dnorm(0,0.1)
} ")

library("rjags")
useme <- complete.cases(cbind(Obs,Exp,x))
jags1 <- jags.model("qpoissonexample.txt", data=list(y=Obs[useme], Exp=Exp[useme], x=x[useme],n=sum(useme)) )
#jags1 <- jags.model("qpoissonexample.txt", data=list(y=Obs[useme], Exp=Exp[useme], x=x[useme]-mean(x[useme]),n=sum(useme)) )
update(jags1, 10000)
store1 <- coda.samples(jags1, c("beta","alpha"), n.iter=100000) 
summary( store1 )
#plot( store1 )

confint.default(poismod)[2,]
confint.default(quasimod)[2,]
#confint.default( glm(Obs~offset(log(Exp))+I(x-5.288),family=quasipoisson(link=log), na.action=na.omit) , corr=TRUE)

mean(as.mcmc(store1)[,"beta[2]"])+ c(-1,1)*1.96*sd(as.mcmc(store1)[,"beta[2]"])
mean(as.mcmc(store1)[,"beta[2]"])+ sqrt(mean(as.mcmc(store1)[,"alpha"]))*c(-1,1)*1.96*sd(as.mcmc(store1)[,"beta[2]"])

