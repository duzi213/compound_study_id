
library(gsDesign)
#install.packages('gsDesign')
#install.packages('pwr')
library(pwr)

x <- getwd()

######

pwr.2p.test(h = -0.15, n = NULL, sig.level = 0.05, power = 0.8, alternative = "greater")

pwr.2p.test(h = 0.15, n = NULL, sig.level = 0.05, power = 0.8,alternative = "greater")
pwr.2p.test(h = ES.h(p1 = 0.55, p2 = 0.50), sig.level = 0.05, power = .80)


# design a 4-analysis trial using a Hwang-Shih-DeCani spending function 
# for both lower and upper bounds 
x <- gsDesign(k=4, sfu=sfHSD, sfupar=-2, sfl=sfHSD, sflpar=1)

#######events calculation based on a ppt

(1.645+1.036)^2/(0.5*0.5*(log(0.7))^2)
(1.645+0.842)^2/(0.5*0.5*(log(0.6))^2)
(1.645+1.01)^2/(0.5*0.5*(log(0.6))^2)
(1.645+1.01)^2/(0.5*0.5*(log(0.7))^2)
(1.645+1.01)^2/(0.5*0.5*(log(0.5))^2)

##non-inferiority
pA=0.96
pB=0.87
delta=0.08
kappa=2
alpha=0.05
beta=0.20
(nB=(pA*(1-pA)/kappa+pB*(1-pB))*((qnorm(1-alpha)+qnorm(1-beta))/(pA-pB-delta))^2)
ceiling(nB) # 25
z=(pA-pB-delta)/sqrt(pA*(1-pA)/nB/kappa+pB*(1-pB)/nB)
(Power=pnorm(z-qnorm(1-alpha))+pnorm(-z-qnorm(1-alpha)))


alpha <- 0.05
belta <- 0.2
f <- (qnorm(1-alpha)+qnorm(1-belta))^2
pi1 <- 0.87
pi2 <- 0.96
margin <- 0.08
n1 <- 2*f*(pi1*(1-pi1)+pi2*(1-pi2))/((pi1-pi2-margin)^2)


alpha <- 0.05
belta <- 0.2
f <- (qnorm(1-alpha)+qnorm(1-belta))^2
pi1 <- 0.17
pi2 <- 0.14
margin <- 0.2
nn <- 2*f*(pi1*(1-pi1)+pi2*(1-pi2))/((pi1-pi2-margin)^2)

alpha <- 0.05
belta <- 0.2
f <- (qnorm(1-alpha)+qnorm(1-belta))^2
pi1 <- 87
pi2 <- 96
margin <- 7
n <- f*(pi1*(100-pi1)+pi2*(100-pi2))/((pi1-pi2-margin)^2)


alpha <- 0.05
belta <- 0.2
f <- (qnorm(1-alpha)+qnorm(1-belta))^2
pi1 <- 87
pi2 <- 96
dropout <- 0.15
margin <- seq(from=7, to= 8.3, by= 0.1)
ne <- ceiling( f*(2*pi1*(100-pi1)+pi2*(100-pi2))/((pi1-pi2-margin)^2))
nc <- ceiling(ne/2)
n_total <- ne+nc
n_adj <- ceiling(n_total/(1-dropout))


plot(margin,n_adj,ylab='total pts')


alpha <- 0.05
belta <- 0.2
f <- (qnorm(1-alpha)+qnorm(1-belta))^2
pi1 <- 87
pi2 <- 96
dropout <- 0.15
margin <-7
ne <- ceiling( f*(2*pi1*(100-pi1)+pi2*(100-pi2))/((pi1-pi2-margin)^2))
nc <- ceiling(ne/2)
n_total <- ne+nc
n_adj <- ceiling(n_total/(1-dropout))




alpha <- 0.025
belta <- 0.1
f <- (qnorm(1-alpha)+qnorm(1-belta))^2
pi1 <- 89  ##control
pi2 <- 87 ##experimental  DTG+3TC
margin <- 10
n2 <- f*(pi1*(100-pi1)+pi2*(100-pi2))/((pi1-pi2-margin)^2)
n2

pa <- 0.96
pb <- 0.87
delta=-0.07
kappa=2
alpha=0.05
beta=0.2
Nc=ceiling((pa*(1-pa)/kappa+pb*(1-pb))*f/((pa-pb-delta)^2))



##########way 2
pa <- 0.96
pb <- 0.87
delta=-0.07
kappa=2
alpha=0.05
beta=0.2
nb=ceiling((pa*(1-pa)/kappa+pb*(1-pb))*f/((pa-pb-delta)^2))

##kappa= na/nb

plot(margin,n_adj,ylab='total pts')



help(pwr.2p.test)
# print the design
x

# since sfHSD is the default for both sfu and sfl,
# this could have been written as
x <- gsDesign(k=4, sfupar=-2, sflpar=1)

# print again
x

# plot the spending function using many points to obtain a smooth curve
# show default values of gamma to see how the spending function changes
# also show gamma=1 which is supposed to approximate a Pocock design
t <- 0:100/100
plot(t,  sfHSD(0.025, t, -4)$spend,
     xlab="Proportion of final sample size", 
     ylab="Cumulative Type I error spending", 
     main="Hwang-Shih-DeCani Spending Function Example", type="l")
lines(t, sfHSD(0.025, t, -2)$spend, lty=2)
lines(t, sfHSD(0.025, t, 1)$spend, lty=3)
legend(x=c(.0, .375), y=.025*c(.8, 1), lty=1:3, 
       legend=c("gamma= -4", "gamma= -2", "gamma= 1"))



##############12/14 #################

#####
w <- -log(.5) / 6
ww <- -log(.5) / 6*0.7

x <- nSurvival(lambda1=-log(.5) / 6, lambda2=-log(.5) / 6*0.7,alpha=0.025,sided=1,
               eta=-log(.95)/12, Tr=30 ,Ts=36, type="rr", entry="unif")

x <- nSurvival(lambda1=w, lambda2=ww,alpha=0.025,sided=1,ratio =1,
               eta=-log(.95)/12, Tr=30 ,Ts=36, type="rr", entry="unif")


x$Sample.size
x$Num.events


######################12.20 #######################
######  Power function #########################
#The power function for a one-sample z-test can be calculated using R.

pow.z.test <- function(alpha,mu1,mu0,sigma,n){
  arg1 <- qnorm(1-alpha/2)-(mu1-mu0)/(sigma/sqrt(n))
  arg2 <- -1*qnorm(1-alpha/2)-(mu1-mu0)/(sigma/sqrt(n))
  1-pnorm(arg1)+pnorm(arg2)
}

power.prop.test(p1 = 0.5,p2 = 0.7,power = 0.817,sig.level=0.05,alternative=c("one.sided"))
power.prop.test(p1 = 0.46,p2 = 0.21,power = 0.802,sig.level=0.1,alternative=c("one.sided"))


####12.29 another way for sample size for binomial endpoint

n.fix <- FarrMannSS(p1=.5, p2=.7, beta=.183, outtype=1)


##################1.5.2020######################

library(powerSurvEpi)
library(survival)
data(Oph)
res <- powerCT(formula = Surv(times, status) ~ group, dat = Oph,
               nE = 200, nC = 200, RR = 0.7, alpha = 0.05)
# Table 14.24 on page 809 of Rosner (2006)
print(round(res$mat.lambda, 4))
# Table 14.12 on page 787 of Rosner (2006)
print(round(res$mat.event, 4))
# the power
print(round(res$power, 2))

##############another package PropCIs used in 081##################
install.packages('PropCIs')
help(diffscoreci)

##############events calcuation##################

#e=4*[z(2/alpha)+zbeta]^2/[ln(hr)]^2

##create functon for events calculation 

# pow <- function(x, y) {
#   # function to print x raised to the power y
#   result <- x^y
#   print(paste(x,"raised to the power", y, "is", result))
# }

events <- (qnorm(1-0.025)+qnorm(0.929))^2*4/((log(0.6))^2)
events <- (qnorm(1-0.025)+qnorm(0.989))^2*4/((log(0.55))^2)
events <- (qnorm(1-0.0125)+qnorm(0.925))^2*4/((log(0.7))^2)
events
