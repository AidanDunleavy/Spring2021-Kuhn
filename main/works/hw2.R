#' ---
#' title: HW 2
#' author: Aidan Dunleavy
#' date: 01/26/2021
#' ---
library(AlgDesign)
library(BsMD)
library(car)
library(DoE.base)
library(FrF2)
library(GAD)
library(gmodels)
library(leaps)
library(lme4)
library(lsmeans)
library(mixexp)
library(multcomp)
library(Vdgraph)
library(pracma)
library(zip)
library(later)
library(htmltools)
library(digest)
library(promises)
library(classInt)
library(units)
library(agricolae)
library(daewr)
#' PROBLEM 1
#' 
#' (a)
D <- expand.grid(Beat = X3_10_1$Beat[1:3],
                 Course = X3_10_1$Course[1:3])
D[] <- lapply(D, factor)
Score = X3_10_1$Score
ScoreData = data.frame(cbind(D,Score))
ScoreData
with(ScoreData, 
     (interaction.plot(
       Course, 
       Beat, 
       Score, 
       type = "b", 
       pch = c(18,24,22),
       leg.bty = "o",
       main = "Interaction Plot of Course length and Beat",
       xlab = "Course Length",
       ylab = "Score"))
)
help("interaction.plot")
#' PROBLEM 5

y <- flight.time <- c(
        5.2, 5.1, 5.3, 5.3, 5.2, 5.2, 5.0, 5.2,
        5.2, 5.1, 5.1, 5.0, 5.2, 5.1, 5.1, 5.2,
        5.1, 5.1, 5.0, 5.0, 5.1, 5.0, 4.9, 5.1,
        4.8, 4.9, 4.8, 5.0, 4.8, 4.9, 4.9, 4.9
); y; class(y)
x <- wing.length <- factor(c(rep(4,8), rep(4.75, 8), rep(5.5,8), rep(6,8))); x; class(x)
helicopter <- as.numeric(1:32); helicopter; class(helicopter)
paper.heli <- data.frame(helicopter, wing.length, flight.time); paper.heli
sd(y)
library(daewr)
rmin <- 2 # smallest number of replicates
rmax <- 2 # largest number of replicates
alpha <- .05
sigma <- 0.32
Delta <- 1.0
nlev <- c(4,4)
nreps <- c(rmin:rmax)
result <- Fpower2(alpha, nlev, nreps, Delta, sigma)
options(digits = 5)
result
help(Fpower1)

if (is.null(alpha)|is.null(nlev)|is.null(nreps)|is.null(Delta)|is.null(sigma))
        stop("you must supply alpha, nlev, nreps, Delta and sigma")
if(length(nlev)<2)
        stop ("nlev must be a two component vecto containing levels of the 1st and 2nd factors")
a <- nlev[1]
b <- nlev[2]
cssb <- (Delta^2)/2
ncb <- a*(nreps*cssb)/(sigma^2)
cssa<-(Delta^2)/2
nca<- b*(nreps*cssa)/(sigma^2)
dfa<- a-1
dfb<- b-1
df2<-(nreps-1)*b*a
powera <- 1-pf(Fcrit(alpha,dfa,df2),dfa,df2,nca)
powerb <- 1-pf(Fcrit(alpha,dfb,df2),dfa,df2,nca)
result <-cbind(nreps,df2,powera,powerb)
result

library(daewr) # 2k design, 2 replicates
# Table 3
Des <- expand.grid( F = c(-1, 1), # 2^3 design
                    C = c(-1, 1), # 2 x 2 x 2  = 8 
                    D = c(-1, 1)) # 8 x 2 replicates = 16

y <- c(10.50, 2.01, 4.16, 1.58, 5.26, 1.58, 4.10, 1.54, 
       10.65, 1.59, 4.14, 1.57, 5.27, 1.54, 4.05, 1.48); y
DesY <- data.frame(cbind(Des,y)); DesY 
m.water <- lm( y ~ F*C*D, data = DesY) # model m.water
summary(m.water)
m.water <- aov( y ~ F*C*D, data = DesY) # model m.water
summary(m.water) # notice estimate = 0.5 effect, except intercept
# Table 4
Des.inv <- expand.grid( F = c(-1, 1), # 2^3 design
                        C = c(-1, 1), # 2 x 2 x 2  = 8 
                        D = c(-1, 1)) # 8 x 2 replicates = 16
y.inv <- 1/y; y.inv
DesY.inv <- data.frame(cbind(Des.inv,y.inv)); DesY.inv 
m.water.inv <- lm( y.inv ~ F*C*D, data = DesY.inv) # model m.water.inverse
summary(m.water.inv) # notice estimate = 0.5 effect, except intercept

flame=rep(c("low","high"),4)
cover=rep(c("no","no","yes","yes"),2)
dish=c("small","small","small","small","large","large","large","large")
data = c(10.50, 2.01, 4.16, 1.58, 5.26, 1.58, 4.10, 1.54, 
         10.65, 1.59, 4.14, 1.57, 5.27, 1.54, 4.05, 1.48)
dataset=data.frame(cbind(dish,flame,cover,data))
m<-aov(data=dataset,data~dish*flame*cover)
summary(m)


#' PROBLEM 9
D <- expand.grid(silane = c(0.1,0.9),
                 gas = c(40,220),
                 pressure = c(300,1200),
                 temp = c(300,460),
                 power = c(10,60))
D[] <- lapply(D, factor)

D
Refract = c(1.92, 3.06, 1.96, 3.33, 1.87, 2.62, 1.97, 2.96, 1.94, 3.53, 2.06, 3.75, 1.96, 3.14, 2.15, 3.43, 1.95, 3.16, 2.01, 3.43, 1.88, 2.14, 1.98, 2.81, 1.97, 3.67, 2.09, 3.73, 1.98, 2.99, 2.19, 3.39)
refractdata = data.frame(cbind(D,Refract))
m <- aov(data = refractdata, Refract ~ silane*gas*pressure*temp*power)
plot(m)
model = lm(data = refractdata, Refract ~ silane*gas*pressure*temp*power)
model
plot(refractdata$Refract)

m <- aov(data = refractdata, Refract ~ gas*pressure*temp*power)
summary(m)
plot(m)

