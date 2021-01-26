# Lecture Notes 2
# Chapter 3 Factorial Designs

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

# 3.1 Introduction 
# 3.2 Classical One at a Time versus Factorial Plans
# 3.3 Interpreting Interactions

D <- expand.grid( temp = c(0, 10, 20, 30), # design
                  noise = c(60, 90, 120) ) # combinations of factors
D <- rbind(D, D); D # row-bind twice, for 2 replicates

set.seed(2591) # ensures everyone has same randomization
D <- D[order(sample(1:24)), ] # randomize mice to temp x noise treatments
ROCDes <- D[c( "temp", "noise" )]; ROCDes # label columns of factors
write.csv(ROCDes, file = "ROCDes.csv", row.names = FALSE) # save

D <- expand.grid( temp = c(0, 10, 20, 30), # design
                  noise = c(60, 90, 120) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
D <- rbind(D, D); D # row-bind twice, for 2 replicates
D[] <- lapply(D, factor) # convert to factors
ROC <- c(
  10.3, 1.8, 1.2, 12.4, 
  9.1, 12.1, 6.5, 16.1, 
  6.1, 5.1, 1.2, 18.1, 
  7.2, 9.8, 8.1, 15.1,
  5.4, 4.2, 4.1, 17.2,
  2.1, 6.2, 2.1, 19.1); 
ROCdata <- data.frame(cbind(D,ROC)); ROCdata # data.frame

m1 <- aov(ROC ~ temp * noise, data = ROCdata)
summary(m1) # ANOVA temp x noise; remember to factorize temp, noise

model.tables(m1, type = "means", se = T) # various means

ct1 <- c(-3, 1, 1, 1) # temp contrast tested: 1st = avg(2,3,4)?
ct2 <- c(0, 2, -1, -1) # contrast for ortho, j-1 = 4 - 1 = 3 conditions
ct3 <- c(0, 0, 1, -1) # dummy contrast for ortho, 3 conditions
dot(ct1,ct2); dot(ct1,ct3); dot(ct2,ct3) # all ortho because dot prod zero
ctm <- cbind( ct1, ct2, ct3 ) # temp contrast matrix

cn1 <- c(-1, 0, 1) # noise contrast tested: 1st = 3rd?
cn2 <- c(1, -2, 1) # dummy contrast for ortho and j-1 = 2 - 1 = 1 
dot(cn1,cn2) # dot product zero because c1 orthog to c2
cnm <- cbind( cn1, cn2 ) # noise contrast matrix

m2 <- aov(ROC ~ temp * noise, 
          contrasts = list( temp = ctm, noise = cnm ), 
          data = ROCdata); summary(m2)
library(gmodels)
# 12 coefficients: mu, t1, t2, t3, n1, n2,
# 12 coefficients: tn11, tn12, tn21, tn22, tn31, tn32
c <- rbind('temp 10 vs ave(20,30,40)' 
           = c(0,1,0,0,0,0,0,0,0,0,0,0), 
           'noise 60 vs 120' 
           = c(0,0,0,0,1,0,0,0,0,0,0,0) )
estimable(m2,c)

# interaction plots of ROC on temp x noise
par( mfrow = c(1,2) )
with(ROCdata, 
     (interaction.plot(
        temp, 
        noise, 
        ROC, 
        type = "b", 
        pch = c(18,24,22,21), 
        leg.bty = "o",
        main = "Interaction Plot of Temperature and Noise",
        xlab = "Temperature",
        ylab = "ROC"))
     )

temp <- ROCdata$temp
with(ROCdata, 
     (interaction.plot(
       noise, 
       temp, 
       ROC, 
       type = "b", 
       pch = c(18,24,22), 
       leg.bty = "o",
       main = "Interaction Plot of Temperature and Noise",
       xlab = "Noise",
       ylab = "ROC"))
)
par(mfrow = c(1,1))

# number of replicates
library(daewr) # 2-factor treated as one factor
rmin <- 12 # smallest number of replicates
rmax <- 20 # largest number of replicates
sigma <- sqrt(9.61) # sqrt of msE
alpha <- .05
Delta <- 5 # ROC differs by 5
nlev <- 12 # 4 temp x 3 noise =12
nreps <- c(rmin:rmax)
power <- Fpower1(alpha, nlev, nreps, Delta, sigma)
options(digits = 5)
power

library(daewr) # 2-factor treated as 2-factor
rmin <- 2 # smallest number of replicates
rmax <- 5 # largest number of replicates
alpha <- .05
sigma <- sqrt(9.61) # sqrt of msE
Delta <- 5 # ROC differs by 5
nlev <- c(4,3) # 4 temp x 3 noise
nreps <- c(rmin:rmax)
power2 <- Fpower2(alpha, nlev, nreps, Delta, sigma)
options(digits = 5)
power2

ROCdatam <- ROCdata # unbalanced design
ROCdatam[24, 3] <- NA; ROCdatam

library(car) # balanced case
m2 <- lm( ROC ~ temp*noise, 
          data = ROCdata, 
          contrasts = list( 
            temp = contr.sum, # makes sure contrasts sum to zero
            noise = contr.sum )
)
Anova( m2, type="II" )

library(car) # unbalanced case
m2 <- lm( ROC ~ temp*noise, 
            data = ROCdatam, 
            contrasts = list( 
              temp = contr.sum, # makes sure contrasts sum to zero
              noise = contr.sum )
            )
Anova( m2, type="II" )

library(lsmeans)
lsmeans(m2,~ temp) # marginal adjusted means for temp
lsmeans(m2,~noise) # marginal adjusted means for noise


library(daewr) # one-replication per cell with numeric factors
ROCcellmeans <- tapply( # two replicates averaged to one per cell
  ROCdata$ROC, 
  list(ROCdata$temp, ROCdata$noise), mean )
dim(ROCcellmeans) <- NULL
D <- expand.grid( temp = c(0, 10, 20, 30), # design
                  noise = c(60, 90, 120) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
ROCcells <- data.frame( D, ROCcellmeans )
mnr <- lm(ROCcellmeans ~ temp*noise, data=ROCcells )
anova(mnr)

tempc <- as.ordered(ROCcells$temp) # order factor levels
noisec <- as.ordered(ROCcells$noise)
tempLin <- contr.poly(tempc)[tempc,".L"] # create residuals from
noiseLin <- contr.poly(noisec)[noisec,".L"] # linear x linear partion of interaction
mROCbo <-lm(ROCcellmeans ~ tempc + noisec + tempLin:noiseLin, data=ROCcells)
anova(mROCbo)

ROCpred <- predict(mROCbo, # predict new ROC from linear portion of temp, noise
    newdata = data.frame(tempc, noisec, tempLin, noiseLin))
pred.means.ROC <- aggregate( # combine into new "linear" data
  ROCpred, 
  by = list(tempc = tempc, noisec = noisec), 
  "mean")
temperature <- pred.means.ROC$tempc
interaction.plot( # interaction plot with "linear" data
  pred.means.ROC$noisec, 
  temperature, 
  pred.means.ROC$x, # "x" is default label for means
  type="b", 
  pch = c(18,24,22), 
  leg.bty ="o",
  xlab = "noise", 
  ylab = "predicted ROC")

library(daewr) # one-replication per cell with categorical factors
D <- expand.grid( temp = c("T1", "T2", "T3", "T4"), # design
                  noise = c("N1", "N2", "N3") ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
ROCcells.categorical <- data.frame( ROCcellmeans, D )
Tukey1df(ROCcells.categorical) 


# 3.4 Creating a Two-Factor Factorial Plain in R
# 3.5 Analysis of a Two-Factor Factorial in R
# 3.6 Factorial Designs with Multiple Factors--CRFD

D <- expand.grid( A = c(1, 2, 3, 4), # binomial, logit link
                  B = c(1, 2, 3), # 4 x 3 x 2 = 24 configurations
                  C = c(1, 2)) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
# as.integer(runif(24,800,900)) # model number of sump pumps tried 
sump.pump <- c( # 
  869, 890, 871, 892, 866, 833, 
  885, 898, 874, 825, 852, 868, 
  888, 841, 889, 869, 856, 887, 
  826, 826, 891, 883, 898, 891); 
# rpois(24,50) # model number of defects 
defects <- c( # 
  61, 65, 72, 69, 60, 62, 
  52, 45, 74, 53, 66, 56,
  64, 56, 65, 66, 44, 45, 
  52, 57, 49, 53, 56, 44); 
spdata <- data.frame(cbind(D,sump.pump,defects)); spdata # data.frame

mb <- glm( # glm: binomial distribution, logit link
  cbind(defects, sump.pump) ~ A * B * C, # 3-factor
  data = spdata, 
  family = binomial ) # mle fit of model
summary(mb)
anova(update(mb, # ANOVA type III sums of squares
  .~ A + B + C + A:B + A:C + B:C + A:B:C ), 
  test = "Chisq") # chi-square tests

sp.prop <- spdata$defects / spdata$sump.pump # interaction plots
spdata.p <- data.frame(spdata,sp.prop); spdata.p
par (mfrow = c(1,2) )
spdata.p1 <- subset(spdata.p, C == 1)
interaction.plot(
  spdata.p1$A, 
  spdata.p1$B, 
  spdata.p1$sp.prop, type = "l",
  legend=FALSE, 
  ylim = c(.04,.09), 
  main = "Water Hard = 0.5", 
  xlab = "Water Temperature", 
  ylab = "Proportion Defect Sump Pumps")
spdata.p2 <- subset(spdata.p, C == 2)
interaction.plot(
  spdata.p2$A, 
  spdata.p2$B, 
  spdata.p2$sp.prop, type = "l",
  legend=FALSE, 
  ylim = c(.04,.09), 
  main = "Water Hard = 0.8", 
  xlab = "Water Temperature", 
  ylab = "Proportion Defect Sump Pumps")
par (mfrow = c(1,1) )

# 3.7 Two-Level Factorials

D <- expand.grid( A = c(0, 30), # 2^3 design
                  B = c(60, 120), # 2 x 2 x 2 = 8 configurations
                  C = c(0.5, 0.8)) # 8 x 2 replicates = 16
D[] <- lapply(D, factor) # convert to factors

defects <- c( # 
  61, 69, 60, 45, 64, 66, 44, 57, 
  60, 61, 52, 41, 56, 61, 51, 63); 
twokdata <- data.frame(cbind(D,defects)); twokdata # data.frame

library(daewr)
library(FrF2) # 2^k factorial design
m2k <- lm( 
  defects ~ A*B*C, 
  data=twokdata, 
  contrast=list( # convert factor levels in +/-
    A=contr.FrF2, 
    B=contr.FrF2, 
    C=contr.FrF2))
summary(m2k)

C_Water.Hard=twokdata$C # interaction plot
with(twokdata, 
     (interaction.plot(
      A, C_Water.Hard, defects, 
      type = "b",
      pch = c(24,22),
      leg.bty = "o",
      xlab = "Water Temperature",
      ylab = "Defects")))

library(daewr) # 2k design, 1 replicate
Des <- expand.grid( A = c(-1, 1), # 2^4 design
                  B = c(-1, 1), # 2 x 2 x 2 x 2 = 16 configurations
                  C = c(-1, 1),
                  D = c(-1, 1)) # 16 x 1 replicates = 16
Des[] <- lapply(Des, factor) # convert to factors

defects <- c( # 
  61, 69, 60, 45, 64, 66, 44, 57, 
  60, 61, 52, 41, 56, 61, 51, 63); 
twok1rdata <- data.frame(cbind(Des,defects)); twok1rdata # data twok1r
m2k1r <- lm( defects ~ A*B*C*D, data = twok1rdata) # model m2k1r
summary(m2k1r)

par (mfrow = c(1,2) ) # plots to identify significant effects
fullnormal(coef(m2k1r)[-1],alpha=.025) # normal prob plot
LGB( coef(m2k1r)[-1], rpt = TRUE) # half-normal plot, with description
par (mfrow = c(1,1) )

library(daewr) # interaction plot
with(twok1rdata, 
 (interaction.plot( 
  A, B, defects, 
  type = "b", 
  pch = c(18,24), 
  main = "Interaction Plot of On-Offs by Water Temperature",
  xlab = "Water Temperature", 
  ylab = "Number Defects"))
 )

library(BsMD) # more plots to detect significant regression coefficients
par( mfrow = c(2,1) )
LenthPlot(m2k1r, main = "Lenth Plot of Effects") # lenth plot
X <- model.matrix(m2k1r)[ , 2:16] # 
defects <- twok1rdata$defects
twok1r.BsProb <- BsProb( 
  X = X, y = defects, 
  blk = 0, 
  mFac = 15, # 15 factors
  mInt = 1, 
  p = 0.2, 
  g = 2.49, 
  ng = 1, 
  nMod = 10)
plot( twok1r.BsProb, main = "Bayes Plot of Effects" )
par (mfrow = c(1,1) )

# 3.8 Verifying Assumptions of the Model

library(daewr) # 2k design, 1 replicate
Des <- expand.grid( A = c(-1, 1), # 2^4 design
                  B = c(-1, 1), # 2 x 2 x 2 x 2 = 16 configurations
                  C = c(-1, 1),
                  D = c(-1, 1)) # 16 x 1 replicates = 16

y <- defects <- c(61, 69, 60, 45, 64, 66, 44, 57, 
  60, 61, 52, 41, 56, 61, 51, 63); y
DesY <- data.frame(cbind(Des,y)); 
names(DesY)[5] <- "y"; DesY # data twok1r

m2k1r <- lm( y ~ A*B*C*D, data = DesY)
summary(m2k1r) # check for gaps
fullnormal(coef(m2k1r)[-1], alpha=.2)
data(DesY)
Gaptest(DesY) # check for outliers


# 3.9 Review of Important Concepts

#############################################
# Quiz 1 help with R code
#############################################

# 3.2.b pp 108-112
D <- expand.grid( stop = c(1, 2, 3), # stop 3 levels
                  start = c(1, 2, 3)) # start 3 levels
D <- rbind(D, D); D # row-bind twice, for 2 replicates
D[] <- lapply(D, factor); D # convert to factors
y <- distance <- c(5.4, 4.2, 6.6, 7.0, 4.4, 7.6, 4.9, 6.5, 4.2,
                   6.7, 4.3, 7.9, 6.6, 5.1, 6.8, 5.7, 4.7, 5.6); y; class(y)
catapult <- data.frame( D, y ); catapult
set.seed(2591) # ensures everyone has same randomization
D <- D[order(sample(1:18)), ]; D # randomize mice to stop x start treatments

# 3.2.c pp 108-112
rmin <- 2 # smallest number of replicates
rmax <- 6 # largest number of replicates
alpha <- .05
sigma <- sqrt(12) # sqrt of sigma^2 = 12
Delta <- 10 # detects 10 inches in cell means
nlev <- 9 # 3 stop x 3 start
nreps <- c(rmin:rmax)
power <- Fpower1(alpha, nlev, nreps, Delta, sigma)
options(digits = 5)
power

# 3.2.d pp 108-112
rmin <- 2 # smallest number of replicates
rmax <- 5 # largest number of replicates
alpha <- .05
sigma <- sqrt(12) # sqrt of sigma^2 = 12
Delta <- 24 # detects 24 inches in marginal means
nlev <- c(3,3) # 3 stop x 3 start
nreps <- c(rmin:rmax)
power2 <- Fpower2(alpha, nlev, nreps, Delta, sigma)
options(digits = 5)
power2

# 3.2.f pp 108-112
D <- expand.grid( stop = c(1, 2, 3), # stop 3 levels
                  start = c(1, 2, 3)) # start 3 levels
D <- rbind(D, D); D # row-bind twice, for 2 replicates
D[] <- lapply(D, factor); D # convert to factors
y <- distance <- c(5.4, 4.2, 6.6, 7.0, 4.4, 7.6, 4.9, 6.5, 4.2,
  6.7, 4.3, 7.9, 6.6, 5.1, 6.8, 5.7, 4.7, 5.6); y; class(y)
catapult <- data.frame( D, y ); catapult
mcat <- lm(y ~ start*stop, data=catapult )
anova(mcat)

# 3.2.g pp 108-112
D <- expand.grid( stop = c(1, 2, 3), # stop 3 levels
                  start = c(1, 2, 3)) # start 3 levels
D <- rbind(D, D); D # row-bind twice, for 2 replicates
D[] <- lapply(D, factor); D # convert to factors
y <- distance <- c(5.4, 4.2, 6.6, 7.0, 4.4, 7.6, 4.9, 6.5, 4.2,
                   6.7, 4.3, 7.9, 6.6, 5.1, 6.8, 5.7, 4.7, 5.6); y; class(y)
with(D, 
     (interaction.plot(
       stop, 
       start, 
       y, 
       type = "b", 
       pch = c(18,24,22,21), 
       leg.bty = "o",
       main = "Interaction Plot of Stop and Start",
       xlab = "Stop",
       ylab = "Distance (inches)"))
)

# 3.4.b pp 108-112
##### Power Calculation for three way ANOVA ###########
# Argument list
# alpha the significance level of the test.
# nlev vector containing the number of levels of the factors. 
# nreps the number of replicates in each combination of factor levels.
# Delta the size of a practical difference in three marginal factor level means.
# sigma the standard deviation of the experimental error.
############################################################
Fpower3 <- function(alpha=NULL, nlev=NULL, nreps=NULL, Delta=NULL, sigma=NULL)
{
  if (is.null(alpha)|is.null(nlev)|is.null(nreps)|is.null(Delta)|is.null(sigma))
    stop("you must supply alpha, nlev, nreps, Delta and sigma")
  if(length(nlev)<2)
    stop ("nlev must be a three component vector containing levels of 1st-3rd factors")
  a <- nlev[1]; a # brand
  b <- nlev[2]; b # power
  c <- nlev[3]; c # time
  cssa <-(Delta^2)/2; cssa
  nca <- b*c*(nreps*cssa)/(sigma^2); nca
  dfa <- a-1; dfa
  cssb <- (Delta^2)/2; cssb
  ncb <- a*c*(nreps*cssb)/(sigma^2); ncb
  dfb <- b-1; dfb
  cssc <-(Delta^2)/2; cssc
  ncc <- a*b*(nreps*cssa)/(sigma^2); ncc
  dfc <- c-1; dfc
  df2 <- (nreps-1)*b*a*c; df2
  powera <- 1-pf(Fcrit(alpha,dfa,df2),dfa,df2,nca); powera
  powerb <- 1-pf(Fcrit(alpha,dfb,df2),dfb,df2,ncb); powerb
  powerc <- 1-pf(Fcrit(alpha,dfc,df2),dfc,df2,ncc); powerc
  dat <- cbind(alpha,a,b,c,nreps,df2,sigma,powera,powerb,powerc); dat
  return(dat) 
}
rmin <- 2 # smallest number of replicates
rmax <- 6 # largest number of replicates
alpha <- .05
sigma <- 0.14 
Delta <- 0.25 # detects 0.25 = 25\% diff marginal proportions
nlev <- c(3,2,2) # 3 brand (a) x 2 time (b) x 2 power (c)
nreps <- c(rmin:rmax)
power3 <- Fpower3(alpha, nlev, nreps, Delta, sigma)
options(digits = 5)
power3

# 3.4.c pp 108-112
D <- expand.grid( brand = c(1, 2, 3), # brand 3 levels
                  power = c(1, 2), # power 2 levels
                  time = c(1, 2)) # time 2 levels
D <- rbind(D, D, D); D # row-bind thrice, for 3 replicates
D[] <- lapply(D, factor); D # convert to factors
y <- prop.edible <- c(0.95, 0.94, 0.95, 0.91, 0.90, 0.89,
                      0.85, 0.84, 0.85, 0.81, 0.80, 0.79,
                      0.97, 0.94, 0.95, 0.90, 0.92, 0.89,
                      0.84, 0.84, 0.87, 0.80, 0.80, 0.78,
                      0.95, 0.93, 0.99, 0.92, 0.90, 0.88,
                      0.85, 0.85, 0.85, 0.81, 0.81, 0.80); y; class(y)
popcorn <- data.frame( D, y ); popcorn
set.seed(2591) # ensures everyone has same randomization
D <- D[order(sample(1:36)), ]; D # randomize treatments

# 3.4.d pp 108-112
D <- expand.grid( brand = c(1, 2, 3), # brand 3 levels
                  power = c(1, 2), # power 2 levels
                  time = c(1, 2)) # time 2 levels
D <- rbind(D, D, D); D # row-bind thrice, for 3 replicates
D[] <- lapply(D, factor); D # convert to factors
y <- prop.edible <- c(0.95, 0.94, 0.95, 0.51, 0.50, 0.49, 
                      0.85, 0.84, 0.85, 0.41, 0.40, 0.41, 
                      0.97, 0.94, 0.96, 0.50, 0.52, 0.50, 
                      0.84, 0.84, 0.85, 0.40, 0.40, 0.39, 
                      0.95, 0.93, 0.95, 0.52, 0.50, 0.49, 
                      0.85, 0.85, 0.85, 0.41, 0.41, 0.42); y; class(y)
popcorn <- data.frame( D, y ); popcorn
mpopcorn <- lm(prop.edible ~ brand*power*time, data=popcorn )
anova(mpopcorn)

# 3.4.e pp 108-112
D <- expand.grid( brand = c(1, 2, 3), # brand 3 levels
                  power = c(1, 2), # power 2 levels
                  time = c(1, 2)) # time 2 levels
D <- rbind(D, D, D); D # row-bind thrice, for 3 replicates
D[] <- lapply(D, factor); D # convert to factors
y <- prop.edible <- c(0.95, 0.94, 0.95, 0.51, 0.50, 0.49, 
                      0.85, 0.84, 0.85, 0.41, 0.40, 0.41, 
                      0.97, 0.94, 0.96, 0.50, 0.52, 0.50, 
                      0.84, 0.84, 0.85, 0.40, 0.40, 0.39, 
                      0.95, 0.93, 0.95, 0.52, 0.50, 0.49, 
                      0.85, 0.85, 0.85, 0.41, 0.41, 0.42); y; class(y)
popcorn <- data.frame( D, y ); popcorn
mpopcorn <- lm(prop.edible ~ brand*power*time, data=popcorn )
anova(mpopcorn)

popcorn.t1 <- subset(popcorn, time == 1); popcorn.t1
power <- popcorn.t1$power
par (mfrow = c(1,2) )
interaction.plot(
  popcorn.t1$brand, 
  power, 
  popcorn.t1$y, type = "l",
  legend=TRUE, 
  ylim = c(.40,.99), 
  main = "Time = short", 
  xlab = "Brand", 
  ylab = "Proportion Popcorn Edible")
popcorn.t2 <- subset(popcorn, time == 2)
interaction.plot(
  popcorn.t2$brand, 
  power, 
  popcorn.t2$y, type = "l",
  legend=TRUE, 
  ylim = c(.40,.99), 
  main = "Time = long", 
  xlab = "Brand", 
  ylab = "Proportion Popcorn Edible")
par (mfrow = c(1,1) )

# 3.6.a pp 108-112
library(daewr) # one-replication per cell with numeric factors
D <- expand.grid( CPPU = c(0, 0.5, 1.0, 10.0), # design
                  dip.time = c(30, 60, 90) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
y <- growth <- c(92.5, 92.9, 91.3,
                 97.8, 94.9, 101.3,
                 97.0, 98.5, 101.6,
                 103.4, 102.9, 98.6); y
plant <- data.frame( D, growth )
m.plant <- lm(growth ~ CPPU*dip.time, data=plant )
anova(m.plant)

CPPUc <- as.ordered(plant$CPPU) # order CPPU factor levels
CPPULin <- contr.poly(CPPUc)[CPPUc,".L"] # linear part CPPU
CPPUQuad <- contr.poly(CPPUc)[CPPUc,".Q"] # quadratic part CPPU
CPPUCubic <- contr.poly(CPPUc)[CPPUc,".C"] # cubic part CPPU

dip.timec <- as.ordered(plant$dip.time) # order dip.time factor levels
dip.timeLin <- contr.poly(dip.timec)[dip.timec,".L"] # linear part dip.time
dip.timeQuad <- contr.poly(dip.timec)[dip.timec,".Q"] # quadratic part dip.time

m.plant.ll <-lm(growth ~ CPPUc + dip.timec 
                + CPPULin:dip.timeLin, data=plant)
anova(m.plant.ll) # linear x linear 
m.plant.ql <-lm(growth ~ CPPUc + dip.timec 
                + CPPUQuad:dip.timeLin, data=plant)
anova(m.plant.ql) # quadratic x linear 
m.plant.cl <-lm(growth ~ CPPUc + dip.timec 
                + CPPUCubic:dip.timeLin, data=plant)
anova(m.plant.cl) # cubic x linear 
m.plant.lq <-lm(growth ~ CPPUc + dip.timec 
                + CPPULin:dip.timeQuad, data=plant)
anova(m.plant.lq) # linear x quadratic 
m.plant.qq <-lm(growth ~ CPPUc + dip.timec 
                + CPPUQuad:dip.timeQuad, data=plant)
anova(m.plant.qq) # quadratic x quadratic 
m.plant.cq <-lm(growth ~ CPPUc + dip.timec 
                + CPPUCubic:dip.timeQuad, data=plant)
anova(m.plant.cq) # cubic x quadratic 

# 3.7 pp 108-112
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

# 3.8.a pp 108-112
library(daewr) # 2k design, 1 replicate
Des <- expand.grid( A = c(-1, 1), # 2^4 design
                    B = c(-1, 1), # 2 x 2 x 2 x 2 = 16 
                    C = c(-1, 1),
                    D = c(-1, 1)) # 16 x 1 replicates = 16

y <- c(45, 41, 78, 67, 50, 39, 95, 66, 
       47, 43, 95, 69, 40, 51, 87, 72); y
DesY <- data.frame(cbind(Des,y)); DesY 
m2k1r <- lm( y ~ A*B*C*D, data = DesY) # model m2k1r
summary(m2k1r)

par (mfrow = c(1,2) ) # plots to identify significant effects
fullnormal(coef(m2k1r)[-1],alpha=.025) # normal prob plot
LGB( coef(m2k1r)[-1], rpt = TRUE) # half-normal plot, with description
par (mfrow = c(1,1) )

# 3.8.b pp 108-112
library(daewr) # 2k design, 1 replicate
Des <- expand.grid( A = c(-1, 1), # 2^4 design
                    B = c(-1, 1), # 2 x 2 x 2 x 2 = 16 
                    C = c(-1, 1),
                    D = c(-1, 1)) # 16 x 1 replicates = 16

y <- c(45, 41, 78, 67, 50, 39, 95, 66, 
       47, 43, 95, 69, 40, 51, 87, 72); y
DesY <- data.frame(cbind(Des,y)); DesY 
m2k1r <- lm( y ~ A*B*C*D, data = DesY) # model m2k1r
summary(m2k1r)

par (mfrow = c(1,2) ) # plots to identify significant effects
fullnormal(coef(m2k1r)[-1],alpha=.025) # normal prob plot
LGB( coef(m2k1r)[-1], rpt = TRUE) # half-normal plot, with description
par (mfrow = c(1,1) )

# 3.8.c pp 108-112
library(BsMD) # Lenth and Bayes plots
par( mfrow = c(2,1) )
LenthPlot(m2k1r, main = "Lenth Plot of Effects") # lenth plot
X <- model.matrix(m2k1r)[ , 2:16] # 
y <- twok1rdata$y
twok1r.BsProb <- BsProb( 
  X = X, y = y, 
  blk = 0, 
  mFac = 15, # 15 factors
  mInt = 1, 
  p = 0.2, 
  g = 2.49, 
  ng = 1, 
  nMod = 10)
plot( twok1r.BsProb, main = "Bayes Plot of Effects" )
par (mfrow = c(1,1) )

