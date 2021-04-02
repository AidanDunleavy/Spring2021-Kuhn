# Lecture Notes 9
# Chapter 9.1-9.6 Crossover and Repeated Measures Designs

library(EMSaov)
library(AlgDesign)
library(BsMD)
library(car)
library(partitions)
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

# 9.1 two treatment AB, BA crossover design withOUT carryover
mouse <- factor(c(
  1, 4, 8, 6,
  1, 4, 8, 6,
  2, 3, 7, 5, 9,
  2, 3, 7, 5, 9
)); mouse
temp <- factor(c(
  "A", "A", "A", "A",
  "B", "B", "B", "B",
  "B", "B", "B", "B", "B",
  "A", "A", "A", "A", "A"
)); temp
period <- factor(c(
  1, 1, 1, 1,
  2, 2, 2, 2,
  1, 1, 1, 1, 1,
  2, 2, 2, 2, 2
)); period
roc <- c(
  5.9, 5.1, 5.3, 5.1,
  7.2, 8.0, 8.1, 8.3,
  7.1, 7.6, 6.9, 7.3, 7.2,
  6.1, 5.5, 5.6, 5.9, 5.8
); roc
data.cross <- data.frame(mouse,temp,period,roc); data.cross
m <- lm( roc ~ mouse + period + temp, data = data.cross,
         contrasts = list(period = contr.sum, mouse = contr.sum, temp = contr.sum))
summary.aov(m, type = "III" )
library(lsmeans)
lsmeans(m, pairwise ~ temp)
# 9.1 residual plots
par( mfrow = c(1,2) ) # plots slightly better with transformed data
plot( m, which = 1 ) # should be constant for different fitted values
plot( m, which = 2 ) # normal prob plot for normality should be linear
par(mfrow = c(1,1))
# 9.1 sample size estimation
alpha <- 0.05
sigma2 <- 0.235 # estimate of variance of residual error
delta <- 5 # required difference in means for different temperatures
n <- seq( 2, 10, by = 1)
stderrdiff <- sqrt( 2 * sigma2 / n)
df <- n - 2
t1 <- qt( 1 - alpha / 2, df )
gamma <- delta / stderrdiff
power <- 1 - pt(t1, df, gamma)
data.frame( alpha = alpha, n = n, delta = delta, power = power)

# 9.2 two treatment AB, BA crossover design with carryover
#   assuming carryover effects equal one another
group <- factor(c(
  1, 1, 1, 1,
  1, 1, 1, 1,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2
)); period
mouse <- factor(c(
  1, 4, 8, 6,
  1, 4, 8, 6,
  2, 3, 7, 5, 9,
  2, 3, 7, 5, 9
)); mouse
temp <- factor(c(
  "A", "A", "A", "A",
  "B", "B", "B", "B",
  "B", "B", "B", "B", "B",
  "A", "A", "A", "A", "A"
)); temp
period <- factor(c(
  1, 1, 1, 1,
  2, 2, 2, 2,
  1, 1, 1, 1, 1,
  2, 2, 2, 2, 2
)); period
roc <- c(
  5.9, 5.1, 5.3, 5.1,
  7.2, 8.0, 8.1, 8.3,
  7.1, 7.6, 6.9, 7.3, 7.2,
  6.1, 5.5, 5.6, 5.9, 5.8
); roc
data.cross <- data.frame(group,mouse,temp,period,roc); data.cross
library(lme4)
c1 <- c( .5, -.5)
m <- lmer( roc ~ 1 + group + (1|mouse:group) + period + temp, 
  contrasts = list(group = c1, period = c1, temp = c1),
  data = data.cross)
summary(m)

# 9.3 two treatment AB, BA crossover design which estimates carryover
group <- factor(c(
  1, 1, 1, 1,
  1, 1, 1, 1,
  1, 1, 1, 1,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2,
  2, 2, 2, 2, 2
)); period
mouse <- factor(c(
  1, 4, 8, 6,
  1, 4, 8, 6,
  1, 4, 8, 6,
  2, 3, 7, 5, 9,
  2, 3, 7, 5, 9,
  2, 3, 7, 5, 9
)); mouse
temp <- factor(c(
  "A", "A", "A", "A",
  "B", "B", "B", "B",
  "B", "B", "B", "B",
  "B", "B", "B", "B", "B",
  "A", "A", "A", "A", "A",
  "A", "A", "A", "A", "A"
)); temp
period <- factor(c(
  1, 1, 1, 1,
  2, 2, 2, 2,
  3, 3, 3, 3,
  1, 1, 1, 1, 1,
  2, 2, 2, 2, 2,
  3, 3, 3, 3, 3
)); period
carry <- factor(c(
  "none", "none", "none", "none",
  "A", "A", "A", "A",
  "B", "B", "B", "B",
  "none", "none", "none", "none", "none",
  "B", "B", "B", "B", "B",
  "A", "A", "A", "A", "A"
)); period
roc <- c(
  5.9, 5.1, 5.3, 5.1,
  7.2, 8.0, 8.1, 8.3,
  7.5, 8.1, 8.1, 8.0,
  7.1, 7.6, 6.9, 7.3, 7.2,
  6.1, 5.5, 5.6, 5.9, 5.8,
  5.7, 5.8, 5.6, 5.8, 5.5
); roc
data.cross <- data.frame(group,mouse,period,carry,temp,roc); data.cross
m <- lm( roc ~ mouse + period + temp + carry, 
  data = data.cross, 
  contrasts = list(mouse = contr.sum, period = contr.sum, temp = contr.sum, carry = contr.sum))
summary.aov(m, type = "III", singular.ok = TRUE)
library(lsmeans)
lsmeans(m, ~ temp)

# 9.4 multiple treatment crossover designs
library(crossdes)
wdes3 <- williams(3); wdes3
rownames(wdes3) <- paste("seqGroup", 1:6, sep = "")
colnames(wdes3) <- paste("Period", 1:3, sep = "")
wdes3
square <- factor(c(
  "I", "I", "I", "I", "I", "I", 
  "I", "I", "I", "I", "I", "I",
  "I", "I", "I", "I", "I", "I",
  "II", "II", "II", "II", "II", "II",
  "II", "II", "II", "II", "II", "II",
  "II", "II", "II", "II", "II", "II"
)); square
group <- factor(c(
  1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3,
  1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3
)); group
mouse <- factor(c(
  1, 2, 1, 2, 1, 2,
  3, 4, 3, 4, 3, 4,
  5, 6, 5, 6, 5, 6,
  7, 8, 7, 8, 7, 8,
  9, 10, 9, 10, 9, 10,
  11, 12, 11, 12, 11, 12
)); mouse
period <- factor(c(
  1, 1, 2, 2, 3, 3,
  1, 1, 2, 2, 3, 3,
  1, 1, 2, 2, 3, 3,
  1, 1, 2, 2, 3, 3,
  1, 1, 2, 2, 3, 3,
  1, 1, 2, 2, 3, 3
)); period
temp <- factor(c(
  "A", "A", "B", "B", "C", "C", 
  "B", "B", "C", "C", "A", "A",
  "C", "C", "A", "A", "B", "B",
  "C", "C", "B", "B", "A", "A",
  "A", "A", "C", "C", "B", "B",
  "B", "B", "A", "A", "C", "C"
)); temp
carry <- factor(c(
  "0", "0", "A", "A", "B", "B", 
  "0", "0", "B", "B", "C", "C",
  "0", "0", "C", "C", "A", "A",
  "0", "0", "C", "C", "B", "B",
  "0", "0", "A", "A", "C", "C",
  "0", "0", "B", "B", "A", "A"
)); carry
roc <- c(
  5.6, 5.3, 6.7, 6.1, 7.3, 7.3,
  6.5, 6.4, 7.2, 7.5, 5.8, 5.5,
  7.7, 7.7, 5.3, 5.8, 6.5, 6.5,
  7.0, 7.3, 6.6, 6.7, 5.6, 5.9,
  5.8, 5.2, 7.2, 7.3, 6.7, 6.0,
  6.9, 6.3, 5.8, 5.4, 7.0, 7.6
); roc
data.cross <- data.frame(square, group, 
  mouse, period, temp, carry, roc); data.cross
m <- lm(roc ~ mouse + period + temp  + carry, 
  data = data.cross, 
  contrasts = list(mouse = contr.sum, 
    period = contr.sum, temp = contr.sum, carry = contr.sum))
library(car)
summary.aov(m, type = "III", singular.ok = TRUE)
with(data.cross, tapply(roc, temp, mean))
sqrt(with(data.cross, tapply(roc, temp, var)))
with(data.cross, tapply(roc, carry, mean))
sqrt(with(data.cross, tapply(roc, carry, var)))
m <- mean(data.cross$roc); m # grand mean roc

# 9.5 repeated measures design, ignoring correlation between time periods
temp <- factor(c(
  "0", "0", "0", "0", 
  "0", "0", "0", "0", 
  "0", "0", "0", "0",
  "10", "10", "10", "10", 
  "10", "10", "10", "10", 
  "10", "10", "10", "10",
  "20", "20", "20", "20", 
  "20", "20", "20", "20", 
  "20", "20", "20", "20"
)); temp
mouse <- factor(c(
  1, 1, 1, 1, 
  2, 2, 2, 2, 
  3, 3, 3, 3,
  1, 1, 1, 1, 
  2, 2, 2, 2, 
  3, 3, 3, 3,
  1, 1, 1, 1, 
  2, 2, 2, 2, 
  3, 3, 3, 3
)); mouse
day <- factor(c(
  1, 2, 3, 4, 
  1, 2, 3, 4, 
  1, 2, 3, 4,
  1, 2, 3, 4, 
  1, 2, 3, 4, 
  1, 2, 3, 4,
  1, 2, 3, 4, 
  1, 2, 3, 4, 
  1, 2, 3, 4
)); day
roc <- c(
  7.6, 6.3, 5.7, 4.1, 
  7.3, 6.3, 5.5, 4.4, 
  7.2, 6.5, 5.8, 4.5,
  6.7, 5.7, 4.3, 3.8, 
  6.5, 5.5, 4.0, 3.3, 
  6.6, 5.7, 4.6, 3.9,
  5.8, 4.2, 3.2, 2.3, 
  5.7, 4.0, 3.9, 2.3, 
  5.8, 4.4, 3.0, 2.6
); roc
d.repeat <- data.frame(temp, mouse, day, roc); d.repeat
temperature <- d.repeat$temp
interaction.plot(d.repeat$day, 
      temperature, 
      d.repeat$roc, 
      type="b", pch=c(18,24,22),
      leg.bty="o",
      main="Trends in Mean ROCs over Time",
      xlab="Day",
      ylab="Average ROC")
library(GAD)
temp <- as.fixed(d.repeat$temp)
day <- as.fixed(d.repeat$day)
mouse <- as.random(d.repeat$mouse)
ROC <- d.repeat$roc
m.repeat <- lm(ROC ~ temp + mouse%in%temp + day + temp*day, data = d.repeat)
gad(m.repeat) # analysis using gad
library(lme4)
m.repeat <- lmer(ROC ~ 1 + temp*day + (1|mouse:temp), data = d.repeat)
anova(m.repeat) # analysis using lmer
# rearrange data for mauchly sphericity test
roc1 <- d.repeat$roc[c(1,5,9,13,17,21,25,29,33)]
roc2 <- d.repeat$roc[c(2,6,10,14,18,22,26,30,34)]
roc3 <- d.repeat$roc[c(3,7,11,15,19,23,27,31,35)]
roc4 <- d.repeat$roc[c(4,8,12,16,20,24,28,32,36)]
temp <- as.factor(d.repeat$temp[c(1:3, 13:15, 25:27)])
d.mult <- data.frame(temp, roc1, roc2, roc3, roc4); d.mult
# multivariate ROCs-for-4-days regressed on temperature analysis
m.mult <- lm(cbind(roc1, roc2, roc3, roc4) ~ temp, data = d.mult) 
# various multivariate tests to check for correlation between days
#  if correlation significant, then cannot perform repeated measures
library(car)
time <- factor(c(1:4)); time
idata <- data.frame(time); idata
m.i <- Anova(m.mult, idata = idata, idesign = ~ time)
summary(m.i, multivariate=FALSE) 

# 9.6 repeated measures design, accommodating correlation between time periods
#  summarizing data by regression, then using predicted values of regression
temp <- factor(c(
  "0", "0", "0", "0", 
  "0", "0", "0", "0", 
  "0", "0", "0", "0",
  "10", "10", "10", "10", 
  "10", "10", "10", "10", 
  "10", "10", "10", "10",
  "20", "20", "20", "20", 
  "20", "20", "20", "20", 
  "20", "20", "20", "20"
)); temp
mouse <- factor(c(
  1, 1, 1, 1, 
  2, 2, 2, 2, 
  3, 3, 3, 3,
  1, 1, 1, 1, 
  2, 2, 2, 2, 
  3, 3, 3, 3,
  1, 1, 1, 1, 
  2, 2, 2, 2, 
  3, 3, 3, 3
)); mouse
day <- factor(c(
  1, 2, 3, 4, 
  1, 2, 3, 4, 
  1, 2, 3, 4,
  1, 2, 3, 4, 
  1, 2, 3, 4, 
  1, 2, 3, 4,
  1, 2, 3, 4, 
  1, 2, 3, 4, 
  1, 2, 3, 4
)); day
roc <- c(
  7.6, 6.3, 5.7, 4.1, 
  7.3, 6.3, 5.5, 4.4, 
  7.2, 6.5, 5.8, 4.5,
  6.7, 5.7, 4.3, 3.8, 
  6.5, 5.5, 4.0, 3.3, 
  6.6, 5.7, 4.6, 3.9,
  5.8, 4.2, 3.2, 2.3, 
  5.7, 4.0, 3.9, 2.3, 
  5.8, 4.4, 3.0, 2.6
); roc
day.n <- as.numeric(day); day.n
r1 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[1:4,]))["as.numeric(day)"]
r2 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[5:8,]))["as.numeric(day)"]
r3 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[9:12,]))["as.numeric(day)"]
r4 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[13:16,]))["as.numeric(day)"]
r5 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[17:20,]))["as.numeric(day)"]
r6 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[21:24,]))["as.numeric(day)"]
r7 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[25:28,]))["as.numeric(day)"]
r8 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[29:32,]))["as.numeric(day)"]
r9 <- coef(lm(log(roc) ~ as.numeric(day), data = d.repeat[33:36,]))["as.numeric(day)"]
roc.est <- c(r1,r2,r3,r4,r5,r6,r7,r8,r9); roc.est
roc1 <- d.repeat$roc[c(1,5,9,13,17,21,25,29,33)]
roc2 <- d.repeat$roc[c(2,6,10,14,18,22,26,30,34)]
roc3 <- d.repeat$roc[c(3,7,11,15,19,23,27,31,35)]
roc4 <- d.repeat$roc[c(4,8,12,16,20,24,28,32,36)]
temp <- as.factor(d.repeat$temp[c(1:3, 13:15, 25:27)])
mouse <- 1:9
d.mult <- data.frame(mouse, temp, roc1, roc2, roc3, roc4, roc.est); d.mult
m.est <- aov(roc.est ~ temp, data = d.mult); summary(m.est)

# linear regression of growth curve plot of ROCs of mouse 7 by day 
data(d.mult)
dm7 <- d.mult[d.mult$mouse == 7, ]; dm7
lc <- log(c(dm7[1,3],dm7[1,4],dm7[1,5],dm7[1,6]))
day <- c(1,2,3,4)
fl <- lm(lc~day)
plot(lc~day, 
  ylab="log ROC", 
  main="Log ROC by Day for Mouse 7")
abline(coef(fl))

###################################################
# Homework 9 and Test 5 questions
###################################################

# 9.2.a,b,c pp 379-382
# Cross-over design with extra-period
square <- factor(c(
  "I", "I", "I", "I", "I", "I", "I", "I", 
  "I", "I", "I", "I", "I", "I", "I", "I",
  "I", "I", "I", "I", "I", "I", "I", "I",
  "II", "II", "II", "II", "II", "II", "II", "II",
  "II", "II", "II", "II", "II", "II", "II", "II",
  "II", "II", "II", "II", "II", "II", "II", "II"
)); square
group <- factor(c(
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3,
  1, 1, 1, 1, 1, 1, 1, 1,
  2, 2, 2, 2, 2, 2, 2, 2,
  3, 3, 3, 3, 3, 3, 3, 3
)); group
cow <- factor(c(
  1, 1, 1, 1, 2, 2, 2, 2,
  3, 3, 3, 3, 4, 4, 4, 4,
  5, 5, 5, 5, 6, 6, 6, 6,
  7, 7, 7, 7, 8, 8, 8, 8,
  9, 9, 9, 9, 10, 10, 10, 10,
  11, 11, 11, 11, 12, 12, 12, 12
)); cow
period <- factor(c(
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4,
  1, 2, 3, 4, 1, 2, 3, 4
)); period
diet <- factor(c(
  "A", "B", "C", "C", "A", "B", "C", "C",
  "B", "C", "A", "A", "B", "C", "A", "A",
  "C", "A", "B", "B", "C", "A", "B", "B",
  "A", "C", "B", "B", "A", "C", "B", "B",
  "B", "A", "C", "C", "B", "A", "C", "C",
  "C", "B", "A", "A", "C", "B", "A", "A"
)); diet
carry <- factor(c(
  "0", "A", "B", "C", "0", "A", "B", "C", 
  "0", "B", "C", "A", "0", "B", "C", "A",
  "0", "C", "A", "B", "0", "C", "A", "B",
  "0", "A", "C", "B", "0", "A", "C", "B",
  "0", "B", "A", "C", "0", "B", "A", "C",
  "0", "C", "B", "A", "A", "0", "C", "B"
)); carry
fcm <- c(
  38.66, 37.43, 34.39, 31.30, 25.72, 26.13, 23.35, 18.69,
  48.85, 46.88, 41.99, 39.61, 30.80, 29.29, 26.41, 23.16,
  34.64, 32.27, 28.50, 27.13, 25.35, 26.00, 23.86, 19.92,
  35.19, 33.50, 28.41, 25.12, 21.80, 23.91, 21.69, 17.55,
  32.90, 33.12, 27.52, 25.10, 21.37, 21.97, 19.38, 16.57,
  30.40, 29.50, 26.70, 23.09, 22.84, 20.97, 18.59, 16.10
); fcm
data <- data.frame(square, group, cow, period, diet, carry, fcm); data
m <- lm(fcm ~ cow + period + diet + carry, 
  data = data, 
  contrasts = list(cow = contr.sum, period = contr.sum, 
  diet = contr.sum, carry = contr.sum))
library(car)
summary.aov(m, type = "I", singular.ok = TRUE)
summary.aov(m, type = "III", singular.ok = TRUE)

with(data, tapply(fcm, diet, mean))
sqrt(with(data, tapply(fcm, diet, var)))
with(data, tapply(fcm, carry, mean))
sqrt(with(data, tapply(fcm, carry, var)))
m <- mean(data$fcm); m # grand mean roc

# 9.4.a,b pp 379-382
# Repeated measures design
library(daewr)
residue # concentration of herbicide (residue) data
m.mult <- lm(cbind(X1, X2, X3, X4, X5) ~ temp * moisture * soil, data = residue) 
# various multivariate tests to check for correlation between days
#  if correlation significant, then cannot perform repeated measures
library(car)
time <- factor(c(1:5)); time # 5 time periods, matches 5 days
idata <- data.frame(time); idata
m.i <- Anova(m.mult, type = "III", idata = idata, idesign = ~ time)
summary(m.i, multivariate=FALSE) 

# 9.4.c pp 379-382
# Repeated measures design
# not worrying about correlation between days
temp <- c(
  10, 10, 10, 10, 10,  
  10, 10, 10, 10, 10,  
  10, 10, 10, 10, 10,  
  10, 10, 10, 10, 10,  
  
  30, 30, 30, 30, 30,  
  30, 30, 30, 30, 30,  
  30, 30, 30, 30, 30,  
  30, 30, 30, 30, 30,  
  
  10, 10, 10, 10, 10,  
  10, 10, 10, 10, 10,  
  10, 10, 10, 10, 10,  
  10, 10, 10, 10, 10,  
  
  30, 30, 30, 30, 30,  
  30, 30, 30, 30, 30,  
  30, 30, 30, 30, 30,  
  30, 30, 30, 30, 30  
); temp
moist <- c(
  "L", "L", "L", "L", "L",  
  "L", "L", "L", "L", "L",  
  "H", "H", "H", "H", "H",  
  "H", "H", "H", "H", "H",  
  
  "L", "L", "L", "L", "L",  
  "L", "L", "L", "L", "L",  
  "H", "H", "H", "H", "H",  
  "H", "H", "H", "H", "H",  
  
  "L", "L", "L", "L", "L",  
  "L", "L", "L", "L", "L",  
  "H", "H", "H", "H", "H",  
  "H", "H", "H", "H", "H",  
  
  "L", "L", "L", "L", "L",  
  "L", "L", "L", "L", "L",  
  "H", "H", "H", "H", "H",  
  "H", "H", "H", "H", "H"  
); moist
soil <- c(
  "C", "C", "C", "C", "C",  
  "P", "P", "P", "P", "P",  
  "C", "C", "C", "C", "C",  
  "P", "P", "P", "P", "P",  
  
  "C", "C", "C", "C", "C",  
  "P", "P", "P", "P", "P",  
  "C", "C", "C", "C", "C",  
  "P", "P", "P", "P", "P",  
  
  "C", "C", "C", "C", "C",  
  "P", "P", "P", "P", "P",  
  "C", "C", "C", "C", "C",  
  "P", "P", "P", "P", "P",  
  
  "C", "C", "C", "C", "C",  
  "P", "P", "P", "P", "P",  
  "C", "C", "C", "C", "C",  
  "P", "P", "P", "P", "P"  
); soil
day <- c(
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60,  
  0, 7, 14, 30, 60  
); day
pot <- c(
  1, 1, 1, 1, 1,  
  1, 1, 1, 1, 1,  
  1, 1, 1, 1, 1,  
  1, 1, 1, 1, 1,  
  
  1, 1, 1, 1, 1,  
  1, 1, 1, 1, 1,  
  1, 1, 1, 1, 1,  
  1, 1, 1, 1, 1,  
  
  2, 2, 2, 2, 2,  
  2, 2, 2, 2, 2,  
  2, 2, 2, 2, 2,  
  2, 2, 2, 2, 2,  
  
  2, 2, 2, 2, 2,  
  2, 2, 2, 2, 2,  
  2, 2, 2, 2, 2,  
  2, 2, 2, 2, 2  
); pot
conc <- c(
  0.77, 0.77, 0.76, 0.74, 0.72,  
  0.78, 0.76, 0.75, 0.72, 0.66,  
  0.76, 0.74, 0.71, 0.66, 0.57, 
  0.78, 0.76, 0.74, 0.70, 0.63,  
  0.77, 0.74, 0.71, 0.65, 0.54,  
  0.79, 0.73, 0.68, 0.58, 0.42,  
  0.78, 0.73, 0.69, 0.60, 0.46,  
  0.78, 0.71, 0.65, 0.53, 0.36,  
  0.77, 0.76, 0.75, 0.72, 0.66,  
  0.77, 0.76, 0.74, 0.71, 0.65,  
  0.78, 0.77, 0.75, 0.72, 0.67,  
  0.77, 0.75, 0.73, 0.68, 0.60,  
  0.79, 0.75, 0.72, 0.65, 0.54,  
  0.78, 0.74, 0.69, 0.60, 0.45,  
  0.79, 0.72, 0.65, 0.53, 0.35,  
  0.78, 0.70, 0.63, 0.49, 0.31
); conc
library(GAD)
data <- data.frame(temp, moist, soil, day, pot, conc); data
T <- as.fixed(data$temp)
M <- as.fixed(data$moist)
S <- as.fixed(data$soil)
D <- as.fixed(data$day)
P <- as.random(data$pot)
C <- data$conc
mdl <- lm(C ~ T + M + S + T*M + T*S + M*S + T*M*S + P%in%(T:M:S) 
          + D + T*D + M*D + S*D + T*M*D + T*S*D + M*S*D + T*M*S*D, data = data); 
library(car)
gad(mdl) # repeated measures analysis using gad, not worrying about correlation



