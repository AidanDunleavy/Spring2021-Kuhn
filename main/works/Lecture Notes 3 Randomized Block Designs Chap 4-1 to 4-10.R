# Lecture Notes 3
# Chapter 4 Randomized Block Design

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

# 4.1 Introduction
# 4.2 Creating an RCB in R

# random assignment for the block
library(agricolae)
treat <- c(0, 10, 20, 30) 
# seed 11 gives a single randomization for ROC data
outdesign <- design.rcbd(treat, 8, seed = 11); rcb 
rcb <- outdesign$book; rcb
# ROC of mice for different temps, block over mice
levels(rcb$block) <- c("mouse1", "mouse2", "mouse3", "mouse4",
                       "mouse5", "mouse6", "mouse7", "mouse8"); rcb

# 4.3 Model for RCB
# 4.4 An Example of an RCB

# input data for ROC for temp factor, with mouse block
D <- expand.grid( temp = c(0, 10, 20, 30), # design
                  mouse = c(1, 2, 3, 4,
                            5, 6, 7, 8) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
ROC <- c(
  10.3, 1.8, 1.2, 12.4, 
  9.1, 12.1, 6.5, 16.1, 
  6.1, 5.1, 1.2, 18.1, 
  7.2, 9.8, 8.1, 15.1,
  5.4, 4.2, 4.1, 17.2,
  7.0, 9.8, 8.1, 14.1,
  5.4, 4.3, 4.1, 17.2,
  2.1, 6.2, 2.1, 19.1); 
ROCdata <- data.frame(cbind(D,ROC)); ROCdata # data.frame

library(daewr) # RCB = factorial without interaction
m.roc1 <- aov( ROC ~ mouse + temp, data = ROCdata )
summary(m.roc1)

# investigate linear, quadratic, cubic contrast trends
contrasts(ROCdata$temp) <- contr.poly(4) 
m.roc2 <- aov( ROC ~ mouse + temp, data = ROCdata )
summary.aov(m.roc2 ,
    split = list( # splits out contrasts for temp factor
      temp = list("Linear" = 1,
                  "Quadratic" = 2,
                  "Cubic" = 3) ) )

# plot of trends, particularly quadratic trend
argue <- do.call("cbind", split(ROCdata$ROC, ROCdata$mouse)) # store list 
y <- apply(argue, 1, mean ) # determine ROC means for each mouse
x <- as.double( levels(ROCdata$temp) ) 
plot( x, y, xlab = "temperature", ylab = "average ROC" )
xx <- seq( 0, 30, 10 )
temp.lin <- lm( y ~ poly( x, 1) )
lines(xx, predict( temp.lin, data.frame( x = xx) ))

# 4.5 Determining the Number of Blocks

D <- expand.grid( temp = c(0, 10, 20, 30), # design
                  mouse = c(1, 2, 3, 4,
                            5, 6, 7, 8) ) # combinations of factors
D[] <- lapply(D, factor) # convert to factors
ROC <- c(
  10.3, 1.8, 1.2, 12.4, 
  9.1, 12.1, 6.5, 16.1, 
  6.1, 5.1, 1.2, 18.1, 
  7.2, 9.8, 8.1, 15.1,
  5.4, 4.2, 4.1, 17.2,
  7.0, 9.8, 8.1, 14.1,
  5.4, 4.3, 4.1, 17.2,
  2.1, 6.2, 2.1, 19.1); 
ROCdata <- data.frame(cbind(D,ROC)); ROCdata # data.frame
m.roc1 <- aov( ROC ~ mouse + temp, data = ROCdata )
summary(m.roc1)

model.tables(m.roc1, type = "means", se = T) # various means, for corrected sums of squares (css)

css <- (6.575 - 8.45625)^2 + (6.662 - 8.45625)^2 + (4.425 - 8.45625)^2 + (16.163 - 8.45625)^2; css
library(daewr)
bmin <- 2
bmax <- 8
alpha <- .05
sigma2 <- 7.29
css <- 82.40341
nu1 <- 4-1
blocks <- c(bmin:bmax)
nu2 <- (blocks - 1) * 4
nc <- (blocks * css) / sigma2
Power <- Fpower( alpha, nu1, nu2, nc )
data.frame(blocks, nu1, nu2, nc, Power)

# 4.6 Factorial Designs in Blocks

D <- expand.grid( 
  temp = c("0", "10", "20", "30"), # temperature factor B
  noise = c("low","high"), # noise factor A
  block = c(1, 2) # mouse
   ) # combinations of factors
D[] <- lapply(D, factor); D # convert to factors
ROC <- c(
  10.3, 1.8, 1.2, 12.4, 6.1, 5.1, 1.2, 18.1,
  9.1, 12.1, 6.5, 16.1, 7.2, 9.8, 8.1, 15.1
); 
ROCdata2 <- data.frame(cbind(D,ROC)); ROCdata2 # data.frame
library(daewr)
m.roc2 <- aov( ROC ~ block + temp * noise, data = ROCdata2)
summary(m.roc2) # categorical factors, different if numerical factors

with(ROCdata2, # 
   (interaction.plot(
     noise, 
     temp, 
     ROC, 
     type = "b", 
     pch = c(18,24,22,21), 
     leg.bty = "o",
     main = "ROC for each Temperature",
     xlab = "Noise",
     ylab = "ROC"))
)

# 4.7 Generalized Complete Block Design

D <- expand.grid( 
  repl = c(1, 2, 3), # replicates
  temp = c("0", "10", "20", "30"), # temperature factor 
  mouse = c(1, 2, 3, 4, 5, 6, 7, 8) # block
) # combinations of factors
D[] <- lapply(D, factor); D # convert to factors
ROC <- c(
  1.4, 6.1, 3.4, 9.4, 8.1, 5.4, 7.4, 9.1, 8.4, 10.4, 16.1, 17.4,
  1.8, 5.1, 4.2, 4.8, 5.1, 7.2, 8.8, 5.1, 7.2, 14.8, 16.1, 19.2,
  1.2, 1.2, 4.1, 1.2, 4.2, 4.1, 1.2, 1.2, 6.1, 11.2, 15.2, 18.1,
  3.4, 5.1, 6.2, 12.4,8.1, 7.2,10.4, 8.1,17.2, 13.4, 18.1, 14.2,
  3.1, 4.2, 3.0, 4.1, 7.2, 7.0, 9.1, 7.7, 9.0, 11.1, 17.2, 16.0,
  2.1, 2.8, 3.8, 2.1, 9.8, 9.8,12.1, 9.8, 9.8, 12.1, 19.8, 17.8,
  1.5, 1.1, 1.1, 6.5, 4.1, 7.1, 6.5, 8.1, 7.1, 10.5, 10.1, 12.1,
  2.1, 5.1, 3.1, 3.1, 5.1, 7.1, 9.1, 8.1, 7.1, 16.1, 15.1, 14.1
); 
ROCdata3 <- data.frame(cbind(D,ROC)); ROCdata3 # data.frame

# GCB, generalized complete block design
library(daewr)
data(ROCdata3); ROCdata3
m.roc3 <- aov(ROC ~ temp + Error(mouse/temp), data = ROCdata3)
summary(m.roc3)

# incorrect method, should use temp x mouse interaction instead
m.roc3a <- aov( ROC ~ temp*mouse, data = ROCdata3)
summary(m.roc3a)

# alternative (correct) method, average replicates 
# TukeyHSD comparison of means
cellmeans <- tapply( ROCdata3$ROC, list(ROCdata3$mouse, ROCdata3$temp), mean); cellmeans
dim(cellmeans) <- NULL
temp <- factor( rep(c(0,10,20,30), each = 8) )
mouse <- factor( rep(c(1,2,3,4,5,6,7,8), 4) )
m.roc4 <- aov( cellmeans ~ mouse + temp ); summary(m.roc4)
model.tables( m.roc4, type = "means" )$tables$temp
TukeyHSD( m.roc4, "temp" )

# block x temp interaction insignificant at alpha = 0.25
# so use additive model instead
m.roc5 <- aov( ROC ~ mouse + temp, data = ROCdata3)
summary(m.roc5)

# 4.8 Two Block Factors LSD

# random assignment of temp for LSD
library(agricolae)
temp <- c("0", "10", "20", "30")
outdesign <- design.lsd( temp, seed = 23)
lsd <- outdesign$book
levels(lsd$row) <- c("Mouse 1", "Mouse 2", "Mouse 3", "Mouse 4")
levels(lsd$col) <- c("Time 1", "Time 2", "Time 3", "Time 4")
lsd

# LSD analysis
ROC <- c(
  10.3, 1.8, 1.2, 12.4, 6.1, 5.1, 1.2, 18.1,
  9.1, 12.1, 6.5, 16.1, 7.2, 9.8, 8.1, 15.1
); 
ROCdata6 <- data.frame(cbind(lsd,ROC)); ROCdata6 # data.frame
colnames(ROCdata6) <- c("plots","mouse","time","temp","ROC"); ROCdata6
m.roc6 <- aov( ROC ~ mouse + time + temp, data = ROCdata6)
summary(m.roc6)

# averages of ROC means
model.tables( m.roc6, type = "means" )$tables$temp
# tukey pairwise comparison of ROC means
TukeyHSD( m.roc6, "temp")

# power of LSD
t <- 4 # number of treatments
b <- 4 # number of blocks (for both block factors)
nu1 <- t - 1 # df for 4 treatment levels
nu2 <- (t - 1) * (t - 2) # df for ANOVA residuals
alpha <- .05
sigma2lsd <- 11.18
css <- 22.36
nc <- (b * css) / sigma2lsd
Power <- Fpower( alpha, nu1, nu2, nc )
data.frame(b, nu1, nu2, nc, Power)

# 4.9 Review of Important Concepts

###################################################
# Test 2 questions
###################################################

# 4.2.b pp 135-140
library(agricolae)
flower <- c("carnation", "daisy", "rose", "tulip")
outdesign <- design.lsd( flower, seed = 23)
lsd <- outdesign$book
levels(lsd$row) <- c("tap", "tap.sugar", "tap.carb", "tap.7up")
levels(lsd$col) <- c("65", "68", "70", "72")
lsd

# 4.3 pp 135-140
D <- expand.grid(barley=c("manchuria","svansota","velvet","trebi","peatland"),year.place=c("1931-1","1932-1","1931-2","1932-2","1931-3","1932-3","1931-4","1932-4","1931-5","1932-5","1931-6","1932-6"))

D[] = lapply(D, factor)

growthlevels=c(81,105.4,119.7,109.7,98.3,80.7,82.3,80.4,87.2,84.2,146.6,142,150.7,191.5,145.7,100.4,115.5,112.2,147.7,108.1,82.3,77.3,78.4,131.3,89.6,103.1,105.1,116.5,139.9,129.6,119.8,121.4,124,140.8,124.8,98.9,61.9,96.2,125.5,75.7,98.9,89,69.1,89.3,104.1,66.4,49.9,96.7,61.9,80.3,86.9,77.1,78.9,101.8,96,67.7,66.7,67.4,91.8,94.1)

GLdata <- data.frame(cbind(D,growthlevels)); GLdata
m.gl <- aov(growthlevels ~ barley + year.place, data=GLdata) 
summary(m.gl) # RCB design

# 4.6.b pp 135-140
library(agricolae)
drug <- c("A", "B", "C", "D", "E")
outdesign <- design.lsd( drug, seed = 23)
lsd <- outdesign$book
levels(lsd$row) <- c("Patient 1", "Patient 2", "Patient 3", "Patient 4", "Patient 5")
levels(lsd$col) <- c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5")
lsd$drug <- c("B","E","A","C","D",
              "D","A","E","B","C",
              "E","B","C","D","A",
              "A","C","D","E","B",
              "C","D","B","A","E"); lsd
sleep.rating <- c(2.92, 2.43, 2.19, 2.71, 2.71, 
               2.86, 1.64, 3.02, 3.03, 3.03, 
               1.97, 2.50, 2.47, 2.65, 1.89, 
#               2.99, 3.39, 3.37, 3.33, 3.71, # 4.6.b..
               1.99, 2.39, 2.37, 2.33, 2.71,
               2.64, 2.31, 2.44, 1.89, 2.78
); 
sleep.data <- data.frame(cbind(lsd,sleep.rating)); sleep.data # data.frame
colnames(sleep.data) <- c("ID","patient","week","drug","sleep.rating"); sleep.data
m.sleep <- aov( sleep.rating ~ patient + week + drug, 
                data = sleep.data)
summary(m.sleep)

# 4.6.c pp 135-140
library(agricolae)
drug <- c("A", "B", "C", "D", "E")
outdesign <- design.lsd( drug, seed = 23)
lsd <- outdesign$book
levels(lsd$row) <- c("Patient 1", "Patient 2", "Patient 3", "Patient 4", "Patient 5")
levels(lsd$col) <- c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5")
lsd$drug <- c("B","E","A","C","D",
              "D","A","E","B","C",
              "E","B","C","D","A",
              "A","C","D","E","B",
              "C","D","B","A","E"); lsd
sleep.rating <- c(2.92, 2.43, 2.19, 2.71, 2.71, 
                  2.86, 1.64, 3.02, 3.03, 3.03, 
                  1.97, 2.50, 2.47, 2.65, 1.89, 
                  1.99, 2.39, 2.37, 2.33, 2.71,
                  2.64, 2.31, 2.44, 1.89, 2.78
); 
sleep.data <- data.frame(cbind(lsd,sleep.rating)); sleep.data # data.frame
colnames(sleep.data) <- c("ID","patient","week","drug","sleep.rating"); sleep.data
m.sleep <- aov( sleep.rating ~ patient + week + drug, 
                data = sleep.data)
summary(m.sleep)

#   c1: placebo vs ave of others contrast
library(gmodels)
model.tables( m.sleep, type = "means" )
fit.contrast( 
  m.sleep, 
  "drug", 
  rbind( "plcb vs ave(B,C,D,E)" = c(1,-1/4,-1/4,-1/4,-1/4) ), 
  conf=0.95 )
cm <- rbind( "plcb vs ave(B,C,D,E)" = c(1,-1/4,-1/4,-1/4,-1/4))
m.sleepc <- aov( 
    sleep.rating ~ patient + week + drug, 
    data = sleep.data,
    contrasts=list( drug=make.contrasts(cm) )
)
summary(m.sleepc, 
      split=list(drug = list( "plcb vs ave(B,C,D,E)" = 1 )))

#   c2: among the four drugs
library(agricolae)
compare <- SNK.test( m.sleep, "drug",  alpha = 0.05 )
print(compare)

# 4.6.d pp 135-140
library(agricolae)
drug <- c("A", "B", "C", "D", "E")
outdesign <- design.lsd( drug, seed = 23)
lsd <- outdesign$book
levels(lsd$row) <- c("Patient 1", "Patient 2", "Patient 3", "Patient 4", "Patient 5")
levels(lsd$col) <- c("Week 1", "Week 2", "Week 3", "Week 4", "Week 5")
lsd$drug <- c("B","E","A","C","D",
              "D","A","E","B","C",
              "E","B","C","D","A",
              "A","C","D","E","B",
              "C","D","B","A","E"); lsd
sleep.rating <- c(2.92, 2.43, 2.19, 2.71, 2.71, 
                  2.86, 1.64, 3.02, 3.03, 3.03, 
                  1.97, 2.50, 2.47, 2.65, 1.89, 
                  1.99, 2.39, 2.37, 2.33, 2.71,
                  2.64, 2.31, 2.44, 1.89, 2.78
); 
sleep.data <- data.frame(cbind(lsd,sleep.rating)); sleep.data # data.frame
colnames(sleep.data) <- c("ID","patient","week","drug","sleep.rating"); sleep.data
m.sleep <- aov( sleep.rating ~ patient + week + drug, 
                data = sleep.data)
summary(m.sleep)

par( mfrow = c(1,2) ) # plots slightly better with transformed data
plot( m.sleep, which = 1 ) # should be constant for different fitted values
plot( m.sleep, which = 2 ) # normal prob plot for normality should be linear
par(mfrow = c(1,1))

# 4.8.a pp 135-140
library(daewr)
t <- 4 # number of treatments
bmin <- 2 # try number of blocks 2 to 5
bmax <- 5
blocks <- c(bmin:bmax)
nu1 <- t - 1
nu2 <- (blocks - 1) * (t - 1)
alpha <- .05
sigma2rcb <- 0.5
css <- 2.0
nc <- (blocks * css) / sigma2rcb
Power <- Fpower( alpha, nu1, nu2, nc )
data.frame(blocks, nu1, nu2, nc, Power)

# 4.8.b pp 135-140
library(daewr)
t <- 4 # number of treatments
b <- 4 # number of blocks (for both block factors)
nu1 <- t - 1 # df for 4 treatment levels
nu2 <- (t - 1) * (t - 2) # df for ANOVA residuals
alpha <- .05
sigma2lsd <- 0.45
css <- 2.0
nc <- (b * css) / sigma2lsd
Power <- Fpower( alpha, nu1, nu2, nc )
data.frame(b, nu1, nu2, nc, Power)
