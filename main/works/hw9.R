data = read.table("https://raw.githubusercontent.com/AidanDunleavy/Spring2021-Kuhn/main/main/lib/Exercise.txt", header = TRUE)
data
datum <- data.frame(pulseTime = rep(1:3, each = 18), Exercise = rep(data$Exercise, 3), Diet = rep(data$Diet, 3), Pulse = c(data$Pulse1, data$Pulse2, data$Pulse3), comboDietEx = rep(1:6,each=3))
datum
data1 <- datum[which(datum$Diet == 1),]
data1 # MEAT!!!!
data2 <- datum[which(datum$Diet == 2),]
data2 # VEGGIES

# exercise 1 = aerobic stair climbing / 2 is racquetball / 3 weight training

Exercise = data1$Exercise
interaction.plot(data1$pulseTime, 
                 Exercise, 
                 data1$Pulse, 
                 type="b", pch=c(18,24,22),
                 leg.bty="o",
                 main="Trends in Pulse over Times taken FOR MEAT!!!",
                 xlab="Pulse Time",
                 ylab="Pulse")


Exercise = data2$Exercise
interaction.plot(data2$pulseTime, 
                 Exercise, 
                 data2$Pulse, 
                 type="b", pch=c(18,24,22),
                 leg.bty="o",
                 main="Trends in Pulse over Times taken FOR VEGGIES!!!",
                 xlab="Pulse Time",
                 ylab="Pulse")



interaction.plot(datum$pulseTime, 
                 datum$comboDietEx, 
                 datum$Pulse, 
                 type="b", pch=c(18,24,22, 23, 21, 19),
                 leg.bty="o",
                 main="Trends in Pulse over Times taken",
                 xlab="Pulse Time",
                 ylab="Pulse")
############################################

data = read.table("https://raw.githubusercontent.com/AidanDunleavy/Spring2021-Kuhn/main/main/lib/Exercise.txt", header = TRUE)
data2 <- cbind(data, subject = rep(1:3, 6))
data2

library(FrF2)
SPFF2 <- FrF2(18, 7, # 18 runs, 7 factors 
              WPs = 6, nfac.WP = 2, # 6 whole plots, 2 whole-plot factors
              factor.names = c("A","B","C","D","E","F","G"), randomize = FALSE)
y <- rnorm(18, 0, 1)
aliases(lm( y ~ (.)^2, data = SPFF2))
print(SPFF2)

library(daewr)
par(mfrow = c(1,2))
fullnormal(Wpeffect, names(Wpeffect), alpha = .10)
fullnormal(Speffect, names(Speffect), alpha = .05)
par(mfrow = c(1,1))


#########################################
data = read.table("https://raw.githubusercontent.com/AidanDunleavy/Spring2021-Kuhn/main/main/lib/Exercise.txt", header = TRUE)

m.mult <- lm(cbind(Pulse1, Pulse2, Pulse3) ~ Exercise + Diet, data = data) 
# various multivariate tests to check for correlation between days
#  if correlation significant, then cannot perform repeated measures
library(car)
time <- factor(c(1:3)); time
idata <- data.frame(time); idata
m.i <- Anova(m.mult, idata = idata, idesign = ~ time)
summary(m.i, multivariate=FALSE) 
#################################################################################

data = read.table("https://raw.githubusercontent.com/AidanDunleavy/Spring2021-Kuhn/main/main/lib/Exercise.txt", header = TRUE)
data
datum <- data.frame(pulseTime = rep(1:3, each = 18), Exercise = rep(data$Exercise, 3), Diet = rep(data$Diet, 3), Pulse = c(data$Pulse1, data$Pulse2, data$Pulse3), comboDietEx = rep(1:6,each=3))
datum

as.numeric(datum$pulseTime)[1:3]
log(datum$Pulse[1:3])
r1 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(1,19,37),]))["as.numeric(pulseTime)"]
r2 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(2,20,38),]))["as.numeric(pulseTime)"]
r3 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(3,21,39),]))["as.numeric(pulseTime)"]
r4 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(4,22,40),]))["as.numeric(pulseTime)"]
r5 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(5,23,41),]))["as.numeric(pulseTime)"]
r6 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(6,24,42),]))["as.numeric(pulseTime)"]
r7 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(7,25,43),]))["as.numeric(pulseTime)"]
r8 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(8,26,44),]))["as.numeric(pulseTime)"]
r9 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(9,27,45),]))["as.numeric(pulseTime)"]
r10 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(10,28,46),]))["as.numeric(pulseTime)"]
r11 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(11,29,47),]))["as.numeric(pulseTime)"]
r12 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(12,30,48),]))["as.numeric(pulseTime)"]
r13 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(13,31,49),]))["as.numeric(pulseTime)"]
r14 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(14,32,50),]))["as.numeric(pulseTime)"]
r15 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(15,33,51),]))["as.numeric(pulseTime)"]
r16 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(16,34,52),]))["as.numeric(pulseTime)"]
r17 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(17,35,53),]))["as.numeric(pulseTime)"]
r18 <- coef(lm(log(Pulse) ~ as.numeric(pulseTime), data = datum[c(18,36,54),]))["as.numeric(pulseTime)"]

roc.est <- c(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18); roc.est

datum <- cbind(datum, roc.est)
datum
m.est <- aov(roc.est ~ Exercise + Diet, data = datum); summary(m.est)











###############################################################################













data = read.table("https://raw.githubusercontent.com/AidanDunleavy/Spring2021-Kuhn/main/main/lib/c9_p1.txt", header = TRUE)
data

subject = c(rep(1:34,2),rep(1:30,2))
resp = c(data$G1Placebo,data$G1Test,na.omit(data$G2Placebo),na.omit(data$G2Test))
resp2 = c(data$G1Placebo,data$G1Test,data$G2Test,data$G2Placebo)
group = c(rep(1,68),rep(2,60))
period = c(rep(1:2, each = 34), rep(1:2,each=30))
treatment = c(rep(1,34), rep(2,64), rep(1,30))


data <- data.frame(subject, resp, group, period, treatment)
group1 <- data[which(data$group == 1),]
group2 <- data[which(data$group == 2),]
group2 <- c(group2[1:30,2], rep(NA,4), group2[31:60,2], rep(NA,4))
group2
c(rep(1,34), rep(2,34), rep(1,30), rep(NA,4), rep(2,30), rep(NA,4))

data <- data.frame(subject = subject[1:68], group1 = group1[1:68,2], group2 = group2, Treatment = rep(1:2, each = 34))
data

m.mult <- lm(cbind(group1, group2) ~ Treatment, data = data) 

library(car)
time <- factor(c(1:2)); time
idata <- data.frame(time); idata
m.i <- Anova(m.mult, idata = idata, idesign = ~ time)
m.i


