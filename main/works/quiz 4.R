# 8.2.a,b pp 345-350
# Split-Plot with CRD in Whole Plots (CRSP) Designs
# 1 whole plot with 3 levels replicated twice, 1 subplot with 3 levels
library(AlgDesign)
wp<-data.frame(A=factor(c(1,2,3)))
wp<-rbind(wp,wp) # 1 whole factor: 3 levels replicated twice: 3 x 2 = 6
sp<-expand.grid(B=factor(c(1,2,3))) # 1 subplot factor: 3 levels
splitP<-optBlock( ~ A*B, 
                  withinData = sp, 
                  blocksizes = rep(3,6), 
                  wholeBlockData=wp)
splitP$Blocks

# 8.2.c,d pp 345-350
# Split-Plot with CRD in Whole Plots (CRSP) Designs
# 2 whole plot with 2 levels replicated twice, 2 subplots with 3 levels
# 2^2 x 3^2 = 36 runs
A <- factor(rep(c(1,1,2,2), each = 9)); A # 36 runs for whole factor A
B <- factor(rep(c(1,2,1,2), each = 9)); B # whole factor B
wp <- factor(rep(c(1,2,3,4), each = 9)); wp
C <- factor(rep(c(1,2,3,1,2,3,1,2,3,1,2,3), each = 3)); C # subplot factor C
D <- factor(rep(c(1,2,3), times = 12)); D # subplot factor D
sp <- factor(rep(c(1,2,3,4,5,6,7,8,9), times = 3)); sp
splitP <- data.frame(cbind(wp,sp,A,B,C,D)); splitP # split plot design

# 8.4.a RCB in whole plots RBSP
# 1 blocking variable (2), 1 whole plot factor (2 x 2), 1 subplot factor (2)
library(FrF2)
spdes <- FrF2(64, 5, # 64 runs, 5 factors A, B, S1, S2, S3, all 2 levels
              WPs = 8, # 8 whole plots, 2 manufacturer blocks each with 4 wplots 
              nfac.WP = 2, # number of whole plot factors (2)  
              factor.names = (c("A", "B", "S1", "S2", "S3")))
spdes <- spdes[,c(3,1,2,4,5,6)]; spdes # reorder columns WP3 (block) first 
spdes <- data.frame(spdes[order(spdes$WP3,spdes$A,spdes$B,spdes$S1,spdes$S2,spdes$S3),]); spdes
library(plyr) # rename "WP3" as "block"
spdes <- rename(spdes, c("WP3"="block")); spdes
wplots <- factor(rep(c(1,2,3,4,1,2,3,4), each = 8)); wplots
spdes1 <- cbind(spdes, wplots); spdes1 # add wplots
spdes2 <- spdes1[,c(1,7,2,3,4,5,6)]; spdes2 # reorder columns wplots second 

# 8.6.a,b,c RCB in whole plots RBSP
# 288 runs total
R <- factor(c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), 
              rep(5,4), rep(6,4), rep(7,4), rep(8,4))); R
Run <- factor(rep(R,9)); Run
FP <- factor(rep(c(1,2,3,4), times = 72)); FP
site <- factor(rep(c(1,2,3,4,5,6,7,8,9), each = 32)); site
# wafer site subplot factor, 9 levels
thick <- c(90.1,	91.9,	88.1,	90,	90.7,	88.6,	90.2,	90.8,	89.4,	89.7,	86.6,	93.2,	87.8,	86.6,	91.9,	89.1,	91.8,	89.3,	90,	90.2,	90.3,	91.1,	92.4,	94.1,	90.3,	92.7,	87,	91.8,	89,	89.9,	89,	88.9,
           90.1,	93.3,	90.8,	93.1,	90.8,	89.1,	90.4,	92.6,	90,	90.1,	94.9,	93.9,	93.2,	92.4,	93.5,	92.1,	90.4,	94.5,	92,	90.4,	91.1,	89.8,	91.7,	91.5,	91.2,	89.3,	94,	91.8,	89.8,	90.6,	89.8,	91.6,
           92.8,	94.1,	91.5,	92.7,	90.3,	91.5,	90.9,	92.6,	93,	92.1,	91,	91.7,	91.7,	90.9,	97.9,	94.6,	91.7,	94.6,	95,	93.4,	93.3,	91.5,	91.6,	95.3,	93,	90.9,	95.8,	91.6,	89,	90.4,	89,	90.7,
           87.8,	89.1,	88.2,	91.6,	92.7,	89.5,	94.7,	88.4,	90.4,	88.6,	89,	90.3,	85.6,	90.9,	90.1,	92,	91.8,	95.8,	92.7,	92.4,	93.5,	91.5,	91.1,	92.8,	89.7,	90.2,	91.7,	94.7,	90.5,	91.8,	90.5,	91.4,
           88.2,	91.4,	90.5,	89.2,	88.4,	86.6,	91.3,	92.4,	90.4,	90,	90.9,	90.5,	90.3,	91.4,	87.7,	89.6,	89,	93,	88.5,	88.8,	87.2,	90.6,	88,	93.4,	88.1,	88.8,	89.7,	92.7,	90.1,	88.3,	90.1,	88.7,
           88.2,	92.2,	92.3,	92.6,	89,	93.4,	91.3,	89.9,	89.9,	92.6,	92.3,	93,	87.9,	90.4,	92.1,	92.4,	90,	91.7,	91.3,	91.7,	88.1,	93.1,	92.4,	92.2,	91,	92.5,	88.7,	92.5,	88.6,	93.1,	87.6,	93.1,
           90.4,	87.5,	87.4,	87,	89.1,	89.9,	90,	89.9,	91.6,	89.2,	90.5,	89.7,	89.1,	89.7,	89,	92.9,	88.9,	89.2,	90,	89.4,	90.1,	88.9,	88.7,	89.4,	89.7,	89.9,	90.7,	90.1,	90.5,	88.4,	90.5,	88.6,
           92.1,	91.2,	92.6,	93.2,	92.6,	91.8,	91.6,	91.9,	92.6,	92.5,	93.6,	92.5,	93.2,	92.6,	92,	96.2,	93.8,	93.3,	92.1,	96.7,	91.9,	92.5,	92.9,	94.5,	95,	94.2,	94.9,	94.9,	91.3,	92.1,	90.3,	92.10,
           91.8	,91.1,	92.4,	95.2,	92.8,	92.3,	92,	94.1,	93,	93,	93.6,	94.6,	90.9,	92.7,	93.4,	96.1,	92.3,	95.2,	93.9,	92.5,	94.5,	92.4,	92.6,	95.4,	95.4,	93.6,	91.4,	92.8,	93.3,	93.1,	94.3,	93.1)
oxide.data <- data.frame(Run,FP,site,thick); oxide.data
library(lme4)
spmodel <- lmer( thick ~ (1|Run) + FP + (1|Run:FP) + site + FP:site, data = oxide.data)
summary(spmodel) # summary random effects
anova(spmodel) # summary fixed effects
spmodel <- aov(thick ~ Run*FP + site + FP:site + Error(Run:FP/site), data = oxide.data ) 
summary(spmodel, type = 3) # error degrees of freedom for FP and site factors
1 - pf(7.0218,3,21) # p-value for whole plot FP factor
1 - pf(24.0271,8,224) # p-value for subplot site factor



par( mfrow = c(1,1) )
with(oxide.data,
     (interaction.plot(
       FP,
       site,
       thick,
       type = "b",
       pch = c(18,22,6,19,1,2,3,4,5),
       leg.bty = "o",
       main = "FP by site Interaction Plot",
       xlab = "FP", 
       ylab = "Quality"))
)
' this is my answer'


spmodel2 <- lmer( thick ~ (1|Run) + FP + (1|Run:FP) + site + FP:site, data = oxide.data)
res <- resid(spmodel2)
plot(fitted(spmodel2), res,
     ylab="Residual", 
     xlab="Predicted Value", 
     main="Residuals by Predicted for Y")
abline(0,0)


qqnorm(res, 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="QQ Plot of Residuals for Y")
abline(0,3/2)
